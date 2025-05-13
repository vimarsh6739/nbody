

// filename: example.cu

#include "main.h"
#include "octree.h"
#include <cuda_runtime.h>
#include <iostream>

// same as CPU mac
__device__ bool dMAC(float target_x, float target_y, float target_z, float x,
                     float y, float z, float mass, float MAC_PARAM) {

  float distance =
      sqrtf((target_x - x) * (target_x - x) + (target_y - y) * (target_y - y) +
            (target_z - z) * (target_z - z));
  float source_mass = mass;

  bool result = distance > MAC_PARAM;

  return distance / source_mass > MAC_PARAM;
}

// CUDA kernel
__global__ void MACKernel(Body *bodies, Body *newBodies, float dt,
                          int n /*nbodies*/, float SOFTENING, int dft_size,
                          DFTNode *dft, float MAC_PARAM) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {
    // iterate over all bodies (targets)

    float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;
    for (int j = 0; j < dft_size; j++) {
      if (dft[j].isLeaf) {
        for (int b = 0; b < dft[j].nBodies; b++) {

          // todo: need to fix vectors in DFTNode...
          int bIndex = dft[j].bodies[b];
          float dx = bodies[bIndex].x - bodies[idx].x;
          float dy = bodies[bIndex].y - bodies[idx].y;
          float dz = bodies[bIndex].z - bodies[idx].z;
          float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
          float invDist = 1.0f / sqrtf(distSqr);
          float invDist3 = invDist * invDist * invDist;

          Fx += dx * bodies[bIndex].m * invDist3;
          Fy += dy * bodies[bIndex].m * invDist3;
          Fz += dz * bodies[bIndex].m * invDist3;
        }
      } else if (dMAC(bodies[idx].x, bodies[idx].y, bodies[idx].z, dft[j].x,
                      dft[j].y, dft[j].z, dft[j].mass, MAC_PARAM)) {

        float dx = dft[j].x - bodies[idx].x;
        float dy = dft[j].y - bodies[idx].y;
        float dz = dft[j].z - bodies[idx].z;
        float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
        float invDist = 1.0f / sqrtf(distSqr);
        float invDist3 = invDist * invDist * invDist;

        Fx += dx * dft[j].mass * invDist3;
        Fy += dy * dft[j].mass * invDist3;
        Fz += dz * dft[j].mass * invDist3;
        j = dft[j].autorope - 1; // -1 because we increment j in the for loop
      }
    }

    newBodies[idx].vx += dt * Fx;
    newBodies[idx].vy += dt * Fy;
    newBodies[idx].vz += dt * Fz;
  }
}

// CUDA kernel
__global__ void DSKernel(Body *bodies, float dt,
                         int n /*nbodies*/, float SOFTENING) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {

    double Fx = 0.0;
    double Fy = 0.0;
    double Fz = 0.0;

    // compute force, update velocities
    for (int j = 0; j < n; ++j) {
      double dx = bodies[j].x - bodies[idx].x;
      double dy = bodies[j].y - bodies[idx].y;
      double dz = bodies[j].z - bodies[idx].z;
      double distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
      double invDist = 1.0f / sqrtf(distSqr);
      double invDist3 = invDist * invDist * invDist;

      Fx += dx * bodies[j].m * invDist3;
      Fy += dy * bodies[j].m * invDist3;
      Fz += dz * bodies[j].m * invDist3;
    }

    // update velocity
    bodies[idx].vx += (double)dt * Fx;
    bodies[idx].vy += (double)dt * Fy;
    bodies[idx].vz += (double)dt * Fz;
  }
}

// CUDA kernel
__global__ void DFTKernel(Body *bodies, DFTNode *dft, Body *newBodies, float dt,
                          int n /*nbodies*/, int dft_size, float SOFTENING) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {

    float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;

    for (int j = 0; j < dft_size; j++) {
      if (dft[j].isLeaf) {
        for (int bIndex : dft[j].bodies) {
          float dx = bodies[bIndex].x - bodies[idx].x;
          float dy = bodies[bIndex].y - bodies[idx].y;
          float dz = bodies[bIndex].z - bodies[idx].z;
          float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
          float invDist = 1.0f / sqrtf(distSqr);
          float invDist3 = invDist * invDist * invDist;

          Fx += dx * bodies[bIndex].m * invDist3;
          Fy += dy * bodies[bIndex].m * invDist3;
          Fz += dz * bodies[bIndex].m * invDist3;
        }
      }
    }
    newBodies[idx].vx += dt * Fx;
    newBodies[idx].vy += dt * Fy;
    newBodies[idx].vz += dt * Fz;
  }
}

// Utility function to check for CUDA errors
void checkCuda(cudaError_t result, const char *msg) {
  if (result != cudaSuccess) {
    std::cerr << "CUDA error: " << msg << " - " << cudaGetErrorString(result)
              << std::endl;
    exit(EXIT_FAILURE);
  }
}

int launchDS(Body *bodies, float dt, int N /*nbodies*/) {
  size_t size = N * sizeof(Body);

  // Device vectors
  Body *d_bodies;

  checkCuda(cudaMalloc(&d_bodies, size), "Allocating d_bodies");
  checkCuda(cudaMemcpy(d_bodies, bodies, size, cudaMemcpyHostToDevice),
            "Copying bodies to d_bodies");


  // Launch kernel (1D grid of 256 threads per block)
  int threadsPerBlock = 256;
  int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
  DSKernel<<<blocksPerGrid, threadsPerBlock>>>(d_bodies, dt, N,
                                               SOFTENING);

  // Check for kernel launch errors
  checkCuda(cudaGetLastError(), "Launching addKernel");

  cudaDeviceSynchronize();
  // Copy result back to host
  checkCuda(cudaMemcpy(bodies, d_newBodies, size, cudaMemcpyDeviceToHost),
            "Copying result back to host");

  // Cleanup
  cudaFree(d_bodies);
  cudaFree(d_newBodies);

  return 0;
}

int launchDFT(Body *bodies, float dt, int N /*nbodies*/,
              std::vector<DFTNode> dft) {
  size_t size = N * sizeof(Body);
  int ndft = dft.size();
  size_t dft_size = dft.size() * sizeof(DFTNode);

  // Device vectors
  Body *d_bodies;
  Body *d_newBodies;
  DFTNode *d_dft;

  checkCuda(cudaMalloc(&d_bodies, size), "Allocating d_bodies");
  checkCuda(cudaMalloc(&d_newBodies, size), "Allocating d_newBodies");
  checkCuda(cudaMalloc(&d_dft, dft_size), "Allocating d_dft");

  checkCuda(cudaMemcpy(d_bodies, bodies, size, cudaMemcpyHostToDevice),
            "Copying bodies to d_bodies");
  checkCuda(cudaMemcpy(d_newBodies, d_bodies, size, cudaMemcpyDeviceToDevice),
            "D2D copy of bodies");

  checkCuda(cudaMemcpy(d_dft, dft.data(), dft_size, cudaMemcpyHostToDevice),
            "Copying dft to d_dft");

  // Launch kernel (1D grid of 256 threads per block)
  int threadsPerBlock = 256;
  int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
  DFTKernel<<<blocksPerGrid, threadsPerBlock>>>(d_bodies, d_dft, d_newBodies,
                                                dt, N, ndft, SOFTENING);

  // Check for kernel launch errors
  checkCuda(cudaGetLastError(), "Launching addKernel");

  cudaDeviceSynchronize();
  // Copy result back to host
  checkCuda(cudaMemcpy(bodies, d_newBodies, size, cudaMemcpyDeviceToHost),
            "Copying result back to host");

  // Cleanup
  cudaFree(d_bodies);
  cudaFree(d_newBodies);
  cudaFree(d_dft);

  return 0;
}

int launchMAC(Body *bodies, float dt, int N /*nbodies*/,
              std::vector<DFTNode> dft, float MAC_PARAM) {
  size_t size = N * sizeof(Body);
  int ndft = dft.size();
  size_t dft_size = dft.size() * sizeof(DFTNode);

  // Device vectors
  Body *d_bodies;
  Body *d_newBodies;
  DFTNode *d_dft;

  checkCuda(cudaMalloc(&d_bodies, size), "Allocating d_bodies");
  checkCuda(cudaMalloc(&d_newBodies, size), "Allocating d_newBodies");
  checkCuda(cudaMalloc(&d_dft, dft_size), "Allocating d_dft");

  checkCuda(cudaMemcpy(d_bodies, bodies, size, cudaMemcpyHostToDevice),
            "Copying bodies to d_bodies");
  checkCuda(cudaMemcpy(d_newBodies, d_bodies, size, cudaMemcpyDeviceToDevice),
            "D2D copy of bodies");

  checkCuda(cudaMemcpy(d_dft, dft.data(), dft_size, cudaMemcpyHostToDevice),
            "Copying dft to d_dft");

  // Launch kernel (1D grid of 256 threads per block)
  int threadsPerBlock = 256;
  int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
  MACKernel<<<blocksPerGrid, threadsPerBlock>>>(
      d_bodies, d_newBodies, dt, N, SOFTENING, ndft, d_dft, MAC_PARAM);

  // Check for kernel launch errors
  checkCuda(cudaGetLastError(), "Launching addKernel");

  cudaDeviceSynchronize();
  // Copy result back to host
  checkCuda(cudaMemcpy(bodies, d_newBodies, size, cudaMemcpyDeviceToHost),
            "Copying result back to host");

  // Cleanup
  cudaFree(d_bodies);
  cudaFree(d_newBodies);
  cudaFree(d_dft);

  return 0;
}
