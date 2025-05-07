

// filename: example.cu

#include "main.h"
#include "octree.h"
#include <cuda_runtime.h>
#include <iostream>

// CUDA kernel
__global__ void DSKernel(Body *bodies, Body *newBodies, float dt,
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
    newBodies[idx].vx += (double)dt * Fx;
    newBodies[idx].vy += (double)dt * Fy;
    newBodies[idx].vz += (double)dt * Fz;
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
  Body *d_newBodies;

  checkCuda(cudaMalloc(&d_bodies, size), "Allocating d_bodies");
  checkCuda(cudaMalloc(&d_newBodies, size), "Allocating d_newBodies");
  checkCuda(cudaMemcpy(d_bodies, bodies, size, cudaMemcpyHostToDevice),
            "Copying bodies to d_bodies");
  checkCuda(cudaMemcpy(d_newBodies, d_bodies, size, cudaMemcpyDeviceToDevice), "D2D copy of bodies");

  // Launch kernel (1D grid of 256 threads per block)
  int threadsPerBlock = 256;
  int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
  DSKernel<<<blocksPerGrid, threadsPerBlock>>>(d_bodies, d_newBodies, dt, N,
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
