

// filename: example.cu

#include "gpu-core.h"
#include "octree.h"
#include <cuda_runtime.h>
#include <iostream>

// CUDA kernel
__global__ void DSKernel(Body *bodies, float dt, int n /*nbodies*/,
                         float SOFTENING) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {

    float Fx = 0.0;
    float Fy = 0.0;
    float Fz = 0.0;

    // compute force, update velocities
    for (int j = 0; j < n; ++j) {

      float dx = bodies[j].x - bodies[idx].x;
      float dy = bodies[j].y - bodies[idx].y;
      float dz = bodies[j].z - bodies[idx].z;
      float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
      float invDist = 1.0f / sqrtf(distSqr);
      float invDist3 = invDist * invDist * invDist;

      Fx += dx * bodies[j].m * invDist3;
      Fy += dy * bodies[j].m * invDist3;
      Fz += dz * bodies[j].m * invDist3;
    }

    // update velocity
    bodies[idx].vx += dt * Fx;
    bodies[idx].vy += dt * Fy;
    bodies[idx].vz += dt * Fz;
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
  DSKernel<<<blocksPerGrid, threadsPerBlock>>>(d_bodies, dt, N, SOFTENING);

  // Check for kernel launch errors
  checkCuda(cudaGetLastError(), "Launching addKernel");

  // Copy result back to host
  checkCuda(cudaMemcpy(bodies, d_bodies, size, cudaMemcpyDeviceToHost),
            "Copying result back to host");

  // Cleanup
  cudaFree(d_bodies);

  return 0;
}
