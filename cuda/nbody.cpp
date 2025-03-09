#include "ctimer.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>

#define SOFTENING 1e-9f
#define BLOCK_SIZE 256  // CUDA thread block size

typedef struct {
  float x, y, z, vx, vy, vz, m;  // each body now has a mass 'm'
} Body;

// CUDA error checking helper function
#define checkCudaErrors(call) do {                                 \
  cudaError_t err = call;                                          \
  if (err != cudaSuccess) {                                        \
    fprintf(stderr, "CUDA error in %s:%d: %s\n",                   \
            __FILE__, __LINE__, cudaGetErrorString(err));          \
    exit(EXIT_FAILURE);                                            \
  }                                                                \
} while (0)

// Randomize the positions and velocities in the range [-1, 1],
// and assign a mass (every 7th float) in the range [0.1, 1.0].
void randomizeBodies(float *data, int n) {
  for (int i = 0; i < n; i++) {
    if ((i + 1) % 7 == 0)
      data[i] = 0.1f + 0.9f * (rand() / (float)RAND_MAX);
    else
      data[i] = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
  }
}

// CUDA kernel to compute body forces
__global__ void bodyForceKernel(Body *p, float dt, int n) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (i < n) {
    float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;
    
    for (int j = 0; j < n; j++) {
      float dx = p[j].x - p[i].x;
      float dy = p[j].y - p[i].y;
      float dz = p[j].z - p[i].z;
      float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
      float invDist = rsqrtf(distSqr);  // GPU optimized inverse square root
      float invDist3 = invDist * invDist * invDist;

      Fx += dx * p[j].m * invDist3;
      Fy += dy * p[j].m * invDist3;
      Fz += dz * p[j].m * invDist3;
    }
    
    p[i].vx += dt * Fx;
    p[i].vy += dt * Fy;
    p[i].vz += dt * Fz;
  }
}

// CUDA kernel to integrate positions
__global__ void integratePositionsKernel(Body *p, float dt, int n) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    p[i].x += p[i].vx * dt;
    p[i].y += p[i].vy * dt;
    p[i].z += p[i].vz * dt;
  }
}

int main(const int argc, const char **argv) {
  srand(2025);
  int nBodies = 10000;
  if (argc > 1)
    nBodies = atoi(argv[1]);

  const float dt = 0.01f; // time step
  const int nIters = 10;  // simulation iterations

  // Allocate memory for nBodies on CPU
  int bytes = nBodies * sizeof(Body);
  Body *h_bodies = (Body *)malloc(bytes);

  // Initialize positions, velocities, and masses on CPU
  randomizeBodies((float *)h_bodies, 7 * nBodies);

  // Allocate memory on GPU
  Body *d_bodies;
  checkCudaErrors(cudaMalloc(&d_bodies, bytes));

  // Copy data from CPU to GPU once at the beginning
  checkCudaErrors(cudaMemcpy(d_bodies, h_bodies, bytes, cudaMemcpyHostToDevice));

  // Define CUDA grid dimensions
  const int blockSize = BLOCK_SIZE;
  const int gridSize = (nBodies + blockSize - 1) / blockSize;

  double totalTime = 0.0;
  ctimer_t timer;

  for (int iter = 1; iter <= nIters; iter++) {
    ctimer_start(&timer);

    // Compute interbody forces
    bodyForceKernel<<<gridSize, blockSize>>>(d_bodies, dt, nBodies);
    checkCudaErrors(cudaGetLastError());

    // Integrate positions
    integratePositionsKernel<<<gridSize, blockSize>>>(d_bodies, dt, nBodies);
    checkCudaErrors(cudaGetLastError());

    // Synchronize GPU to make sure computation is done before stopping the timer
    checkCudaErrors(cudaDeviceSynchronize());

    ctimer_stop(&timer);
    ctimer_measure(&timer);
    double tElapsed = timespec_sec(timer.elapsed);

    if (iter > 1) // skip the first iteration as warm-up
      totalTime += tElapsed;

#ifndef SHMOO
    printf("Iteration %d: %.3f seconds\n", iter, tElapsed);
#endif
  }

  // Copy results back to CPU once at the end
  checkCudaErrors(cudaMemcpy(h_bodies, d_bodies, bytes, cudaMemcpyDeviceToHost));

  double avgTime = totalTime / (double)(nIters - 1);
#ifdef SHMOO
  printf("%d, %0.3f\n", nBodies, 1e-9 * nBodies * nBodies / avgTime);
#else
  printf("Average time for iterations 2 through %d: %.3f seconds.\n", nIters, avgTime);
  printf("%d Bodies: average %0.3f Billion Interactions / second\n", nBodies,
         1e-9 * nBodies * nBodies / avgTime);
#endif

  // Cleanup
  free(h_bodies);
  checkCudaErrors(cudaFree(d_bodies));

  return 0;
}