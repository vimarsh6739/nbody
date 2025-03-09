#include "ctimer.h"
#include "philox_engine.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>

#define SOFTENING 1e-9f
#define BLOCK_SIZE 256 

typedef struct {
  float x, y, z, vx, vy, vz, m; 
} Body;

#define checkCudaErrors(call) do {                                 \
  cudaError_t err = call;                                          \
  if (err != cudaSuccess) {                                        \
    fprintf(stderr, "CUDA error in %s:%d: %s\n",                   \
            __FILE__, __LINE__, cudaGetErrorString(err));          \
    exit(EXIT_FAILURE);                                            \
  }                                                                \
} while (0)

void randomizeBodies(float *data, int n) {
  for (int i = 0; i < n; i++) {
    if ((i + 1) % 7 == 0) {
      data[i] = 0.1f + 0.9f * philox_random_float();
    } else {
      data[i] = 2.0f * philox_random_float() - 1.0f;
    }
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
  philox_seed(2025);
  int nBodies = 10000;
  if (argc > 1)
    nBodies = atoi(argv[1]);

  const float dt = 0.01f;
  const int nIters = 10; 

  int bytes = nBodies * sizeof(Body);
  Body *h_bodies = (Body *)malloc(bytes);

  randomizeBodies((float *)h_bodies, 7 * nBodies);

  Body *d_bodies;
  checkCudaErrors(cudaMalloc(&d_bodies, bytes));

  checkCudaErrors(cudaMemcpy(d_bodies, h_bodies, bytes, cudaMemcpyHostToDevice));

  const int blockSize = BLOCK_SIZE;
  const int gridSize = (nBodies + blockSize - 1) / blockSize;

  double totalTime = 0.0;
  ctimer_t timer;

  for (int iter = 1; iter <= nIters; iter++) {
    ctimer_start(&timer);

    bodyForceKernel<<<gridSize, blockSize>>>(d_bodies, dt, nBodies);
    checkCudaErrors(cudaGetLastError());

    integratePositionsKernel<<<gridSize, blockSize>>>(d_bodies, dt, nBodies);
    checkCudaErrors(cudaGetLastError());

    checkCudaErrors(cudaDeviceSynchronize());

    ctimer_stop(&timer);
    ctimer_measure(&timer);
    double tElapsed = timespec_sec(timer.elapsed);

    if (iter > 1)
      totalTime += tElapsed;

#ifndef SHMOO
    printf("Iteration %d: %.3f seconds\n", iter, tElapsed);
#endif
  }

  checkCudaErrors(cudaMemcpy(h_bodies, d_bodies, bytes, cudaMemcpyDeviceToHost));

  double avgTime = totalTime / (double)(nIters - 1);
#ifdef SHMOO
  printf("%d, %0.3f\n", nBodies, 1e-9 * nBodies * nBodies / avgTime);
#else
  printf("Average time for iterations 2 through %d: %.3f seconds.\n", nIters, avgTime);
  printf("%d Bodies: average %0.3f Billion Interactions / second\n", nBodies,
         1e-9 * nBodies * nBodies / avgTime);
#endif

  free(h_bodies);
  checkCudaErrors(cudaFree(d_bodies));

  return 0;
}