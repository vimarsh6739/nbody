#include "ctimer.h"
#include "philox_engine.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>
#include <random>

#define SOFTENING 1e-9f
#define BLOCK_SIZE 256

typedef struct {
  float x, y, z, vx, vy, vz, m;
} Body;

#define CUDA_CALL(call)                                                  \
  do {                                                                         \
    cudaError_t err = call;                                                    \
    if (err != cudaSuccess) {                                                  \
      fprintf(stderr, "CUDA error in %s:%d: %s\n", __FILE__, __LINE__,         \
              cudaGetErrorString(err));                                        \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  } while (0)

void randomizeBodies(PhiloxEngine &rng, Body *bodies, int n) {
  std::uniform_real_distribution<float> x_dis(0, 1);
  std::uniform_real_distribution<float> v_dis(0, 1);

  for (int i = 0; i < n; i++) {
    Body &body = bodies[i];
    body.x = 1 * x_dis(rng);
    body.y = 1 * x_dis(rng);
    body.z = 1 * x_dis(rng);

    body.vx = (v_dis(rng)) * 2 - 1;
    body.vy = (v_dis(rng)) * 2 - 1;
    body.vz = (v_dis(rng)) * 2 - 1;
    body.m = (v_dis(rng)) + 0.15f;
  }
}

__global__ void bodyForceKernel(Body *body, float dt, int n) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n)
    return;
  float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;

  for (int j = 0; j < n; j++) {
    float dx = body[j].x - body[i].x;
    float dy = body[j].y - body[i].y;
    float dz = body[j].z - body[i].z;

    float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
    float invDist = rsqrtf(distSqr); // GPU optimized inverse square root
    float invDist3 = invDist * invDist * invDist;

    // Compute force components
    Fx += dx * body[j].m * invDist3;
    Fy += dy * body[j].m * invDist3;
    Fz += dz * body[j].m * invDist3;
  }

  // update velocities
  body[i].vx += dt * Fx;
  body[i].vy += dt * Fy;
  body[i].vz += dt * Fz;
}

// CUDA kernel to integrate positions
__global__ void integratePositionsKernel(Body *p, float dt, int n) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < n) {
    p[tid].x += p[tid].vx * dt;
    p[tid].y += p[tid].vy * dt;
    p[tid].z += p[tid].vz * dt;
  }
}

int main(const int argc, const char **argv) {
  PhiloxEngine rng(2025); // set the rng with the seed
  int nBodies = 10000;
  if (argc > 1)
    nBodies = atoi(argv[1]);

  const float dt = 0.01f;
  const int nIters = 15;

  int bytes = nBodies * sizeof(Body);
  Body *h_bodies = (Body *)calloc(nBodies, sizeof(Body));

  randomizeBodies(rng, h_bodies, nBodies);

  Body *d_bodies;
  CUDA_CALL(cudaMalloc(&d_bodies, bytes));

  CUDA_CALL(
      cudaMemcpy(d_bodies, h_bodies, bytes, cudaMemcpyHostToDevice));

  const int blockSize = BLOCK_SIZE;
  const int gridSize = (nBodies + blockSize - 1) / blockSize;

  double totalTime = 0.0;
  ctimer_t timer;

  for (int iter = 1; iter <= nIters; iter++) {
    ctimer_start(&timer);

    bodyForceKernel<<<gridSize, blockSize>>>(d_bodies, dt, nBodies);
    CUDA_CALL(cudaGetLastError());

    integratePositionsKernel<<<gridSize, blockSize>>>(d_bodies, dt, nBodies);
    CUDA_CALL(cudaGetLastError());

    CUDA_CALL(cudaDeviceSynchronize());

    ctimer_stop(&timer);
    ctimer_measure(&timer);
    double tElapsed = timespec_sec(timer.elapsed);

    if (iter > 1)
      totalTime += tElapsed;

#ifndef SHMOO
    printf("Iteration %d: %.3f seconds\n", iter, tElapsed);
#endif
  }

  CUDA_CALL(
      cudaMemcpy(h_bodies, d_bodies, bytes, cudaMemcpyDeviceToHost));

  double avgTime = totalTime / (double)(nIters - 1);
#ifdef SHMOO
  printf("%d, %0.3f\n", nBodies, 1e-9 * nBodies * nBodies / avgTime);
#else
  printf("Average time for iterations 2 through %d: %.3f seconds.\n", nIters,
         avgTime);
  printf("%d Bodies: average %0.3f Billion Interactions / second\n", nBodies,
         1e-9 * nBodies * nBodies / avgTime);
#endif

  free(h_bodies);
  CUDA_CALL(cudaFree(d_bodies));

  return 0;
}
