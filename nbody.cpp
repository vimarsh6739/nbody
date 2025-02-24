#include "ctimer.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <cilk/cilk.h>

#define SOFTENING 1e-9f

typedef struct {
  float x, y, z, vx, vy, vz, m;  // each body now has a mass 'm'
} Body;

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

void bodyForceRange(Body *p, float dt, int start, int end, int n) {
  for (int i = start; i < end; i++) {
    float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;
    for (int j = 0; j < n; j++) {
      float dx = p[j].x - p[i].x;
      float dy = p[j].y - p[i].y;
      float dz = p[j].z - p[i].z;
      float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
      float invDist = 1.0f / sqrtf(distSqr);
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

void bodyForce(Body *p, float dt, int n) {
  int mid = n / 2;
  cilk_scope {
    cilk_spawn bodyForceRange(p, dt, 0, mid, n);
    bodyForceRange(p, dt, mid, n, n);
  }
}

void integratePositionsRange(Body *p, float dt, int start, int end) {
  for (int i = start; i < end; i++) {
    p[i].x += p[i].vx * dt;
    p[i].y += p[i].vy * dt;
    p[i].z += p[i].vz * dt;
  }
}

void integratePositions(Body *p, float dt, int n) {
  int mid = n / 2;
  cilk_scope {
    cilk_spawn integratePositionsRange(p, dt, 0, mid);
    integratePositionsRange(p, dt, mid, n);
  }
}

int main(const int argc, const char **argv) {
  srand(2025);
  int nBodies = 10000;
  if (argc > 1)
    nBodies = atoi(argv[1]);

  const float dt = 0.01f; // time step
  const int nIters = 10;  // simulation iterations

  // Allocate memory for nBodies (each Body has 7 floats)
  int bytes = nBodies * sizeof(Body);
  float *buf = (float *)malloc(bytes);
  Body *p = (Body *)buf;

  // Initialize positions, velocities, and masses (7 values per body)
  randomizeBodies(buf, 7 * nBodies);

  double totalTime = 0.0;
  ctimer_t timer;

  for (int iter = 1; iter <= nIters; iter++) {
    ctimer_start(&timer);

    // Compute interbody forces using cilk_scope with spawn.
    bodyForce(p, dt, nBodies);

    // Integrate positions concurrently by dividing the work.
    integratePositions(p, dt, nBodies);

    ctimer_stop(&timer);
    ctimer_measure(&timer);
    double tElapsed = timespec_sec(timer.elapsed);

    if (iter > 1) // skip the first iteration as warm-up
      totalTime += tElapsed;

#ifndef SHMOO
    printf("Iteration %d: %.3f seconds\n", iter, tElapsed);
#endif
  }

  double avgTime = totalTime / (double)(nIters - 1);
#ifdef SHMOO
  printf("%d, %0.3f\n", nBodies, 1e-9 * nBodies * nBodies / avgTime);
#else
  printf("Average time for iterations 2 through %d: %.3f seconds.\n", nIters, avgTime);
  printf("%d Bodies: average %0.3f Billion Interactions / second\n", nBodies,
         1e-9 * nBodies * nBodies / avgTime);
#endif

  free(buf);
  return 0;
}
