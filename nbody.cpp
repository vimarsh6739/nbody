#include "ctimer.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
using namespace std;

#define SOFTENING 1e-9f

typedef struct {
  float x, y, z, vx, vy, vz, m; // added mass field 'm'
} Body;

// Randomize positions and velocities in the range [-1, 1].
// For the mass (every 7th float), randomize in the range [0.1, 1.0].
void randomizeBodies(float *data, int n) {
  for (int i = 0; i < n; i++) {
    // Every 7th element corresponds to the mass of a body.
    if ((i + 1) % 7 == 0) {
      data[i] = 0.1f + 0.9f * (rand() / (float)RAND_MAX);
    } else {
      data[i] = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    }
  }
}

// Compute gravitational force on each body incorporating the mass of the other bodies.
void bodyForce(Body *p, float dt, int n) {
  for (int i = 0; i < n; i++) {
    float Fx = 0.0f;
    float Fy = 0.0f;
    float Fz = 0.0f;

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

int main(const int argc, const char **argv) {
  srand(2025);
  int nBodies = 10000;
  if (argc > 1)
    nBodies = atoi(argv[1]);

  const float dt = 0.01f; // time step
  const int nIters = 10;  // simulation iterations

  int bytes = nBodies * sizeof(Body);
  float *buf = (float *)malloc(bytes);
  Body *p = (Body *)buf;

  randomizeBodies(buf, 7 * nBodies);

  double totalTime = 0.0;
  ctimer_t timer;
  
  for (int iter = 1; iter <= nIters; iter++) {
    ctimer_start(&timer);
    bodyForce(p, dt, nBodies); // compute interbody forces

    for (int i = 0; i < nBodies; i++) { // integrate position
      p[i].x += p[i].vx * dt;
      p[i].y += p[i].vy * dt;
      p[i].z += p[i].vz * dt;
    }
    
    ctimer_stop(&timer);
    ctimer_measure(&timer);
    const double tElapsed = timespec_sec(timer.elapsed); 

    if (iter > 1) { // First iter is warm up
      totalTime += tElapsed;
    }

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
