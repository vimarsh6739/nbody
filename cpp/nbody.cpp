#include "ctimer.h"
#include "octree.h"
#include "philox_engine.h"
#include <assert.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>

#define SOFTENING 1e-9f
#define POSMAX 100

// Randomize the positions and velocities in the range [-1, 1],
// and assign a mass (every 7th float) in the range [0.1, 1.0].
int randomizeBodies(PhiloxEngine &rng, Body *bodies, int n) {
  int maxKeyLength = 0;
  std::uniform_real_distribution<float> x_dis(0, 1);
  std::uniform_real_distribution<float> v_dis(0, 1);
  for (int i = 0; i < n; i++) {
    Body &body = bodies[i];
    body.x = x_dis(rng);
    body.y = x_dis(rng);
    body.z = x_dis(rng);

    body.key = getKeyNoPrepend(body);
    int keylength = binaryLength(body.key);
    if (keylength > maxKeyLength) {
      maxKeyLength = keylength;
    }

    body.vx = (v_dis(rng)) * 2 - 1;
    body.vy = (v_dis(rng)) * 2 - 1;
    body.vz = (v_dis(rng)) * 2 - 1;
    body.m = (v_dis(rng)) + 0.1f;
  }

  // prepend all keys
  int prepend = 1 << maxKeyLength;
  for (int i = 0; i < n; i++) {
    bodies[i].key += prepend;
  }

  return maxKeyLength + 1;
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
  bodyForceRange(p, dt, 0, mid, n);
  bodyForceRange(p, dt, mid, n, n);
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
  integratePositionsRange(p, dt, 0, mid);
  integratePositionsRange(p, dt, mid, n);
}

int main(const int argc, const char **argv) {
  printf("NBody using Octree\n");
  PhiloxEngine rng(2025);
  int nBodies = 2;
  if (argc > 1)
    nBodies = atoi(argv[1]);

  const float dt = 0.01f; // time step
  const int nIters = 10;  // simulation iterations

  printf("Simulating %d bodies\n", nBodies);
  // Allocate memory for nBodies
  Body *p = new Body[nBodies];

  // Initialize positions, velocities, and masses (7 values per body)
  int maxkeylength = randomizeBodies(rng, p, nBodies);
  printf("Randomized bodies with maxkeylength = %d\n", maxkeylength);

  Octree *octree = new Octree(maxkeylength);
  for (int i = 0; i < nBodies; i++) {
    octree->insert(p[i]);
  }

  octree->printTree(octree->root, 0);

  std::vector<DFTNode> dft;
  octree->buildDFT(dft);

  //   double totalTime = 0.0;
  //   ctimer_t timer;

  //   for (int iter = 1; iter <= nIters; iter++) {
  //     ctimer_start(&timer);

  //     // Compute interbody forces using cilk_scope with spawn.
  //     bodyForce(p, dt, nBodies);

  //     // Integrate positions concurrently by dividing the work.
  //     integratePositions(p, dt, nBodies);

  //     ctimer_stop(&timer);
  //     ctimer_measure(&timer);
  //     double tElapsed = timespec_sec(timer.elapsed);

  //     if (iter > 1) // skip the first iteration as warm-up
  //       totalTime += tElapsed;

  // #ifndef SHMOO
  //     printf("Iteration %d: %.3f seconds\n", iter, tElapsed);
  // #endif
  //   }

  //   double avgTime = totalTime / (double)(nIters - 1);
  // #ifdef SHMOO
  //   printf("%d, %0.3f\n", nBodies, 1e-9 * nBodies * nBodies / avgTime);
  // #else
  //   printf("Average time for iterations 2 through %d: %.3f seconds.\n",
  //   nIters,
  //          avgTime);
  //   printf("%d Bodies: average %0.3f Billion Interactions / second\n",
  //   nBodies,
  //          1e-9 * nBodies * nBodies / avgTime);
  // #endif

  delete[] p;
  return 0;
}
