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

bool USE_TREE = true; // use tree or list
bool USE_BH = true;   // use barnes hut or not
float MAC_PARAM = 1;  // MAC parameter

bool MAC(float target_x, float target_y, float target_z, float x, float y,
         float z, float mass) {

  float distance =
      sqrtf((target_x - x) * (target_x - x) + (target_y - y) * (target_y - y) +
            (target_z - z) * (target_z - z));
  float source_mass = mass;

  bool result = distance > MAC_PARAM;

  return distance > MAC_PARAM;
}

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

    body.index = i;

    body.key = getKeyNoPrepend(body);
    int keylength = binaryLength(body.key);
    if (keylength > maxKeyLength) {
      maxKeyLength = keylength;
    }

    body.vx = (v_dis(rng)) * 2 - 1;
    body.vy = (v_dis(rng)) * 2 - 1;
    body.vz = (v_dis(rng)) * 2 - 1;
    body.m = (v_dis(rng)) + 0.1f;

    FILE *file = fopen("bodiesBefore.txt", "a");
    if (file) {
      fprintf(file, "Body %d: x=%f, y=%f, z=%f\n", i, body.x, body.y, body.z);
      fclose(file);
    } else {
      printf("Error opening file for writing.\n");
    }
  }

  // prepend all keys
  int prepend = 1 << maxKeyLength;
  for (int i = 0; i < n; i++) {
    bodies[i].key += prepend;
  }

  return maxKeyLength + 1;
}

int MACInteractionsDFT(std::vector<DFTNode> dft, Body *bodies, int target,
                       float dt, int nBodies) {
  // iterate over all bodies (targets)
  int nInteractions = 0;
  float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;
  for (int j = 0; j < dft.size(); j++) {
    if (dft[j].isLeaf) {
      for (int bIndex : dft[j].bodies) {
        nInteractions++;
        float dx = bodies[bIndex].x - bodies[target].x;
        float dy = bodies[bIndex].y - bodies[target].y;
        float dz = bodies[bIndex].z - bodies[target].z;
        float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
        float invDist = 1.0f / sqrtf(distSqr);
        float invDist3 = invDist * invDist * invDist;

        Fx += dx * bodies[bIndex].m * invDist3;
        Fy += dy * bodies[bIndex].m * invDist3;
        Fz += dz * bodies[bIndex].m * invDist3;
      }
    } else if (MAC(bodies[target].x, bodies[target].y, bodies[target].z,
                   dft[j].x, dft[j].y, dft[j].z, dft[j].mass)) {
      nInteractions++;
      float dx = dft[j].x - bodies[target].x;
      float dy = dft[j].y - bodies[target].y;
      float dz = dft[j].z - bodies[target].z;
      float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
      float invDist = 1.0f / sqrtf(distSqr);
      float invDist3 = invDist * invDist * invDist;

      Fx += dx * dft[j].mass * invDist3;
      Fy += dy * dft[j].mass * invDist3;
      Fz += dz * dft[j].mass * invDist3;
      j = dft[j].autorope - 1; // -1 because we increment j in the for loop
    }
  }

  bodies[target].vx += dt * Fx;
  bodies[target].vy += dt * Fy;
  bodies[target].vz += dt * Fz;
  return nInteractions;
}

int allInteractionsDFT(std::vector<DFTNode> dft, Body *bodies, int target,
                       float dt, int nBodies) {
  // iterate over all bodies (targets)
  int nInteractions = 0;
  float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;
  for (int j = 0; j < dft.size(); j++) {
    if (dft[j].isLeaf) {
      for (int bIndex : dft[j].bodies) {
        nInteractions++;
        float dx = bodies[bIndex].x - bodies[target].x;
        float dy = bodies[bIndex].y - bodies[target].y;
        float dz = bodies[bIndex].z - bodies[target].z;
        float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
        float invDist = 1.0f / sqrtf(distSqr);
        float invDist3 = invDist * invDist * invDist;

        Fx += dx * bodies[bIndex].m * invDist3;
        Fy += dy * bodies[bIndex].m * invDist3;
        Fz += dz * bodies[bIndex].m * invDist3;
      }
    }
  }
  bodies[target].vx += dt * Fx;
  bodies[target].vy += dt * Fy;
  bodies[target].vz += dt * Fz;
  return nInteractions;
}

int allInteractionsDS(Body *bodies, int target, float dt, int nBodies) {
  int nInteractions = 0;

  // iterate over all bodies (targets)
  float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;
  for (int j = 0; j < nBodies; j++) {
    nInteractions++;
    float dx = bodies[j].x - bodies[target].x;
    float dy = bodies[j].y - bodies[target].y;
    float dz = bodies[j].z - bodies[target].z;
    float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
    float invDist = 1.0f / sqrtf(distSqr);
    float invDist3 = invDist * invDist * invDist;

    Fx += dx * bodies[j].m * invDist3;
    Fy += dy * bodies[j].m * invDist3;
    Fz += dz * bodies[j].m * invDist3;
  }
  bodies[target].vx += dt * Fx;
  bodies[target].vy += dt * Fy;
  bodies[target].vz += dt * Fz;
  return nInteractions;
}
void bodyForce(Body *p, float dt, int n, std::vector<DFTNode> dft) {
  int mid = n / 2;
  int totalInteractions = 0;
  for (int i = 0; i < n; i++) {
    if (USE_TREE && USE_BH)
      totalInteractions += MACInteractionsDFT(dft, p, i, dt, n);
    else if (USE_TREE)
      totalInteractions += allInteractionsDFT(dft, p, i, dt, n);
    else
      totalInteractions += allInteractionsDS(p, i, dt, n);
  }

  printf("Total interactions: %d\n", totalInteractions);
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

void nbodyIterate(Body *p, float dt, int nBodies, std::vector<DFTNode> dft,
                  int nIters) {
  double totalTime = 0.0;
  ctimer_t timer;
  for (int iter = 1; iter <= nIters; iter++) {
    ctimer_start(&timer);

    // Compute interbody forces using cilk_scope with spawn.
    bodyForce(p, dt, nBodies, dft);

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
  printf("Average time for iterations 2 through %d: %.3f seconds.\n", nIters,
         avgTime);
  printf("%d Bodies: average %0.3f Billion Interactions / second\n", nBodies,
         1e-9 * nBodies * nBodies / avgTime);
#endif
}

void printBodiesToFile(Body *p, int nBodies) {
  FILE *file;
  if (USE_TREE)
    if (USE_BH)
      file = fopen("bodiesAfterBH.txt", "a");
    else
      file = fopen("bodiesAfterDFT.txt", "a");
  else
    file = fopen("bodiesAfterDS.txt", "a");

  if (file) {
    for (int i = 0; i < nBodies; i++) {
      fprintf(file, "Body %d: x=%f, y=%f, z=%f\n", i, p[i].x, p[i].y, p[i].z);
    }
    fclose(file);
  } else {
    printf("Error opening file for writing.\n");
  }
}

void checkAccuracy(Body *p, Body *orig, int nBodies) {
  printf("-----------------------------------------------\n");
  printf("ACCURACY CHECK AGAINST DS\n");

  USE_TREE = false;
  USE_BH = false;

  nbodyIterate(orig, 0.01f, nBodies, std::vector<DFTNode>(), 1);

  float maxError = 0.0f;
  for (int i = 0; i < nBodies; i++) {

    float dx = p[i].x - orig[i].x;
    float dy = p[i].y - orig[i].y;
    float dz = p[i].z - orig[i].z;
    float error = sqrtf(dx * dx + dy * dy + dz * dz);
    if (error > maxError) {
      maxError = error;
    }
  }
  printf("Max error: %f\n", maxError);
}

int main(const int argc, const char **argv) {
  printf("NBody Simulation\n");

  int nBodies = 2000;
  if (argc > 1)
    nBodies = atoi(argv[1]);

  if (argc > 2) {
    nBodies = atoi(argv[1]);
    const char *method = argv[2];
    if (strcmp(method, "DS") == 0) {
      USE_TREE = false;
    } else if (strcmp(method, "DFT") == 0) {
      USE_TREE = true;
      USE_BH = false;
    } else if (strcmp(method, "MAC") == 0) {
      USE_TREE = true;
      USE_BH = true;
    } else {
      printf("Unknown method: %s\n", method);
      return -1;
    }
  }

  // TODO: other parameters we might want to control:
  // - SHIFT_DIGITS (controls how many keys/buckets are created)
  // - MAC_PARAM (controls the MAC parameter)
  // - SOFTENING (controls the softening factor)

  printf(" - %s method \n - %d bodies\n",
         USE_TREE ? (USE_BH ? "MAC" : "DFT") : "DS", nBodies);
  printf("-----------------------------------------------\n");

  const float dt = 0.01f; // time step
  const int nIters = 1;   // simulation iterations

  // Allocate memory for nBodies
  Body *p = new Body[nBodies];

  // Initialize positions, velocities, and masses (7 values per body)
  PhiloxEngine rng(2025);
  int maxkeylength = randomizeBodies(rng, p, nBodies);

  // make a copy of bodies for accuracy check
  Body *pCopy = new Body[nBodies];
  for (int i = 0; i < nBodies; i++) {
    pCopy[i] = p[i];
  }

  std::vector<DFTNode> dft;
  Octree *octree = new Octree(maxkeylength);

  if (USE_TREE) {
    int nUniqueLeaves = 0;
    for (int i = 0; i < nBodies; i++) {
      nUniqueLeaves += octree->insert(p[i]);
    }

    octree->buildDFT(dft, p);
    printf("Octree built with %d buckets (leaves with unique key)\n",
           nUniqueLeaves);
    // octree->printTree(octree->root, 0);
  }

  nbodyIterate(p, dt, nBodies, dft, nIters);
  checkAccuracy(p, pCopy, nBodies);

  delete[] p;
  delete[] pCopy;
  delete octree;
  return 0;
}
