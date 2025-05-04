#include "ctimer.h"
#include "octree.h"
#include "philox_engine.h"
#include <assert.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
// #include <cilk.h>

#define SOFTENING 1e-9f
#define POSMAX 100

bool USE_TREE = false; // use tree or list
bool USE_BH = true;    // use barnes hut or not
bool PRINT_TIME = true;
float MAC_PARAM = .5; // MAC parameter

bool MAC(float target_x, float target_y, float target_z, float x, float y,
         float z, float mass) {

  float distance =
      sqrtf((target_x - x) * (target_x - x) + (target_y - y) * (target_y - y) +
            (target_z - z) * (target_z - z));
  float source_mass = mass;

  bool result = distance > MAC_PARAM;

  return distance / source_mass > MAC_PARAM;
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

void allInteractionsDS(Body *bodies, float dt, int nBodies) {
  int n = nBodies;

  for (int i = 0; i < n; ++i) {

    float Fx = 0.0;
    float Fy = 0.0;
    float Fz = 0.0;

    // compute force, update velocities
    for (int j = 0; j < n; ++j) {

      float dx = bodies[j].x - bodies[i].x;
      float dy = bodies[j].y - bodies[i].y;
      float dz = bodies[j].z - bodies[i].z;
      float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
      float invDist = 1.0f / sqrtf(distSqr);
      float invDist3 = invDist * invDist * invDist;

      Fx += dx * bodies[j].m * invDist3;
      Fy += dy * bodies[j].m * invDist3;
      Fz += dz * bodies[j].m * invDist3;
    }

    // update velocity
    bodies[i].vx += dt * Fx;
    bodies[i].vy += dt * Fy;
    bodies[i].vz += dt * Fz;
  }
}

void bodyForce(Body *p, float dt, int n, std::vector<DFTNode> dft) {

  if (USE_TREE && USE_BH) {
  } else if (USE_TREE) {
  } else {
    allInteractionsDS(p, dt, n);
  }
  // // printf("DFT size: %ld\n", dft.size());
  // for (int i = 0; i < n; i++) {
  //   if (USE_TREE && USE_BH)
  //     MACInteractionsDFT(dft, p, i, dt, n);
  //   else if (USE_TREE)
  //     allInteractionsDFT(dft, p, i, dt, n);
  //   else
  //     allInteractionsDS(p, i, dt, n);
  // }
  //
  // printf("Total interactions: %ld\n", totalInteractions);
}

void integratePositions(Body *p, float dt, int start, int end) {
  for (int i = start; i < end; i++) {
    p[i].x += (p[i].vx * dt);
    p[i].y += (p[i].vy * dt);
    p[i].z += (p[i].vz * dt);
  }
}

void nbodyIterate(Body *p, float dt, int nBodies, int nIters,
                  int maxkeylength) {
  double totalTime = 0.0;
  ctimer_t timer;
  for (int iter = 1; iter <= nIters; iter++) {

    std::vector<DFTNode> dft;

    Octree *octree;

    ctimer_start(&timer);

    // Reconstruct octree in every timestep
    if (USE_TREE) {
      octree = new Octree(maxkeylength);
      int nUniqueLeaves = 0;
      for (int i = 0; i < nBodies; i++) {
        nUniqueLeaves += octree->insert(p[i]);
      }

      octree->buildDFT(dft, p);
      printf("Octree built with %d buckets (leaves with unique key)\n",
             nUniqueLeaves);
      // octree->printTree(octree->root, 0);
    }

    bodyForce(p, dt, nBodies, dft);

    if (USE_TREE)
      delete octree;

    // Integrate positions
    integratePositions(p, dt, 0, nBodies);

    ctimer_stop(&timer);
    ctimer_measure(&timer);

    double tElapsed = timespec_sec(timer.elapsed);

    // skip the first iteration as warm-up
    if (iter > 1)
      totalTime += tElapsed;
  }

  double avgTime = totalTime / (double)(nIters - 1);

  if (PRINT_TIME) {
    printf("Total time = %.3f seconds, Average time = %.3f seconds.\n",
           totalTime, avgTime);
    printf("%d Bodies: average %0.3f Billion Interactions / second\n", nBodies,
           1e-9 * nBodies * nBodies / avgTime);
  }
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
  PRINT_TIME = false;
  nbodyIterate(orig, 0.01f, nBodies, 10, 0);
  PRINT_TIME = true;
  // Compute RMS error(standard check for positional accuracy)
  float rmsError = 0.0f;
  for (int i = 0; i < nBodies; i++) {

    float dx = p[i].x - orig[i].x;
    float dy = p[i].y - orig[i].y;
    float dz = p[i].z - orig[i].z;
    rmsError += dx * dx + dy * dy + dz * dz;
  }
  rmsError = std::sqrt((rmsError / nBodies));
  printf("RMS error: %f\n", rmsError);
}

using namespace std;

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
      USE_BH = false;
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

  if (argc > 3) {
    MAC_PARAM = atof(argv[3]);
  }

  // TODO: other parameters we might want to control:
  // - SHIFT_DIGITS (controls how many keys/buckets are created)
  // - MAC_PARAM (controls the MAC parameter)
  // - SOFTENING (controls the softening factor)

  printf(" - %s method \n - %d bodies\n",
         USE_TREE ? (USE_BH ? "MAC" : "DFT") : "DS", nBodies);
  printf("-----------------------------------------------\n");

  const float dt = 0.01f; // time step
  const int nIters = 10;  // simulation iterations

  // Allocate memory for nBodies
  Body *particles = new Body[nBodies];

  // Initialize positions, velocities, and masses (7 values per body)
  PhiloxEngine rng(2025);
  int maxkeylength = randomizeBodies(rng, particles, nBodies);

  // make a copy of bodies for accuracy check
  Body *particles_cp = new Body[nBodies];
  for (int i = 0; i < nBodies; i++) {
    particles_cp[i] = particles[i];
  }

  nbodyIterate(particles, dt, nBodies, nIters, maxkeylength);
  checkAccuracy(particles, particles_cp, nBodies);

  delete[] particles;
  delete[] particles_cp;
  return 0;
}
