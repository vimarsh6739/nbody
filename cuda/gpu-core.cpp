#include "cpu-core.h"
#include "ctimer.h"
#include <assert.h>

#include <random>

#ifdef ENABLE_CUDA
bool CUDA = true;
#else
bool CUDA = false;
#endif

// Define all CLI options

extern int launchDS(Body *bodies, float dt, int nBodies);

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
void randomizeBodies(PhiloxEngine &rng, Body *bodies, int n) {
  std::uniform_real_distribution<float> x_dis(0, 1);
  std::uniform_real_distribution<float> v_dis(0, 1);

  for (int i = 0; i < n; i++) {
    Body &body = bodies[i];
    body.x = x_dis(rng);
    body.y = x_dis(rng);
    body.z = x_dis(rng);

    body.index = i;

    body.key = getKey(body, SHIFT_DIGITS);

    body.vx = (v_dis(rng)) * 2 - 1;
    body.vy = (v_dis(rng)) * 2 - 1;
    body.vz = (v_dis(rng)) * 2 - 1;
    body.m = (v_dis(rng)) + 0.1f;
  }
}

void createOctree(std::vector<DFTNode> &dft, Octree *octree, Body *bodies,
                  int nBodies) {

  int nUniqueLeaves = 0;
  for (int i = 0; i < nBodies; i++) {
    nUniqueLeaves += octree->insert(bodies[i]);
  }
  octree->buildDFT(dft, bodies);
  // printf("Octree built with %d buckets (leaves with unique key)\n",
  //        nUniqueLeaves);
  // octree->printTree(octree->root, 0);
  // octree->printTree(octree->root, 0);
  assert(octree->root->nLeaves == nBodies);
}

void MACInteractionsDFT(Body *bodies, float dt, int nBodies,
                        std::vector<DFTNode> &dft) {

  // iterate over all bodies (targets)

  for (int target = 0; target < nBodies; target++) {
    float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;
    for (int j = 0; j < dft.size(); j++) {
      if (dft[j].isLeaf) {
        for (int bIndex : dft[j].bodies) {

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
  }
}

void allInteractionsDFT(Body *bodies, float dt, int nBodies,
                        std::vector<DFTNode> &dft) {

  for (int target = 0; target < nBodies; target++) {
    float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;

    for (int j = 0; j < dft.size(); j++) {
      if (dft[j].isLeaf) {
        for (int bIndex : dft[j].bodies) {
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
  }
}

void allInteractionsDS(Body *bodies, float dt, int nBodies) {
  int n = nBodies;

  for (int i = 0; i < n; ++i) {
    double Fx = 0.0;
    double Fy = 0.0;
    double Fz = 0.0;

    // compute force, update velocities
    for (int j = 0; j < n; ++j) {

      double dx = bodies[j].x - bodies[i].x;
      double dy = bodies[j].y - bodies[i].y;
      double dz = bodies[j].z - bodies[i].z;
      double distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
      double invDist = 1.0f / sqrtf(distSqr);
      double invDist3 = invDist * invDist * invDist;

      Fx += dx * bodies[j].m * invDist3;
      Fy += dy * bodies[j].m * invDist3;
      Fz += dz * bodies[j].m * invDist3;
    }

    // update velocity
    bodies[i].vx += (double) dt * Fx;
    bodies[i].vy += (double) dt * Fy;
    bodies[i].vz += (double) dt * Fz;
  }
}

float checkAccuracy(Body *p, Body *orig, int nBodies, int nIters) {
  printf("-----------------------------------------------\n");
  printf("ACCURACY CHECK AGAINST DS\n");

  USE_TREE = false;
  USE_BH = false;
  PRINT_TIME = false;
  CUDA = false;
  nbodyIterate(orig, 0.01f, nBodies, nIters);
  PRINT_TIME = true;

  // Compute RMS error(standard check for positional accuracy)
  double rmsError = 0.0f;
  for (int i = 0; i < nBodies; i++) {
    double dx = p[i].x - orig[i].x;
    double dy = p[i].y - orig[i].y;
    double dz = p[i].z - orig[i].z;
    rmsError += dx * dx + dy * dy + dz * dz;
  }

  rmsError = std::sqrt((rmsError / nBodies));

  return rmsError;
}

void nbodyIterate(Body *p, float dt, int nBodies, int nIters) {
  double totalTime = 0.0;
  ctimer_t timer;
  std::vector<DFTNode> dft;
  for (int iter = 1; iter <= nIters; iter++) {
    ctimer_start(&timer);

    if (USE_TREE && iter == 1) {
      Octree *octree = new Octree(MAX_KEY_LENGTH);
      createOctree(dft, octree, p, nBodies);
    }

    int nInteractions = bodyForce(p, dt, nBodies, dft);

    // Integrate positions
    integratePositions(p, dt, 0, nBodies);

    ctimer_stop(&timer);
    ctimer_measure(&timer);

    double tElapsed = timespec_sec(timer.elapsed);

    // skip the first iteration as warm-up
    if (iter > 1)
      totalTime += tElapsed;

    printf("Iteration %d: time = %.3f seconds\n", iter, tElapsed);
  }

  double avgTime = totalTime / (double)(nIters - 1);

  if (PRINT_TIME) {
    printf("Total time = %.3f seconds, Average time = %.3f seconds.\n",
           totalTime, avgTime);
    printf("%d Bodies: average %0.3f Billion Interactions / second\n", nBodies,
           1e-9 * nBodies * nBodies / avgTime);
  }
}

int bodyForce(Body *p, float dt, int n, std::vector<DFTNode> dft) {

  if (USE_TREE && USE_BH)
    MACInteractionsDFT(p, dt, n, dft);
  else if (USE_TREE)
    allInteractionsDFT(p, dt, n, dft);
  else if (CUDA)
    launchDS(p, dt, n);
  else
    allInteractionsDS(p, dt, n);

  return 0;
}

void integratePositions(Body *p, float dt, int start, int end) {
  for (int i = start; i < end; i++) {
    p[i].x += (p[i].vx * dt);
    p[i].y += (p[i].vy * dt);
    p[i].z += (p[i].vz * dt);
  }
}
