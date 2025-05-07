#include "ctimer.h"
#include "cxxopts.hpp"
#include "octree.h"
#include "philox_engine.h"
#include <assert.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <string>
// #include <cilk.h>

// Define all CLI options
float SOFTENING = 1e-9f;
#define POSMAX 100

bool USE_TREE = false; // use tree or list
bool USE_BH = true;    // use barnes hut or not
bool PRINT_TIME = true;
float MAC_PARAM = .5; // MAC parameter

int SHIFT_DIGITS = 16;

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
  int maxKeyLength = sizeof(Key) * 8;
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
  std::cout << "the maximum key length is = " << maxKeyLength << std::endl;

  return maxKeyLength;
}

int MACInteractionsDFT(std::vector<DFTNode> dft, Body *bodies, float dt,
                       int nBodies) {
  // iterate over all bodies (targets)
  int nInteractions = 0;
  for (int target = 0; target < nBodies; target++) {
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
  }
  return nInteractions;
}

int allInteractionsDFT(std::vector<DFTNode> dft, Body *bodies, float dt,
                       int nBodies) {
  // iterate over all bodies (targets)
  int nInteractions = 0;
  for (int target = 0; target < nBodies; target++) {
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
  }
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

  int nInteractions = 0;
  if (USE_TREE && USE_BH)
    nInteractions = MACInteractionsDFT(dft, p, dt, n);
  else if (USE_TREE)
    nInteractions = allInteractionsDFT(dft, p, dt, n);
  else {
    nInteractions = n * n;
    allInteractionsDS(p, dt, n);
  }

  printf("Total interactions: %ld\n", nInteractions);
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
      // printf("Octree built with %d buckets (leaves with unique key)\n",
      //        nUniqueLeaves);
      // octree->printTree(octree->root, 0);

      printf("Octree created with %d nodes\n", octree->root->subTreeSize);
      // octree->printTree(octree->root, 0);
      assert(octree->root->nLeaves == nBodies);
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

void checkAccuracy(Body *p, Body *orig, int nBodies, int nIters) {
  printf("-----------------------------------------------\n");
  printf("ACCURACY CHECK AGAINST DS\n");

  USE_TREE = false;
  USE_BH = false;
  PRINT_TIME = false;
  nbodyIterate(orig, 0.01f, nBodies, nIters, 0);
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

int main(const int argc, const char **argv) {

  cxxopts::Options options(argv[0], "Run n-body simulations on the CPU");
  options.add_options()("n,nbodies", "Number of bodies",
                        cxxopts::value<int>()->default_value("2000"))(
      "m,method",
      "Method type: DS(direct sum), DFT(depth-first traversal) or "
      "MAC(multipole acceptance criterion)",
      cxxopts::value<std::string>()->default_value("MAC"))(
      "p,param", "MAC parameter (only used with MAC method)",
      cxxopts::value<float>()->default_value("0.5"))(
      "i,iterations", "Number of simulation iterations",
      cxxopts::value<int>()->default_value("10"))(
      "t,timestep", "Simulation time step",
      cxxopts::value<float>()->default_value("0.01"))(
      "s,softening", "Softening factor",
      cxxopts::value<float>()->default_value("1e-9"))(
      "q,quiet", "Suppress timing output",
      cxxopts::value<bool>()->default_value("false"))(
      "h,help",
      "Print usage information")("w,shiftwidth", "Octree shift width",
                                 cxxopts::value<int>()->default_value("16"));

  auto result = options.parse(argc, argv);

  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    return 0;
  }

  // parse options
  int nBodies = result["nbodies"].as<int>();

  std::string method = result["method"].as<std::string>();
  float mac_param = result["param"].as<float>();
  int nIters = result["iterations"].as<int>();
  float dt = result["timestep"].as<float>();
  SOFTENING = result["softening"].as<float>();
  PRINT_TIME = !result["quiet"].as<bool>();
  SHIFT_DIGITS = result["shiftwidth"].as<int>();

  // Set simulation method
  if (method == "DS") {
    USE_TREE = false;
    USE_BH = false;
  } else if (method == "DFT") {
    USE_TREE = true;
    USE_BH = false;
  } else if (method == "MAC") {
    USE_TREE = true;
    USE_BH = true;
  } else {
    std::cerr << "Unknown method: " << method << std::endl;
    std::cout << "Valid methods are: DS, DFT, MAC" << std::endl;
    return -1;
  }

  MAC_PARAM = mac_param;
  printf("Beginning nbody simulation...\n");
  printf(" - %s method \n - %d bodies\n",
         USE_TREE ? (USE_BH ? "MAC" : "DFT") : "DS", nBodies);
  printf(" - MAC parameter: %f\n", MAC_PARAM);
  printf(" - Iterations: %d\n", nIters);
  printf(" - Time step: %f\n", dt);
  printf("-----------------------------------------------\n");

  // TODO: other parameters we might want to control:
  // - SHIFT_DIGITS (controls how many keys/buckets are created)
  // - MAC_PARAM (controls the MAC parameter)
  // - SOFTENING (controls the softening factor)

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
  checkAccuracy(particles, particles_cp, nBodies, nIters);

  delete[] particles;
  delete[] particles_cp;
  return 0;
}
