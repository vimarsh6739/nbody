#include "octree.h"
#include <assert.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctimer.h>
#include <cxxopts.hpp>
#include <iostream>
#include <main.h>
#include <new>
#include <random>
#include <string>

#ifdef ENABLE_CILK
#include <cilk/cilk.h>
#endif
#include <random>

float SOFTENING = 1e-9f;

bool USE_TREE = false; // use tree or list
bool USE_BH = true;    // use barnes hut or not
bool PRINT_TIME = true;
float MAC_PARAM = .5; // MAC parameter

int AXIS_RESOLUTION = 16;
int MAX_KEY_LENGTH = sizeof(Key) * 8;

bool MAC(float target_x, float target_y, float target_z, float x, float y,
         float z, float mass) {
  float dx = target_x - x;
  float dy = target_y - y;
  float dz = target_z - z;
  float distance = sqrtf(dx * dx + dy * dy + dz * dz);
  float size = mass;

  return size / distance < MAC_PARAM;
}

int allInteractionsDFT(Body *bodies, float dt, int nBodies,
                       std::vector<DFTNode> &dft) {
  int nInteractions = 0;

  std::vector<float> new_vx(nBodies), new_vy(nBodies), new_vz(nBodies);

#ifdef ENABLE_CILK
  cilk_for(int target = 0; target < nBodies; target++) {
#else
  for (int target = 0; target < nBodies; target++) {
#endif
    float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;

    for (uint j = 0; j < dft.size(); j++) {
      if (dft[j].isLeaf) {
        for (int bIndex : dft[j].bodies) {
          if (bIndex == target)
            continue;
#ifndef ENABLE_CILK
          nInteractions++;
#endif
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

    new_vx[target] = bodies[target].vx + dt * Fx;
    new_vy[target] = bodies[target].vy + dt * Fy;
    new_vz[target] = bodies[target].vz + dt * Fz;
  }

  for (int i = 0; i < nBodies; i++) {
    bodies[i].vx = new_vx[i];
    bodies[i].vy = new_vy[i];
    bodies[i].vz = new_vz[i];
  }

  return nInteractions;
}

void allInteractionsDS(Body *bodies, float dt, int nBodies) {
  std::vector<float> new_vx(nBodies), new_vy(nBodies), new_vz(nBodies);

#ifdef ENABLE_CILK
  cilk_for(int i = 0; i < nBodies; ++i) {
#else
  for (int i = 0; i < nBodies; ++i) {
#endif
    float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;

    for (int j = 0; j < nBodies; ++j) {
      if (i == j)
        continue;

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

    new_vx[i] = bodies[i].vx + dt * Fx;
    new_vy[i] = bodies[i].vy + dt * Fy;
    new_vz[i] = bodies[i].vz + dt * Fz;
  }

  for (int i = 0; i < nBodies; i++) {
    bodies[i].vx = new_vx[i];
    bodies[i].vy = new_vy[i];
    bodies[i].vz = new_vz[i];
  }
}

int MACInteractionsDFT(Body *bodies, float dt, int nBodies,
                       std::vector<DFTNode> &dft) {
  int nInteractions = 0;

  std::vector<float> new_vx(nBodies), new_vy(nBodies), new_vz(nBodies);

#ifdef ENABLE_CILK
  cilk_for(int target = 0; target < nBodies; target++) {
#else
  for (int target = 0; target < nBodies; target++) {
#endif
    float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;

    for (size_t j = 0; j < dft.size(); j++) {
      if (dft[j].isLeaf) {
        for (int bIndex : dft[j].bodies) {
          if (bIndex == target)
            continue;
#ifndef ENABLE_CILK
          nInteractions++;
#endif
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
#ifndef ENABLE_CILK
        nInteractions++;
#endif
        float dx = dft[j].x - bodies[target].x;
        float dy = dft[j].y - bodies[target].y;
        float dz = dft[j].z - bodies[target].z;
        float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
        float invDist = 1.0f / sqrtf(distSqr);
        float invDist3 = invDist * invDist * invDist;

        Fx += dx * dft[j].mass * invDist3;
        Fy += dy * dft[j].mass * invDist3;
        Fz += dz * dft[j].mass * invDist3;

        if (dft[j].autorope < dft.size()) {
          j = dft[j].autorope - 1;
        } else {
          break;
        }
      }
    }

    new_vx[target] = bodies[target].vx + dt * Fx;
    new_vy[target] = bodies[target].vy + dt * Fy;
    new_vz[target] = bodies[target].vz + dt * Fz;
  }

  for (int i = 0; i < nBodies; i++) {
    bodies[i].vx = new_vx[i];
    bodies[i].vy = new_vy[i];
    bodies[i].vz = new_vz[i];
  }

  return nInteractions;
}

void reconstructOctree(Octree *&octree, std::vector<DFTNode> &dft,
                       Body *particles, int nBodies) {

  // regenerate keys for new positions
  for (int i = 0; i < nBodies; ++i) {
    particles[i].key = computeMortonKey(particles[i], AXIS_RESOLUTION);
  }

  // free old octree, create new octree (will optimize later)
  delete octree;
  octree = new Octree(AXIS_RESOLUTION);

  for (int i = 0; i < nBodies; i++) {
    octree->insert(particles[i]);
  }
  octree->buildDFT(dft, particles);

  // printf("Octree built with %d buckets (leaves with unique key)\n",
  //        nUniqueLeaves);
  // octree->printTree(octree->root, 0);
  // octree->printTree(octree->root, 0);
}

// Randomize the positions and velocities in the range [-1, 1],
// and assign a mass (every 7th float) in the range [0.1, 1.1].
void randomizeBodies(PhiloxEngine &rng, Body *bodies, int n) {
  std::uniform_real_distribution<double> x_dis(0, 0.99);
  std::uniform_real_distribution<double> v_dis(0, 1);

  for (int i = 0; i < n; i++) {
    Body &body = bodies[i];
    body.x = x_dis(rng);
    body.y = x_dis(rng);
    body.z = x_dis(rng);

    body.index = i;
    // set the key for the current body
    body.key = computeMortonKey(body, AXIS_RESOLUTION);

    body.vx = (v_dis(rng)) * 2 - 1;
    body.vy = (v_dis(rng)) * 2 - 1;
    body.vz = (v_dis(rng)) * 2 - 1;
    body.m = (v_dis(rng)) + 0.1;
  }
}

float checkAccuracy(Body *p, Body *orig, int nBodies, int nIters) {
  printf("-----------------------------------------------\n");
  printf("ACCURACY CHECK AGAINST DS\n");

  USE_TREE = false;
  USE_BH = false;
  PRINT_TIME = false;
  nbodyIterate(orig, 0.01f, nBodies, nIters);
  PRINT_TIME = true;

  // Compute RMS error(standard check for positional accuracy)
  double rmsError = 0.0;
  for (int i = 0; i < nBodies; i++) {
    double dx = p[i].x - orig[i].x;
    double dy = p[i].y - orig[i].y;
    double dz = p[i].z - orig[i].z;
    rmsError += dx * dx + dy * dy + dz * dz;
  }

  rmsError = std::sqrt((rmsError / nBodies));

  return rmsError;
}

void integratePositions(Body *p, float dt, int start, int end) {
  for (int i = start; i < end; i++) {
    p[i].x += (p[i].vx * dt);
    p[i].y += (p[i].vy * dt);
    p[i].z += (p[i].vz * dt);
  }
}

double nbodyIterate(Body *p, float dt, int nBodies, int nIters) {

  // begin benchmark
  ctimer_t timer;
  ctimer_start(&timer);

  std::vector<DFTNode> dft;
  Octree *octree = nullptr;

  if (USE_TREE) {
    octree = new Octree(AXIS_RESOLUTION);
  }

  for (int iter = 1; iter <= nIters; iter++) {

    if (USE_TREE) {
      reconstructOctree(octree, dft, p, nBodies);
    }

    int nInteractions = bodyForce(p, dt, nBodies, dft);
    integratePositions(p, dt, 0, nBodies);
  }

  if (USE_TREE) {
    delete octree;
  }

  ctimer_stop(&timer);
  ctimer_measure(&timer);
  double tElapsed = timespec_sec(timer.elapsed);

  // Report time
  if (PRINT_TIME) {
    printf("Total time = %.6f seconds\n", tElapsed);
  }

  return tElapsed;
}

int bodyForce(Body *p, float dt, int n, std::vector<DFTNode> dft) {
  if (USE_TREE && USE_BH)
    return MACInteractionsDFT(p, dt, n, dft);
  else if (USE_TREE)
    return allInteractionsDFT(p, dt, n, dft);
  else {
    allInteractionsDS(p, dt, n);
    return n * n;
  }
}

int libMain(const int argc, const char **argv) {

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
      "h,help", "Print usage information")(
      "w,resolution", "Octree bit-resolution (same for x,y,z axis)",
      cxxopts::value<int>()->default_value("16"))("e,error",
                                                  "Check for accuracy");

  auto result = options.parse(argc, argv);

  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    return 0;
  }

  bool errorCheck = false;
  if (result.count("error")) {
    errorCheck = true;
  }

  printf("Error check: %s\n", errorCheck ? "true" : "false");

  // parse options
  int nBodies = result["nbodies"].as<int>();

  std::string method = result["method"].as<std::string>();
  float mac_param = result["param"].as<float>();
  int nIters = result["iterations"].as<int>();
  float dt = result["timestep"].as<float>();
  SOFTENING = result["softening"].as<float>();
  PRINT_TIME = !result["quiet"].as<bool>();
  AXIS_RESOLUTION = result["resolution"].as<int>();

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
    std::cerr << "Valid methods are: DS, DFT, MAC" << std::endl;
    return EXIT_FAILURE;
  }

  MAC_PARAM = mac_param;
  printf("Beginning nbody simulation...\n");
  printf(" - %s method \n - %d bodies\n",
         USE_TREE ? (USE_BH ? "MAC" : "DFT") : "DS", nBodies);
  printf(" - MAC parameter: %f\n", MAC_PARAM);
  printf(" - Iterations: %d\n", nIters);
  printf(" - Time step: %f\n", dt);
  printf("-----------------------------------------------\n");

  // Allocate memory for nBodies
  Body *particles = new Body[nBodies];

  // Initialize positions, velocities, and masses (7 values per body)
  PhiloxEngine rng(2025);
  randomizeBodies(rng, particles, nBodies);

  // make a copy of bodies for accuracy check
  Body *particles_cp = new Body[nBodies];
  for (int i = 0; i < nBodies; i++) {
    particles_cp[i] = particles[i];
  }

  nbodyIterate(particles, dt, nBodies, nIters);

  if (errorCheck) {
    float err = checkAccuracy(particles, particles_cp, nBodies, nIters);
    printf("RMS error: %f\n", err);
  }

  delete[] particles;
  delete[] particles_cp;
  return 0;
}
