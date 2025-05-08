#include <ctimer.h>
#include <cxxopts.hpp>
#include <main.h>
#include <assert.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
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

int SHIFT_DIGITS = 16;
int MAX_KEY_LENGTH = sizeof(Key) * 8;

#ifdef ENABLE_CUDA
extern bool CUDA = true;
#else
extern bool CUDA = false;
#endif

bool MAC(float target_x, float target_y, float target_z, float x, float y,
         float z, float mass) {

  float distance =
      sqrtf((target_x - x) * (target_x - x) + (target_y - y) * (target_y - y) +
            (target_z - z) * (target_z - z));
  float source_mass = mass;

  bool result = distance > MAC_PARAM;

  return distance / source_mass > MAC_PARAM;
}

int allInteractionsDFT(Body *bodies, float dt, int nBodies,
                       std::vector<DFTNode> &dft) {

  // iterate over all bodies (targets)
  int nInteractions = 0;

  Body *bodies_copy = new Body[nBodies];
  for (int i = 0; i < nBodies; ++i) {
    bodies_copy[i] = bodies[i];
  }

#ifdef ENABLE_CILK
  cilk_for(int target = 0; target < nBodies; target++) {
#else
  for (int target = 0; target < nBodies; target++) {
#endif
    float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;

    for (int j = 0; j < dft.size(); j++) {
      if (dft[j].isLeaf) {
        for (int bIndex : dft[j].bodies) {
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
    bodies_copy[target].vx += dt * Fx;
    bodies_copy[target].vy += dt * Fy;
    bodies_copy[target].vz += dt * Fz;
  }

  delete[] bodies;
  bodies = bodies_copy;
  return nInteractions;
}

void allInteractionsDS(Body *bodies, float dt, int nBodies) {
  int n = nBodies;

  Body *bodies_copy = new Body[n];
  for (int i = 0; i < n; ++i) {
    bodies_copy[i] = bodies[i];
  }

#ifdef ENABLE_CILK
  cilk_for(int i = 0; i < n; ++i) {
#else
  for (int i = 0; i < n; ++i) {
#endif
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
    bodies_copy[i].vx += dt * Fx;
    bodies_copy[i].vy += dt * Fy;
    bodies_copy[i].vz += dt * Fz;
  }

  delete[] bodies;
  bodies = bodies_copy;
}

int MACInteractionsDFT(Body *bodies, float dt, int nBodies,
                       std::vector<DFTNode> &dft) {

  // iterate over all bodies (targets)

  // iterate over all bodies (targets)
  int nInteractions = 0;

  Body *bodies_copy = new Body[nBodies];
  for (int i = 0; i < nBodies; ++i) {
    bodies_copy[i] = bodies[i];
  }

#ifdef ENABLE_CILK
  cilk_for(int target = 0; target < nBodies; target++) {
#else
  for (int target = 0; target < nBodies; target++) {
#endif
    float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;
    for (int j = 0; j < dft.size(); j++) {
      if (dft[j].isLeaf) {
        for (int bIndex : dft[j].bodies) {

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
        j = dft[j].autorope - 1; // -1 because we increment j in the for loop
      }
    }

    bodies_copy[target].vx += dt * Fx;
    bodies_copy[target].vy += dt * Fy;
    bodies_copy[target].vz += dt * Fz;
  }

  delete[] bodies;
  bodies = bodies_copy;
  return nInteractions;
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
      "w,shiftwidth", "Octree shift width",
      cxxopts::value<int>()->default_value("16"))("e,error", "Enable flag");

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

