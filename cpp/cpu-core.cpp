#include "octree.h"
#include <assert.h>
#include <cilk/cilk.h>
#include <cmath>
#include <cstdint>
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

float SOFTENING = 1e-9f;

bool USE_TREE = false; // use tree or list
bool USE_BH = true;    // use barnes hut or not
bool PRINT_TIME = true;
float THETA = .5; // MAC parameter

int AXIS_RESOLUTION = 16;
int MAX_KEY_LENGTH = sizeof(Key) * 8;
int bh_ctr = 0;
int ds_ctr = 0;
void BarnesHutDFS(Octree *&tree, Node *&node, Body *&particles, int tid, float &Fx,
                  float &Fy, float &Fz, int level) {
  if (tree->isEmpty(node))
    return;
  // printf("level: %d, tid: %d\n", level, tid);
  // compute approx forces
  float dx = node->cx - particles[tid].x;
  float dy = node->cy - particles[tid].y;
  float dz = node->cz - particles[tid].z;
  float dstSq = dx * dx + dy * dy + dz * dz + SOFTENING;
  float invDist = 1.0f / sqrtf(dstSq);

  float octantSize = 1.0f / (1 << level);
  float MAC = octantSize * invDist;

  if(tree->isLeaf(node)) {
    for(int idx: node->bodyIdx){
      
      if(idx==tid) continue;
      float dx =  particles[idx].x - particles[tid].x;
      float dy =  particles[idx].y - particles[tid].y;
      float dz =  particles[idx].z - particles[tid].z;
      float dstSq = dx * dx + dy * dy + dz * dz + SOFTENING;
      float invDist = 1.0f / sqrtf(dstSq);
      float invDist3 = invDist * invDist * invDist;
      
      Fx += dx * particles[idx].m * invDist3;
      Fy += dy * particles[idx].m * invDist3;
      Fz += dz * particles[idx].m * invDist3;
      bh_ctr++;
    }
  } else if (USE_BH && (MAC < THETA)) {
    // multipole acceptance criterion
    if (node->key != particles[tid].key) {
      float invDist3 = invDist * invDist * invDist;
      Fx += dx * node->tm * invDist3;
      Fy += dy * node->tm * invDist3;
      Fz += dz * node->tm * invDist3;
    }
  } else {
    for (int i = 0; i < 8; ++i) {
      if (node->maskChildren & (1 << i)) {
        BarnesHutDFS(tree, node->children[i], particles, tid, Fx, Fy, Fz,
                    (level + 1));
      }
    }
  }
}

void BarnesHutInteractions(Octree *&tree, Body *&particles, float dt,
                           int nbodies) {
  for (int i = 0; i < nbodies; ++i) {
    float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;

    // perform a truncated DFS to compute Fx,Fy,Fz
    BarnesHutDFS(tree, tree->root, particles, i, Fx, Fy, Fz, 0);

    // update momentum
    particles[i].vx += Fx * dt;
    particles[i].vy += Fy * dt;
    particles[i].vz += Fz * dt;
  }
}

void allInteractionsDS(Body *bodies, float dt, int nBodies) {
  for (int i = 0; i < nBodies; ++i) {
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
      ds_ctr++;
    }

    bodies[i].vx += dt * Fx;
    bodies[i].vy += dt * Fy;
    bodies[i].vz += dt * Fz;
  }
}

void reconstructOctree(Octree *&tree, Body *particles, int nBodies) {
  for (int i = 0; i < nBodies; ++i) {
    // Normalize to [0,1] for Morton key computation
    float norm_x = particles[i].x / 10.0f;
    float norm_y = particles[i].y / 10.0f;
    float norm_z = particles[i].z / 10.0f;
    particles[i].key = computeMortonKey({norm_x, norm_y, norm_z, 0, 0, 0, 1.0f, 0, 0}, AXIS_RESOLUTION);
  }

  // free old octree, create new octree (will optimize later)
  delete tree;
  tree = new Octree(AXIS_RESOLUTION);

  for (int i = 0; i < nBodies; i++) {
    tree->insert(particles[i]);
  }
  uint64_t root_size = 0;
  tree->finalizeStats(tree->root, root_size);

  // tree->printTree(tree->root, 0);
}

void wrap(float &pos, float &vel, float minv , float maxv ) {
  // Check for NaN
  if (std::isnan(pos)) {
    pos = minv;  // Reset to minimum value
    vel = fabsf(vel);  // Make velocity positive
  }
  
  if (pos <= minv) {
    pos = minv;
    vel = fabsf(vel);
  } else if (pos >= maxv) {
    pos = maxv;
    vel = -fabsf(vel);
  }
}

void integratePositions(Body *&p, float dt, int start, int end) {
  for (int i = start; i < end; i++) {
    p[i].x += (p[i].vx * dt);
    p[i].y += (p[i].vy * dt);
    p[i].z += (p[i].vz * dt);
    wrap(p[i].x, p[i].vx);
    wrap(p[i].y, p[i].vy);
    wrap(p[i].z, p[i].vz);
  }
}

void bodyForce(Octree *&tree, Body *p, float dt, int n) {
  if (USE_TREE) {
    BarnesHutInteractions(tree, p, dt, n);
  } else {
    allInteractionsDS(p, dt, n);
  }
}

double nbodyIterate(Body *particles, float dt, int nBodies, int nIters) {

  printf("Beginning nbody simulation...\n");
  printf(" - %s method \n - %d bodies\n",
         USE_TREE ? (USE_BH ? "MAC" : "DFT") : "DS", nBodies);
  printf(" - MAC parameter: %f\n", THETA);
  printf(" - Iterations: %d\n", nIters);
  printf(" - Time step: %f\n", dt);
  printf("-----------------------------------------------\n");

  // begin benchmark
  ctimer_t timer;
  ctimer_start(&timer);

  Octree *octree = nullptr;

  if (USE_TREE) {
    octree = new Octree(AXIS_RESOLUTION);
  }

  for (int iter = 1; iter <= nIters; iter++) {

    if (USE_TREE) {
      reconstructOctree(octree, particles, nBodies);
    }

    bodyForce(octree, particles, dt, nBodies);
    integratePositions(particles, dt, 0, nBodies);
  }

  if (USE_TREE) {
    delete octree;
  }

  ctimer_stop(&timer);
  ctimer_measure(&timer);
  double tElapsed = timespec_sec(timer.elapsed);

  // Report time
  if (PRINT_TIME) {
    printf("Total time = %.6f seconds for %d iterations\n", tElapsed, nIters);
  }

  return tElapsed;
}

// Randomize the positions and velocities in the range [-1, 1],
// and assign a mass (every 7th float) in the range [0.1, 1.1].
void randomizeBodies(PhiloxEngine &rng, Body *bodies, int n) {
  std::uniform_real_distribution<float> x_dis(0.0, 10.0);  // Changed to [0,10]
  std::uniform_real_distribution<float> v_dis(0, 1);

  for (int i = 0; i < n; i++) {
    Body &body = bodies[i];
    body.x = x_dis(rng);
    body.y = x_dis(rng);
    body.z = x_dis(rng);

    body.index = i;
    // Normalize to [0,1] for Morton key computation
    float norm_x = body.x / 10.0f;
    float norm_y = body.y / 10.0f;
    float norm_z = body.z / 10.0f;
    body.key = computeMortonKey({norm_x, norm_y, norm_z, 0, 0, 0, 1.0f, 0, 0}, AXIS_RESOLUTION);

    body.vx = (v_dis(rng)) * 2 - 1;
    body.vy = (v_dis(rng)) * 2 - 1;
    body.vz = (v_dis(rng)) * 2 - 1;
    body.m = (v_dis(rng)) + 0.1;
  }
}

float computeRmsError(Body *p, Body *orig, int nBodies, int nIters, float dt) {

  printf("-----------------------------------------------\n");
  printf("COMPUTING RMS ERROR AGAINST DIRECT SUMMATION\n");

  USE_TREE = false;
  USE_BH = false;
  PRINT_TIME = true;
  nbodyIterate(orig, dt, nBodies, nIters);
  PRINT_TIME = true;

  // Compute RMS error(standard check for positional accuracy)
  double rmsError = 0.0;
  for (int i = 0; i < nBodies; i++) {
    double dx = static_cast<double>(p[i].x - orig[i].x);
    double dy = static_cast<double>(p[i].y - orig[i].y);
    double dz = static_cast<double>(p[i].z - orig[i].z);
    rmsError += dx * dx + dy * dy + dz * dz;
  }

  rmsError = std::sqrt((rmsError / nBodies));

  printf("DS interactions: %d\n", ds_ctr);
  printf("BH interactions: %d\n", bh_ctr);
  return rmsError;
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

  THETA = mac_param;

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
    float err = computeRmsError(particles, particles_cp, nBodies, nIters, dt);
    printf("RMS error: %f\n", err);
  }

  delete[] particles;
  delete[] particles_cp;
  return 0;
}
