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

void BarnesHutDFS(Octree *&tree, Node *&node, Body &particle, float &Fx,
                  float &Fy, float &Fz, int level) {
  if (tree->isEmpty(node))
    return;

  // compute approx forces
  float dx = node->cx - particle.x;
  float dy = node->cy - particle.y;
  float dz = node->cz - particle.z;
  float dstSq = dx * dx + dy * dy + dz * dz + SOFTENING;
  float invDist = 1.0f / sqrtf(dstSq);

  float octantSize = 1.0f / (1 << level);
  float MAC = octantSize * invDist;

  // multipole acceptance criterion
  if (tree->isLeaf(node) || (USE_BH && (MAC < THETA))) {
    if (node->key != particle.key) {
      float invDist3 = invDist * invDist * invDist;
      Fx += dx * node->tm * invDist3;
      Fy += dy * node->tm * invDist3;
      Fz += dz * node->tm * invDist3;
    }
  } else {
    for (int i = 0; i < 8; ++i) {
      // sparse DFS is possible here
      if (node->maskChildren & (1 << i)) {
        BarnesHutDFS(tree, node->children[i], particle, Fx, Fy, Fz, (level + 1));
      }
    }
  }
}

void BarnesHutInteractions(Octree *&tree, Body *&particles, float dt,
                           int nbodies) {
  for (int i = 0; i < nbodies; ++i) {
    float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;

    // perform a truncated DFS to compute Fx,Fy,Fz
    BarnesHutDFS(tree, tree->root, particles[i], Fx, Fy, Fz, 0);

    // update momentum
    particles[i].vx += Fx * dt;
    particles[i].vy += Fy * dt;
    particles[i].vz += Fz * dt;
  }
}

/* int allInteractionsDFT(Body *bodies, float dt, int nBodies,
                       std::vector<DFTNode> &dft) {
  int nInteractions = 0;

  std::vector<float> new_vx(nBodies), new_vy(nBodies), new_vz(nBodies);

  for (int target = 0; target < nBodies; target++) {
    float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;

    for (uint j = 0; j < dft.size(); j++) {
      if (dft[j].isLeaf) {
        for (int bIndex : dft[j].bodies) {
          if (bIndex == target)
            continue;
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
} */

void allInteractionsDS(Body *bodies, float dt, int nBodies) {
  std::vector<float> new_vx(nBodies), new_vy(nBodies), new_vz(nBodies);

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

/* int MACInteractionsDFT(Body *bodies, float dt, int nBodies,
                       std::vector<DFTNode> &dft) {
  int nInteractions = 0;

  std::vector<float> new_vx(nBodies), new_vy(nBodies), new_vz(nBodies);

  cilk_for(int target = 0; target < nBodies; target++) {
    float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;

    for (size_t j = 0; j < dft.size(); j++) {
      if (dft[j].isLeaf) {
        for (int bIndex : dft[j].bodies) {
          if (bIndex == target)
            continue;
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
} */

void reconstructOctree(Octree *&tree, Body *particles, int nBodies) {
  for (int i = 0; i < nBodies; ++i) {
    particles[i].key = computeMortonKey(particles[i], AXIS_RESOLUTION);
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

void wrap(float &pos, float &vel, float minv, float maxv) {
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
    wrap(p[i].x, p[i].vx);
    wrap(p[i].y, p[i].vy);
    wrap(p[i].z, p[i].vz);
    p[i].x += (p[i].vx * dt);
    p[i].y += (p[i].vy * dt);
    p[i].z += (p[i].vz * dt);
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

float computeRmsError(Body *p, Body *orig, int nBodies, int nIters, float dt) {

  printf("-----------------------------------------------\n");
  printf("COMPUTING RMS ERROR AGAINST DIRECT SUMMATION\n");

  USE_TREE = false;
  USE_BH = false;
  PRINT_TIME = false;
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
