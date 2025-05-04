#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "ctimer.h"
#include "philox_engine.h"
#include <cuda_runtime.h>

#include "octree.h"

using namespace std;

#define SOFTENING 1e-5f
#define BLOCK_SIZE 256
#define MAC_PARAM 0.6f

#define CUDA_CALL(call)                                                        \
  do {                                                                         \
    cudaError_t err = call;                                                    \
    if (err != cudaSuccess) {                                                  \
      fprintf(stderr, "CUDA error %s@%d: %s (%d)\n", __FILE__, __LINE__,       \
              cudaGetErrorString(err), err);                                   \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  } while (0)

__device__ inline bool MAC(float tx, float ty, float tz, float nx, float ny,
                           float nz, float node_size_sq, float theta_sq) {
  float dx = nx - tx;
  float dy = ny - ty;
  float dz = nz - tz;
  float r_sq = dx * dx + dy * dy + dz * dz;
  if (r_sq < 1e-18f)
    return false;
  return node_size_sq < theta_sq * r_sq;
}

__global__ void MACForceKernel(Body *bodies, const DFTNode *dftNodes, int nDFT,
                               float dt, float force_scale, int nBodies) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= nBodies)
    return;
  float Fx = 0, Fy = 0, Fz = 0;
  float tx = bodies[i].x, ty = bodies[i].y, tz = bodies[i].z;
  const float theta_sq = MAC_PARAM * MAC_PARAM;

  for (int j = 0; j < nDFT;) {
    DFTNode node = dftNodes[j];
    if (node.nBodies <= 0 || node.mass <= 1e-20f ||
        (node.isLeaf && node.singleBodyIndex == i)) {
      ++j;
      continue;
    }
    bool approx = !node.isLeaf && MAC(tx, ty, tz, node.x, node.y, node.z,
                                      node.size * node.size, theta_sq);
    if (node.isLeaf || approx) {
      float dx = node.x - tx;
      float dy = node.y - ty;
      float dz = node.z - tz;
      float r2 = dx * dx + dy * dy + dz * dz + SOFTENING;
      float invR = rsqrtf(r2);
      float invR3 = invR * invR * invR;
      Fx += dx * node.mass * invR3;
      Fy += dy * node.mass * invR3;
      Fz += dz * node.mass * invR3;
      j = approx ? node.autorope : j + 1;
    } else {
      j++;
    }
  }
  bodies[i].vx += (dt * Fx) * force_scale;
  bodies[i].vy += (dt * Fy) * force_scale;
  bodies[i].vz += (dt * Fz) * force_scale;
}

__global__ void DirectSumForceKernel(Body *bodies, float dt, float force_scale,
                                     int n) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n)
    return;
  float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;
  float tx = bodies[i].x;
  float ty = bodies[i].y;
  float tz = bodies[i].z;
  for (int j = 0; j < n; j++) {
    if (i == j)
      continue;
    float dx = bodies[j].x - tx;
    float dy = bodies[j].y - ty;
    float dz = bodies[j].z - tz;
    float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
    float invDist = rsqrtf(distSqr);
    float invDist3 = invDist * invDist * invDist;
    Fx += dx * bodies[j].m * invDist3;
    Fy += dy * bodies[j].m * invDist3;
    Fz += dz * bodies[j].m * invDist3;
  }
  bodies[i].vx += (dt * Fx) * force_scale;
  bodies[i].vy += (dt * Fy) * force_scale;
  bodies[i].vz += (dt * Fz) * force_scale;
}

__global__ void leapfrogDriftKernel(Body *p, float dt, int n) {
  int id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id < n) {
    p[id].x += p[id].vx * dt;
    p[id].y += p[id].vy * dt;
    p[id].z += p[id].vz * dt;
  }
}

int main(int argc, char **argv) {
  // Set precision for floating point output if needed later
  // cout << fixed << setprecision(6); cerr << fixed << setprecision(6);

  printf("NBody Simulation\n"); // Match old header

  int N = 10000; // Default value
  if (argc > 1) {
    try {
      N = stoi(argv[1]);
    } catch (const std::exception &e) {
      cerr << "Error parsing N=" << argv[1] << ": " << e.what() << endl;
      return 1;
    }
  }
  if (N <= 0) {
    cerr << "Error: N must be > 0" << endl;
    return 1;
  }

  // --- Parameters ---
  const int nIter = 10;
  const float dt = 0.005f;
  // Force method to MAC for this output style
  printf(" - MAC method \n - %d bodies\n", N);
  printf(
      "-----------------------------------------------\n"); // Match separator

  // --- Host Init & Alloc ---
  PhiloxEngine rng(2025);
  size_t bodyBytes = N * sizeof(Body);
  vector<Body> h_bodies(N);
  vector<Body> h_bodies_ref(N);
  randomInit(h_bodies.data(), N, rng);
  memcpy(h_bodies_ref.data(), h_bodies.data(), bodyBytes);

  // --- GPU Alloc ---
  Body *d_bodies = nullptr;
  Body *d_bodies_ref = nullptr;
  DFTNode *d_dft = nullptr;
  CUDA_CALL(cudaMalloc(&d_bodies, bodyBytes));
  CUDA_CALL(cudaMalloc(&d_bodies_ref, bodyBytes));
  size_t dftCapacityBytes = 0;

  // --- Initial Setup ---
  int maxBits = 0;
  if (N > 0) {
    for (const auto &b : h_bodies) {
      maxBits = max(maxBits, binaryLength(mortonNoPrepend(b)));
    }
  }
  Key prepend = (maxBits > 0) ? ((Key)1 << maxBits) : 1;
  int keyBits = maxBits + 1;
  for (auto &b : h_bodies) {
    b.key = mortonNoPrepend(b) + prepend;
  }
  vector<DFTNode> h_dft;
  try {
    OctreeBuilder builder(keyBits, N);
    for (auto &b : h_bodies) {
      builder.insert(b);
    }
    builder.buildDFT(h_dft, h_bodies.data());
  } catch (const exception &e) {
    cerr << "Initial build error: " << e.what() << endl;
    return 1;
  }

  // Print approximate Octree build message
  printf("Octree built, DFT nodes: %zu\n", h_dft.size());

  size_t currentDftBytes = h_dft.size() * sizeof(DFTNode);
  if (currentDftBytes > dftCapacityBytes) {
    if (d_dft) {
      CUDA_CALL(cudaFree(d_dft));
    }
    d_dft = nullptr;
    dftCapacityBytes = 0;
    if (currentDftBytes > 0) {
      CUDA_CALL(cudaMalloc(&d_dft, currentDftBytes));
      dftCapacityBytes = currentDftBytes;
    }
  }
  CUDA_CALL(
      cudaMemcpy(d_bodies, h_bodies.data(), bodyBytes, cudaMemcpyHostToDevice));
  if (currentDftBytes > 0 && d_dft) {
    CUDA_CALL(cudaMemcpy(d_dft, h_dft.data(), currentDftBytes,
                         cudaMemcpyHostToDevice));
  }
  int gridSize = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
  if (d_dft != nullptr || h_dft.empty()) {
    MACForceKernel<<<gridSize, BLOCK_SIZE>>>(d_bodies, d_dft, h_dft.size(), dt,
                                             0.5f, N);
    CUDA_CALL(cudaGetLastError());
    CUDA_CALL(cudaDeviceSynchronize());
  } else {
    cerr << "Error: Cannot calc initial kick.\n";
    return 1;
  }

  // --- Simulation Loop ---
  ctimer_t loop_timer;
  ctimer_start(&loop_timer);

  for (int step = 0; step < nIter; ++step) {
    // 1. Drift
    leapfrogDriftKernel<<<gridSize, BLOCK_SIZE>>>(d_bodies, dt, N);
    CUDA_CALL(cudaGetLastError());
    CUDA_CALL(cudaDeviceSynchronize());

    // 2. Copy Back, Rebuild Tree/DFT
    CUDA_CALL(cudaMemcpy(h_bodies.data(), d_bodies, bodyBytes,
                         cudaMemcpyDeviceToHost));
    for (auto &b : h_bodies) {
      b.key = mortonNoPrepend(b) + prepend;
    }
    h_dft.clear();
    try {
      OctreeBuilder builder(keyBits, N);
      for (auto &b : h_bodies) {
        builder.insert(b);
      }
      builder.buildDFT(h_dft, h_bodies.data());
    } catch (const exception &e) {
      cerr << "Build error step " << step << ": " << e.what() << endl;
      break;
    }

    // 3. Resize & Copy H->D
    currentDftBytes = h_dft.size() * sizeof(DFTNode);
    if (currentDftBytes > dftCapacityBytes) {
      if (d_dft) {
        CUDA_CALL(cudaFree(d_dft));
      }
      d_dft = nullptr;
      dftCapacityBytes = 0;
      if (currentDftBytes > 0) {
        CUDA_CALL(cudaMalloc(&d_dft, currentDftBytes));
        dftCapacityBytes = currentDftBytes;
      }
    }
    if (currentDftBytes > 0 && d_dft) {
      CUDA_CALL(cudaMemcpy(d_dft, h_dft.data(), currentDftBytes,
                           cudaMemcpyHostToDevice));
    }

    // 4. Kick
    if (d_dft != nullptr || h_dft.empty()) {
      MACForceKernel<<<gridSize, BLOCK_SIZE>>>(d_bodies, d_dft, h_dft.size(),
                                               dt, 1.0f, N);
      CUDA_CALL(cudaGetLastError());
      CUDA_CALL(cudaDeviceSynchronize());
    } else {
      cerr << "Error: Skipping kick step " << step << ".\n";
      break;
    }
  }
  ctimer_stop(&loop_timer);
  ctimer_measure(&loop_timer);
  double totalLoopTime = timespec_sec(loop_timer.elapsed);
  CUDA_CALL(cudaMemcpy(h_bodies.data(), d_bodies, bodyBytes,
                       cudaMemcpyDeviceToHost)); // Final copy back

  // Print average time per iteration (excluding initial build/kick)
  // Note: Original code skipped iter 1 for timing avg; this includes all loop
  // time.
  if (nIter > 0) {
    printf("Average time per iteration: %.3f seconds.\n",
           totalLoopTime / nIter);
  }
  // Skip "Billion Interactions" line

  // --- Direct Summation Reference ---
  CUDA_CALL(cudaMemcpy(d_bodies_ref, h_bodies_ref.data(), bodyBytes,
                       cudaMemcpyHostToDevice));
  DirectSumForceKernel<<<gridSize, BLOCK_SIZE>>>(d_bodies_ref, dt, 0.5f, N);
  CUDA_CALL(cudaGetLastError());
  CUDA_CALL(cudaDeviceSynchronize());
  ctimer_t ds_timer;
  ctimer_start(&ds_timer);
  for (int iter = 0; iter < nIter; iter++) {
    leapfrogDriftKernel<<<gridSize, BLOCK_SIZE>>>(d_bodies_ref, dt, N);
    DirectSumForceKernel<<<gridSize, BLOCK_SIZE>>>(d_bodies_ref, dt, 1.0f, N);
    CUDA_CALL(cudaGetLastError());
  }
  CUDA_CALL(cudaDeviceSynchronize());
  ctimer_stop(&ds_timer);
  ctimer_measure(&ds_timer);
  // No DS timing print here to match target output

  CUDA_CALL(cudaMemcpy(h_bodies_ref.data(), d_bodies_ref, bodyBytes,
                       cudaMemcpyDeviceToHost));

  // --- Accuracy Check ---
  printf(
      "-----------------------------------------------\n"); // Match separator
  printf("ACCURACY CHECK AGAINST DS\n");                    // Match header

  double maxError = 0.0;
  double totalError = 0.0;
  for (int i = 0; i < N; ++i) {
    double dx = (double)h_bodies[i].x - (double)h_bodies_ref[i].x;
    double dy = (double)h_bodies[i].y - (double)h_bodies_ref[i].y;
    double dz = (double)h_bodies[i].z - (double)h_bodies_ref[i].z;
    double err = sqrt(dx * dx + dy * dy + dz * dz);
    if (err > maxError)
      maxError = err;
    totalError += err;
  }
  // double avgError = (N > 0) ? totalError / N : 0.0; // Avg error not printed
  // in target

  // Print only max error in the old format
  printf("Max error: %f\n", (float)maxError);

  // --- Cleanup ---
  if (d_bodies)
    cudaFree(d_bodies);
  if (d_bodies_ref)
    cudaFree(d_bodies_ref);
  if (d_dft)
    cudaFree(d_dft);

  return 0; // No final "Done." message
}
