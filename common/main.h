#ifndef MAIN_H
#define MAIN_H

#include "octree.h"
#include "philox_engine.h"
#include <vector>

extern float SOFTENING;
#define POSMAX 100

extern bool USE_TREE;
extern bool USE_BH;
extern bool PRINT_TIME;
extern float THETA;

extern int AXIS_RESOLUTION;
extern int MAX_KEY_LENGTH;

extern bool CUDA;

void randomizeBodies(PhiloxEngine &rng, Body *bodies, int n);

/**
 * @brief Construct/Re-construct octree
 *
 * @param octree
 * @param bodies
 * @param nBodies
 */
void reconstructOctree(Octree *&octree, Body *bodies, int nBodies);

/**
 * @brief Run simulation
 *
 * @param p
 * @param dt
 * @param nBodies
 * @param nIters
 * @return
 */
double nbodyIterate(Body *p, float dt, int nBodies, int nIters);

/**
 * @brief Validate accuracy against all-to-all interactions
 *
 * @param p
 * @param orig
 * @param nBodies
 * @param nIters
 * @return
 */
float computeRmsError(Body *particles, Body *orig, int nBodies, int nIters,
                      float dt);

/**
 * @brief Compute new positions after applying induced momentum
 *
 * @param p
 * @param dt
 * @param start
 * @param end
 */
void integratePositions(Body *particles, float dt, int start, int end);

/**
 * @brief Compute forces using multiple methods
 *
 * @param tree
 * @param particles
 * @param dt
 * @param nbodies
 */
void bodyForce(Octree *&tree, Body *particles, float dt, int nbodies);

/**
 * @brief Compute forces using all-to-all interactions
 *
 * @param bodies
 * @param dt
 * @param nBodies
 */
void allInteractionsDS(Body *bodies, float dt, int nBodies);

/**
 * @brief Depth first iteration
 *
 * @param tree
 * @param node
 * @param particle
 * @param Fx
 * @param Fy
 * @param Fz
 * @param level
 */
void BarnesHutDFS(Octree *&tree, Node *&node, Body &particle, float &Fx,
                  float &Fy, float &Fz, int level);
/**
 * @brief Compute forces using truncated Barnes-Hut computations
 *
 * @param tree
 * @param particles
 * @param dt
 * @param nbodies
 */
void BarnesHutInteractions(Octree *&tree, Body *&particles, float dt,
                           int nbodies);

// int MACInteractionsDFT(Body *bodies, float dt, int nBodies);
// int allInteractionsDFT(Body *bodies, float dt, int nBodies);

/**
 * @brief Driver code
 *
 * @param argc
 * @param argv
 * @return
 */
int libMain(const int argc, const char **argv);

#endif // !MAIN_H
