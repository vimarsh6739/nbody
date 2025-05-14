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
extern float MAC_PARAM;

extern int AXIS_RESOLUTION;
extern int MAX_KEY_LENGTH;

extern bool CUDA;

void randomizeBodies(PhiloxEngine &rng, Body *bodies, int n);
void reconstructOctree(Octree *&octree, Body *bodies, int nBodies);
double nbodyIterate(Body *p, float dt, int nBodies, int nIters);
float checkAccuracy(Body *p, Body *orig, int nBodies, int nIters);
void integratePositions(Body *p, float dt, int start, int end);
void bodyForce(Body *p, float dt, int n);
void allInteractionsDS(Body *bodies, float dt, int nBodies);
int MACInteractionsDFT(Body *bodies, float dt, int nBodies);
int allInteractionsDFT(Body *bodies, float dt, int nBodies);
int libMain(const int argc, const char **argv);

#endif // !MAIN_H
