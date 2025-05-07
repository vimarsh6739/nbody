#include "octree.h"
#include "philox_engine.h"

#define POSMAX 100
extern float SOFTENING;
extern bool USE_TREE;
extern bool USE_BH;
extern bool PRINT_TIME;
extern float MAC_PARAM;

extern int SHIFT_DIGITS;
extern int MAX_KEY_LENGTH;

void randomizeBodies(PhiloxEngine &rng, Body *bodies, int n);
void createOctree(std::vector<DFTNode> &dft, Octree *octree, Body *bodies,
                  int nBodies);
void nbodyIterate(Body *p, float dt, int nBodies, int nIters);
float checkAccuracy(Body *p, Body *orig, int nBodies, int nIters);
void integratePositions(Body *p, float dt, int start, int end);
int bodyForce(Body *p, float dt, int n, std::vector<DFTNode> dft);
