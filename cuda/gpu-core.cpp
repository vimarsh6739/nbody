
#include "main.h"

#include <cmath>

// Define all CLI options

extern int launchDS(Body *bodies, float dt, int nBodies);

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
