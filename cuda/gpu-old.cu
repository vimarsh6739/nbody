#include "main.h"

#include <cmath>

// Define all CLI options

extern int launchDS(Body *bodies, float dt, int nBodies);
extern int launchDFT(Body *bodies, float dt, int nBodies,
                     std::vector<DFTNode> &dft);
extern int launchMAC(Body *bodies, float dt, int nBodies,
                     std::vector<DFTNode> &dft, float mac_param);

int bodyForce(Body *p, float dt, int n, std::vector<DFTNode> dft) {
  if (USE_TREE && USE_BH)
    launchMAC(p, dt, n, dft, MAC_PARAM);
  else if (USE_TREE)
    launchDFT(p, dt, n, dft);
  else if (CUDA)
    launchDS(p, dt, n);
  else
    allInteractionsDS(p, dt, n);

  return 0;
}
