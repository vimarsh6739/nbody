#include "ctimer.h"
#include "main.h"
#include <assert.h>

#ifdef ENABLE_CILK
#include <cilk/cilk.h>
#endif
#include <random>

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

int main(const int argc, const char **argv) {
  return libMain(argc, argv);
}
