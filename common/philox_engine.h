// philox_engine.h
#ifndef PHILOX_ENGINE_H
#define PHILOX_ENGINE_H

#include "philox_rng.h"
#include <cstdint>

struct PhiloxEngine {
  using result_type = uint32_t;
  explicit PhiloxEngine(uint32_t seed) { 
    philox_seed(seed); 
  }
  uint32_t operator()() { 
    return philox_random_uint32(); 
  }
  static constexpr uint32_t min() { return 0; }
  static constexpr uint32_t max() { return UINT32_MAX; }
};

#endif // PHILOX_ENGINE_H

