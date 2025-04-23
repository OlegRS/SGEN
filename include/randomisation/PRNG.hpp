#ifndef __PRNG_HPP__
#define __PRNG_HPP__

#include<random>

typedef std::mt19937_64 PRN_gen; // Mersenne twister

class PRNG {
  PRN_gen gen;
  std::uniform_real_distribution<double> dis;

  // Private constructor (singleton style)
  PRNG(const double& max_rand = 1.0, const unsigned& seed = static_cast<unsigned>(time(nullptr)))
    : gen(seed), dis(0, max_rand) {}

public:
  // Delete copy and move constructors to enforce singleton
  PRNG(const PRNG&) = delete;
  PRNG& operator=(const PRNG&) = delete;
  PRNG(PRNG&&) = delete;
  PRNG& operator=(PRNG&&) = delete;

  static PRNG& instance() {
    static PRNG singleton_instance;
    return singleton_instance;
  }

  static PRNG& set_seed(const unsigned& seed) {
    auto& inst = instance();
    inst.gen.seed(seed);
    return inst;
  }

  static PRNG& set_max(const double& max) {
    auto& inst = instance();
    inst.dis = std::uniform_real_distribution<double>(0, max);
    return inst;
  }
  
  double rand_max() {return dis.max();}
  double operator()() { return dis(gen); }
};

#endif
