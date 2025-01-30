#ifndef __PRNG_HPP__
#define __PRNG_HPP__

#include<random>

typedef std::mt19937_64 PRN_gen; // Mersenne twister PRNG

class PRNG {
  PRN_gen gen;
  std::uniform_real_distribution<double> dis;

public:  
  PRNG(const double& MAX_RAND): gen(time(NULL)), dis(0,MAX_RAND) {}
  PRNG(const double& MAX_RAND, const int seed): gen(seed), dis(0,MAX_RAND) {}
  double rand_max() {return dis.max();}
  double operator()() { return dis(gen); }
};

#endif
