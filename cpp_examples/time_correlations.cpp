#include <vector>
#include "../Neuron.hpp"
#include "../compartments/Soma.hpp"
#include "../compartments/Dendritic_segment.hpp"
#include "../compartments/Spine.hpp"
#include "../engines/analytic/Analytic_engine.hpp"
#include "../randomisation/PRNG.hpp"

using namespace std;

#define N_FORKS 1 // Note that it is (2^N_FORKS - 1)*3 compartments (if 2 synapses on each dend seg)!

#define syn_dec_rate 0.
// #define syn_dec_rate 1.21e-6*3600
// #define syn_dec_rate 1.21e-3*3600

PRNG rnd(1);

void fork_dendrite(Dendritic_segment* ds, size_t depth=0) {
  if (depth < N_FORKS) {
    auto ds1 = new Dendritic_segment(*ds, ds->get_name() + "-1");
    new Spine(*ds1, "s_" + ds1->get_name() + "_1", .6, 6*(1 + rnd()), syn_dec_rate);
    new Spine(*ds1, "s_" + ds1->get_name() + "_2", .6, 6*(1 + rnd()), syn_dec_rate);
    fork_dendrite(ds1, depth+1);

    auto ds2 = new Dendritic_segment(*ds, ds->get_name() + "-2");
    new Spine(*ds2, "s_" + ds2->get_name() + "_1", .6, 6*(1 + rnd()), syn_dec_rate);
    new Spine(*ds2, "s_" + ds2->get_name() + "_2", .6, 6*(1 + rnd()), syn_dec_rate);
    fork_dendrite(ds2, depth+1);
  }
}

int main() {

  Soma soma("soma", 50 /*,Parameters of the soma*/);
  
  ///// Branching neuron
  Dendritic_segment* p_ds = new Dendritic_segment(soma, "d_1");
  new Spine(*p_ds, "s_1_1", .6, 6 + 6*rnd(), syn_dec_rate);
  new Spine(*p_ds, "s_1_2", .6, 6 + 6*rnd(), syn_dec_rate);
  
  // Dendritic_segment* p_ds_small = new Dendritic_segment(*p_ds, "d_2", 1);
  // fork_dendrite(p_ds_small);

  fork_dendrite(p_ds);

  Neuron neuron(soma, "Test_neuron");

  // cout << "----------------- ANALYTIC ENGINE -----------------\n";
  Analytic_engine ae(neuron);

  std::list<double> times;
  double d_tau = .01;
  for(double tau=0; tau<10000; tau+=d_tau)
    times.push_back(tau);

  ae.sem_stationary_time_correlations(times);
  
  return 0;
}
