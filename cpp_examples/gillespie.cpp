#include <vector>
#include "../Neuron.hpp"
#include "../compartments/Soma.hpp"
#include "../compartments/Dendritic_segment.hpp"
#include "../compartments/Spine.hpp"
#include "../engines/analytic/Analytic_engine.hpp"
#include "../engines/stochastic/Gillespie_engine.hpp"

using namespace std;

#define N_AVRG  10000
#define N_FORKS 1 // Note that it is (2^N_FORKS - 1)*3 compartments (if 2 synapses on each dend seg)!
void fork_dendrite(Dendritic_segment* ds, size_t depth=0) {
  if (depth < N_FORKS) {
    auto ds1 = new Dendritic_segment(*ds, ds->get_name() + "-1");
    new Spine(*ds1, "s_" + ds1->get_name() + "_1");
    new Spine(*ds1, "s_" + ds1->get_name() + "_2", .6, 6 + .01);
    fork_dendrite(ds1, depth+1);

    auto ds2 = new Dendritic_segment(*ds, ds->get_name() + "-2");
    new Spine(*ds2, "s_" + ds2->get_name() + "_1");
    new Spine(*ds2, "s_" + ds2->get_name() + "_2", .6, 6 + .02);
    fork_dendrite(ds2, depth+1);
  }
}

int main() {

  PRNG rnd(1);

  std::string file_name = "../../data/gillespie/single_Y_fork_test/single_Y_fork_";
  // std::string file_name = "../../data/soma_only/tests2/fast_gene_activation__fast_protein_decay_";//"../../data/soma_only/slow_proteins_1e-2/sim/g_";

  std::list<double> times;
  for(size_t i=0; i<100000; ++i)
    times.push_back(i*.1);

  std::cerr << file_name + ":\n";
  for(size_t i=0; i<N_AVRG; ++i) {
    std::cerr << i << '\t';

    Soma *p_soma = new Soma("soma" /*,Parameters of the soma*/);

    //    p_soma->set_number_of_gene_copies(1).set_gene_activation_rate(1).set_gene_deactivation_rate(0);//.set_transcription_rate(10);

    ///// Branching neuron
    auto *p_ds = new Dendritic_segment(*p_soma, "d_1");
    new Spine(*p_ds, "s_1_1", .6, 6 + .03);
    new Spine(*p_ds, "s_1_2", .6, 6 + .04);

    fork_dendrite(p_ds);

    Neuron *p_neuron = new Neuron(*p_soma, "simple_neuron");
  
    std::ofstream ofs(file_name + std::to_string(i));
    Gillespie_engine(*p_neuron, rnd).run_Gillespie(times, ofs);

    ofs.close();
  }
  
  return 0;
}
