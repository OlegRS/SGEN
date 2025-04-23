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

PRNG::instance().set_seed(1);

std::string file_name = "../../data/gillespie/all_stationary_time_correlations/stationary_time_correlations_24/TC_";


void fork_dendrite(Dendritic_segment* ds, size_t depth=0) {
  if (depth < N_FORKS) {
    auto ds1 = new Dendritic_segment(*ds, ds->get_name() + "-1");
    new Spine(*ds1, "s_" + ds1->get_name() + "_1", .6, 6 + 6*(.5-PRNG::instance()()));
    new Spine(*ds1, "s_" + ds1->get_name() + "_2", .6, 6 + 6*(.5-PRNG::instance()()));
    fork_dendrite(ds1, depth+1);

    auto ds2 = new Dendritic_segment(*ds, ds->get_name() + "-2");
    new Spine(*ds2, "s_" + ds2->get_name() + "_1", .6, 6 + 6*(.5-PRNG::instance()()));
    new Spine(*ds2, "s_" + ds2->get_name() + "_2", .6, 6 + 6*(.5-PRNG::instance()()));
    fork_dendrite(ds2, depth+1);
  }
}

int main() {

  PRNG rnd_MC(1);

  Soma soma("soma" /*,Parameters of the soma*/);
  
  ///// Branching neuron
  Dendritic_segment* p_ds = new Dendritic_segment(soma, "d_1");
  new Spine(*p_ds, "s_1_1", .6, 6 + 6*(.5-PRNG::instance()()));
  new Spine(*p_ds, "s_1_2", .6, 6 + 6*(.5-PRNG::instance()()));
  fork_dendrite(p_ds);

  Neuron neuron(soma, "Test_neuron");

  std::cout << "Neuron:\n"
            << neuron << std::endl;

  std::list<double> times;
  double d_tau = .01;
  for(double tau=0; tau<10000; tau+=d_tau)
    times.push_back(tau);

  std::cout << "----------------- ANALYTIC ENGINE -----------------\n";
  Analytic_engine ae(neuron);

  ae.sem_stationary_time_correlations(times);


  //////////// Gillespie simulation /////////////
  std::cout << "----------------- GILLESPIE ENGINE -----------------\n";
  for(size_t i=0; i<N_AVRG; ++i) {
    std::cerr << "Writing Gillespie results to: " << file_name + std::to_string(i) << '\n';
    
    std::ofstream ofs_gillespie(file_name + std::to_string(i));
    
    std::cout << neuron << std::endl;
 
    Gillespie_engine(neuron, rnd_MC).run_Gillespie(times, ofs_gillespie);

    ofs_gillespie.close();
  }
  
  return 0;
}
