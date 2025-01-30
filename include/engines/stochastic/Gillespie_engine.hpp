#ifndef __GILLESPIE_ENGINE_HPP__
#define __GILLESPIE_ENGINE_HPP__

#include "events/Event.hpp"
#include "../../Neuron.hpp"
#include "../../compartments/Soma.hpp"
#include "../../randomisation/PRNG.hpp"

class Gillespie_engine {
  // Parameters
  Neuron* p_neuron = NULL;
  unsigned int dim;

  Compartment* initialise_soma();
  void initialise_from(Compartment&);
  inline double draw_delta_t();
  void update_Gillespie();

  std::vector<Event*> p_events;
  size_t ev_ind = 0; // Event index (needed for recursions)
  PRNG& rnd;
  
  
public:
  // Gillespie_engine(Neuron& neuron) :
  //   p_neuron(&neuron),
  //   dim(3 + 2*neuron.p_dend_segments.size() + neuron.p_synapses.size()),
  //   p_events(6 + 4*neuron.p_dend_segments.size() + 3*(neuron.n_SDJ + neuron.n_DDJ) + 2*neuron.n_DSJ),
  //   rnd(1)
  // {}

  Gillespie_engine(Neuron& neuron, PRNG& prng) :
    p_neuron(&neuron),
    dim(3 + 2*neuron.p_dend_segments.size() + neuron.p_synapses.size()),
    p_events(6 + 3*neuron.p_dend_segments.size() + neuron.p_synapses.size() + 4*(neuron.n_SDJ + neuron.n_DDJ) + 2*neuron.n_DSJ),
    rnd(prng)
  { initialise_from(*initialise_soma()); }

  Gillespie_engine& refresh() { initialise_from(*initialise_soma()); return *this;}

  Gillespie_engine& run_Gillespie(const double& time);
  Gillespie_engine& run_Gillespie(const std::list<double>& write_times, std::ofstream&, const double& time_offset=0);
  
};

#endif
