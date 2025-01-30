#ifndef __DEN_SYN_JUNCTION_HPP__
#define __DEN_SYN_JUNCTION_HPP__

#include "../Neuron.hpp"
#include "../compartments/Spine.hpp"

class Neuron::Den_syn_junction : public Junction {
  friend class Analytic_engine;
  
public:
  Den_syn_junction();
  Den_syn_junction(Compartment* p_from, Compartment* p_to) : Junction(p_from, p_to) {
    set_hopping_rate_constants(); }

  Junction& set_hopping_rate_constants();

  Junction::Type type() const {return DEN_SYN;}
};

#endif
