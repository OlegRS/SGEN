#ifndef __SOM_DEN_JUNCTION_HPP__
#define __SOM_DEN_JUNCTION_HPP__

#include "../Neuron.hpp"

class Neuron::Som_den_junction : public Junction {
  friend class Analytic_engine;

public:
  Som_den_junction();
  Som_den_junction(Compartment* p_from, Compartment* p_to) : Junction(p_from, p_to)
  {set_hopping_rate_constants();}

  Junction& set_hopping_rate_constants();

  // std::vector<double> hopping_rates() {return std::vector<double>{fwd_mRNA_hop_rate, bkwd_mRNA_hop_rate, fwd_prot_hop_rate, bkwd_prot_hop_rate};}

  Junction::Type type() const {return SOM_DEN;}
};

#endif
