#ifndef __DEN_DEN_JUNCTION_HPP__
#define __DEN_DEN_JUNCTION_HPP__

#include "../Neuron.hpp"

class Neuron::Den_den_junction : public Junction {
  friend class Analytic_engine;
  
public:
  Den_den_junction();
  Den_den_junction(Compartment* p_from, Compartment* p_to) : Junction(p_from, p_to)
  { set_hopping_rate_constants(); }
  
  Junction& set_hopping_rate_constants();

  Junction::Type type() const {return DEN_DEN;}
};


#endif
