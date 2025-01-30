#ifndef __MORPHOLOGIC_ENGINE_HPP__
#define __MORPHOLOGIC_ENGINE_HPP__

#include "../../Neuron.hpp"
#include "../../compartments/Soma.hpp"
#include "../../compartments/Dendritic_segment.hpp"
#include "../../compartments/Spine.hpp"
#include "../../junctions/Som_den_junction.hpp"
#include "../../junctions/Den_syn_junction.hpp"
#include "../../junctions/Den_den_junction.hpp"

class Morphologic_engine {

  Neuron* p_neuron = NULL;
  
public:

  Morphologic_engine(Neuron& neuron) : p_neuron(&neuron) {}

  Morphologic_engine& rediscretise();

  std::vector<std::vector<std::vector<double>>> segments();

  std::vector<double> volumes();
    
  ~Morphologic_engine() {}
};

#endif
