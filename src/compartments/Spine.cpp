#include "../../include/compartments/Spine.hpp"

Spine::Spine(Compartment &parent, const std::string& name, const double& length, const double& radius, const double& binding_rate, const double& unbinding_rate, const double& d_theta, const double& d_phi) : Compartment(parent, name, length, radius, d_theta, d_phi), protein_binding_rate(binding_rate), protein_unbinding_rate(unbinding_rate) {

  protein_decay_rate = 0;

  parent.p_descendants.push_back(this);
  
  if(parent.p_neuron) {
    p_neuron = parent.p_neuron;
    p_neuron -> p_synapses.push_back(this);
  }
}

Spine::Spine(Compartment &parent, const std::string& name, const double &protein_binding_rate, const double &protein_unbinding_rate, const double &protein_decay_rate, const unsigned int &protein_number) : Compartment(2.5, name), protein_binding_rate(protein_binding_rate), protein_unbinding_rate(protein_unbinding_rate) {
  this->protein_decay_rate = protein_decay_rate;
  
  parent.p_descendants.push_back(this);
  
  if(parent.p_neuron) {
    p_neuron = parent.p_neuron;
    p_neuron -> p_synapses.push_back(this);
  }
}

std::ostream& operator<<(std::ostream& os, const Spine& syn) {
  os << syn.get_name() << ", params: " << "br=" << syn.protein_binding_rate << "; ubr=" << syn.protein_unbinding_rate;
  return os;
}
