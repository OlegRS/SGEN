#include "../../include/compartments/Dendritic_segment.hpp"
#include "../../include/compartments/Soma.hpp"

Dendritic_segment::Dendritic_segment(Compartment& parent, const double& x, const double& y, const double& z, const double& radius, const std::string& name, const double& length, const double& d_theta, const double& d_phi, const double& mRNA_decay_rate, const double& translation_rate, const double& protein_decay_rate, const double& mRNA_diffusion_constant, const double& protein_diffusion_constant, const double& mRNA_forward_trafficking_velocity, const double& mRNA_backward_trafficking_velocity, const double& protein_forward_trafficking_velocity, const double& protein_backward_trafficking_velocity) : Compartment(parent, x, y, z, radius, name, length, d_theta, d_phi, mRNA_decay_rate, translation_rate, protein_decay_rate, mRNA_diffusion_constant, protein_diffusion_constant, mRNA_forward_trafficking_velocity, mRNA_backward_trafficking_velocity, protein_forward_trafficking_velocity, protein_backward_trafficking_velocity) {
  
  parent.p_descendants.push_back(this);

  if(parent.type() == BASAL_DENDRITE || parent.type() == APICAL_DENDRITE)
    ++static_cast<Dendritic_segment&>(parent).n_descending_DS;
  else if(parent.type() == SOMA)
    ++static_cast<Soma&>(parent).n_descending_DS;
  
  if(parent.p_neuron) {
    p_neuron = parent.p_neuron;
    p_neuron -> p_dend_segments.push_back(this);
  }
}

Dendritic_segment::Dendritic_segment(Compartment& parent, const std::string& name, const double& length, const double& radius, const double& d_theta, const double& d_phi, const double& mRNA_decay_rate, const double& translation_rate, const double& protein_decay_rate, const double& mRNA_diffusion_constant, const double& protein_diffusion_constant, const double& mRNA_forward_trafficking_velocity, const double& mRNA_backward_trafficking_velocity, const double& protein_forward_trafficking_velocity, const double& protein_backward_trafficking_velocity, const std::string& placement) : Compartment(parent, name, length, radius, d_theta, d_phi, placement, mRNA_decay_rate, translation_rate, protein_decay_rate, mRNA_diffusion_constant, protein_diffusion_constant, mRNA_forward_trafficking_velocity, mRNA_backward_trafficking_velocity, protein_forward_trafficking_velocity, protein_backward_trafficking_velocity) {
  
  parent.p_descendants.push_back(this);

  if(parent.type() == BASAL_DENDRITE || parent.type() == APICAL_DENDRITE)
    ++static_cast<Dendritic_segment&>(parent).n_descending_DS;
  else if(parent.type() == SOMA)
    ++static_cast<Soma&>(parent).n_descending_DS;
  
  if(parent.p_neuron) {
    p_neuron = parent.p_neuron;
    p_neuron -> p_dend_segments.push_back(this);
  }
}

std::ostream& operator<<(std::ostream& os, const Dendritic_segment& ds) {
  os << ds.get_name() <<", n_descending_DS: "<< ds.n_descending_DS;
  os << "\nIn_junctions:\n";
  for(auto& it_p_in_junc : ds.it_p_in_junctions)
    os << **it_p_in_junc << std::endl;
  os << "\nOut_junctions:\n";
  for(auto& it_p_out_junc : ds.it_p_out_junctions)
    os << **it_p_out_junc << std::endl;
  return os;
}
