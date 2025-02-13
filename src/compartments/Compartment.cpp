#include "../../include/compartments/Compartment.hpp"
#include "../../include/compartments/Soma.hpp"
#include "../../include/compartments/Dendritic_segment.hpp"
#include "../../include/compartments/Spine.hpp"

// Compartment::Compartment(const Compartment &comp) : name(comp.name),
//                                                     p_descendants(comp.p_descendants),
//                                                     p_neuron(comp.p_neuron),
//                                                     iterator(comp.iterator) {}

std::ostream& operator<<(std::ostream &os, const Compartment &comp) {
  if(comp.type() == SPINE)
    os << static_cast<const Spine&>(comp);
  else if(comp.type() == APICAL_DENDRITE || comp.type() == BASAL_DENDRITE)
    os << static_cast<const Dendritic_segment&>(comp);
  else if(comp.type() == SOMA)
    os << static_cast<const Soma&>(comp);
  else {
    os << "\n!!!!! ERROR: UNKNOWN TYPE COMPARTMENT DETECTED !!!!!\n";
    exit(1);
  }

  os << " | n_mRNAs = " << comp.n_mRNAs << ", n_proteins = " << comp.n_proteins;
      
  return os;
}

// Compartment::Compartment(Compartment *parent, const std::string& name) : name(name), p_neuron(parent.p_neuron) {
//   // Construct the compartment and link it to its parent in the tree
//   p_attached_compartments.push_back(parent);
//   parent -> p_attached_compartments.push_back(*this);
//   p_neuron = parent->p_neuron;
//   if(p_neuron == nullptr)
//     return;
//   else
//     add_compartment(*this);
// }
