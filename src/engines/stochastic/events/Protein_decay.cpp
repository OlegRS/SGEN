#include "../../../../include/Neuron.hpp"

void Compartment::Protein_decay::operator()() {
  auto& location = *((Compartment*)p_location);
  location.n_proteins--;
  // Decrementing protein decay rate
  rate -= location.protein_decay_rate;
  location.p_neuron->total_rate -= location.protein_decay_rate;

  for(auto& p_junc : location.it_p_in_junctions) {
    (*p_junc)->prot_hop_backward.rate -= (*p_junc)->bkwd_prot_hop_rate;
    location.p_neuron->total_rate -= (*p_junc)->bkwd_prot_hop_rate;
  }
  for(auto& p_junc : location.it_p_out_junctions) {
    (*p_junc)->prot_hop_forward.rate -= (*p_junc)->fwd_prot_hop_rate;
    location.p_neuron->total_rate -= (*p_junc)->fwd_prot_hop_rate;
  }
}
