#include "../../../../include/Neuron.hpp"

void Compartment::MRNA_creation::operator()() {
  Compartment& location = *((Compartment*)p_location);
  location.n_mRNAs++;
  // Updating mRNA decay rate
  location.mRNA_decay.rate += location.mRNA_decay_rate;
  location.p_neuron->total_rate += location.mRNA_decay_rate;
  // Updating protein creation rate
  location.protein_creation.rate += location.translation_rate;
  location.p_neuron->total_rate += location.translation_rate;
  // Updating hopping rates
  for(auto& p_junc : location.it_p_in_junctions) {
    (*p_junc)->mRNA_hop_backward.rate += (*p_junc)->bkwd_mRNA_hop_rate;
    location.p_neuron->total_rate += (*p_junc)->bkwd_mRNA_hop_rate;
  }
  for(auto& p_junc : location.it_p_out_junctions) {
    (*p_junc)->mRNA_hop_forward.rate += (*p_junc)->fwd_mRNA_hop_rate;
    location.p_neuron->total_rate += (*p_junc)->fwd_mRNA_hop_rate;
  }
}
