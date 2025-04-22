#include "../../include/junctions/Den_syn_junction.hpp"

Junction& Neuron::Den_syn_junction::set_hopping_rate_constants() {
  fwd_prot_hop_rate = static_cast<Spine*>(p_to)->protein_binding_rate;
  bkwd_prot_hop_rate = static_cast<Spine*>(p_to)->protein_unbinding_rate;

  fwd_mRNA_hop_rate = 0;
  bkwd_mRNA_hop_rate = 0;

  return *this;
}
