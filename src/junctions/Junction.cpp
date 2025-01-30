#include "../../include/junctions/Junction.hpp"

std::ostream& operator<<(std::ostream &os, const Junction &junc) {
  os << junc.p_from->name << "--[" << junc.type() << "]-->" << junc.p_to->name
     << ": Forward_mRNA_hopping_rate=" << junc.fwd_mRNA_hop_rate
     << ", Backward_mRNA_hopping_rate=" << junc.bkwd_mRNA_hop_rate
     << ", Forward_protein_hopping_rate=" << junc.fwd_prot_hop_rate
     << ", Backward_protein_hopping_rate=" << junc.bkwd_prot_hop_rate;
  return os;
}
