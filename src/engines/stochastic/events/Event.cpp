#include "../../../../include/engines/stochastic/events/Event.hpp"

Event::Type::operator std::string() const {
  if(id==GENE_ACTIVATION) return "Gene_activation";
  else if(id==GENE_DEACTIVATION) return "Gene_deactivation";
  else if(id==MRNA_CREATION) return "mRNA_creation";
  else if(id==MRNA_DECAY) return "mRNA_decay";
  else if(id==PROTEIN_CREATION) return "Protein_creation";
  else if(id==PROTEIN_DECAY) return "Protein_decay";
  else if(id==MRNA_HOP_FORWARD) return "mRNA_hop_forward";
  else if(id==MRNA_HOP_BACKWARD) return "mRNA_hop_backward";
  else if(id==PROT_HOP_FORWARD) return "Prot_hop_forward";
  else if(id==PROT_HOP_BACKWARD) return "Prot_hop_backward";
  else return "UNKNOWN-TYPE EVENT";
}
