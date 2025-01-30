#include "../../../../include/junctions/Junction.hpp"

void Junction::MRNA_hop_forward::operator()() {
  Junction& location = *((Junction*)p_location);
  location.p_to->mRNA_creation();
  location.p_from->mRNA_decay();
}
