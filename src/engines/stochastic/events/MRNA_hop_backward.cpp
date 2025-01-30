#include "../../../../include/junctions/Junction.hpp"

void Junction::MRNA_hop_backward::operator()() {
  Junction& location = *((Junction*)p_location);
  location.p_from->mRNA_creation();
  location.p_to->mRNA_decay();
}
