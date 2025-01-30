#include "../../../../include/junctions/Junction.hpp"

void Junction::Prot_hop_backward::operator()() {
  auto& location = *((Junction*)p_location);
  location.p_to->protein_decay();
  location.p_from->protein_creation();  
}
