#include "../../../../include/junctions/Junction.hpp"

void Junction::Prot_hop_forward::operator()() {
  auto& location = *((Junction*)p_location);
  location.p_to->protein_creation();
  location.p_from->protein_decay();
}
