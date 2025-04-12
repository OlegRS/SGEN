#include "../../../../include/engines/stochastic/events/Event.hpp"

Event::Type::operator std::string() const {
  switch (static_cast<Event::Type::EventType>(id)) {
  case GENE_ACTIVATION: return "GENE_ACTIVATION";
  case GENE_DEACTIVATION: return "GENE_DEACTIVATION";
  case MRNA_CREATION: return "MRNA_CREATION";
  case MRNA_DECAY: return "MRNA_DECAY";
  case PROTEIN_CREATION: return "PROTEIN_CREATION";
  case PROTEIN_DECAY: return "PROTEIN_DECAY";
  case MRNA_HOP_FORWARD: return "MRNA_HOP_FORWARD";
  case MRNA_HOP_BACKWARD: return "MRNA_HOP_BACKWARD";
  case PROT_HOP_FORWARD: return "PROT_HOP_FORWARD";
  case PROT_HOP_BACKWARD: return "PROT_HOP_BACKWARD";
  default: return "UNKNOWN_TYPE_EVENT";
  }
}
