#include "../../include/junctions/Junction.hpp"

Junction::Type::operator std::string() const {
  if(id==SOM_DEN) return "SOM_DEN_junction";
  else if(id==DEN_DEN) return "DEN_DEN_junction";
  else if(id==DEN_SYN) return "DEN_SYN_junction";
  else return "UNKNOWN-TYPE JUNCTION";
}
