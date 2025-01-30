#ifndef __EVENT_HPP__
#define __EVENT_HPP__

#include <ostream>

struct Event {
  
  struct Type {
#define GENE_ACTIVATION 0
#define GENE_DEACTIVATION 1
#define MRNA_CREATION 2
#define MRNA_DECAY 3
#define PROTEIN_CREATION 4
#define PROTEIN_DECAY 5 
#define MRNA_HOP_FORWARD 6
#define MRNA_HOP_BACKWARD 7
#define PROT_HOP_FORWARD 8
#define PROT_HOP_BACKWARD 9
    
    short id;
  
    Type(const int& type_id) : id(type_id) {}
    operator std::string() const;

    friend bool operator==(const Type& type, const short&& id) {return type.id == id;}
    friend bool operator!=(const Type& type, const short&& id) {return type.id != id;}
    friend std::ostream& operator<<(std::ostream& os, const Type& type) {
      return os << std::string(type);
    }
  };


  Event(void* p_loc=NULL) : p_location(p_loc) {}
  
  void* p_location; // Compartment or junction, where event takes place
  double rate;
  Event& set_rate(double r) {rate=r; return *this;}
  virtual void operator()() = 0;
  virtual Type type() = 0;
};

#endif
