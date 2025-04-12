#ifndef __EVENT_HPP__
#define __EVENT_HPP__

#include <ostream>

struct Event {
  
  class Type {
    short id;

  public:
    enum EventType : short {
       GENE_ACTIVATION = 0,
       GENE_DEACTIVATION,
       MRNA_CREATION,
       MRNA_DECAY,
       PROTEIN_CREATION,
       PROTEIN_DECAY,
       MRNA_HOP_FORWARD,
       MRNA_HOP_BACKWARD,
       PROT_HOP_FORWARD,
       PROT_HOP_BACKWARD
    };
  
    Type(EventType type_id) : id(type_id) {}
    operator std::string() const;

    friend bool operator==(const Type& type, const short& id) {return type.id == id;}
    friend bool operator!=(const Type& type, const short& id) {return type.id != id;}
    friend std::ostream& operator<<(std::ostream& os, const Type& type) {
      return os << std::string(type);
    }
  };


  Event(void* p_loc=nullptr) : p_location(p_loc) {}
  
  void* p_location; // Compartment or junction, where event takes place
  double rate;
  Event& set_rate(double r) {rate=r; return *this;}
  virtual void operator()() = 0;
  virtual Type type() const = 0;

  virtual ~Event() = default;
};

#endif
