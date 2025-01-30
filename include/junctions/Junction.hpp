#ifndef __JUNCTION_HPP__
#define __JUNCTION_HPP__

#include <iostream>
#include <vector>
#include "../compartments/Compartment.hpp"

class Analytic_engine;

class Junction {// Abstract class
  friend class Compartment;
  friend class Analytic_engine;
  friend class Gillespie_engine;
  friend class Morphologic_engine;
protected:

  struct Type {
#define SOM_DEN 0 // Better replace these macros with enum
#define DEN_DEN 1
#define DEN_SYN 2

    short id;
  
    Type(const int& type_id) : id(type_id) {}
    operator std::string() const;

    friend bool operator==(const Type& type, const short&& id) {return type.id == id;}
    friend bool operator!=(const Type& type, const short&& id) {return type.id != id;}
    friend std::ostream& operator<<(std::ostream& os, const Type& type) {
      return os << std::string(type);
    }
  };

  Compartment *p_from = NULL, *p_to = NULL;

  // For stochastic engines
  struct MRNA_hop_forward : public Event {
    MRNA_hop_forward(Junction* p_loc) : Event(p_loc) {}
    Event::Type type() {return MRNA_HOP_FORWARD;}
    void operator()();
  } mRNA_hop_forward;
  
  struct MRNA_hop_backward : public Event {
    MRNA_hop_backward(Junction* p_loc) : Event(p_loc) {}
    Event::Type type() {return MRNA_HOP_BACKWARD;}
    void operator()();
  } mRNA_hop_backward;

  struct Prot_hop_forward : public Event {
    Prot_hop_forward(Junction* p_loc) : Event(p_loc) {}
    Event::Type type() {return PROT_HOP_FORWARD;}
    void operator()();
  } prot_hop_forward;

  struct Prot_hop_backward : public Event {
    Prot_hop_backward(Junction* p_loc) : Event(p_loc) {}
    Event::Type type() {return PROT_HOP_BACKWARD;}
    void operator()();
  } prot_hop_backward;

  double fwd_mRNA_hop_rate, bkwd_mRNA_hop_rate, fwd_prot_hop_rate, bkwd_prot_hop_rate;

public:
  
  Junction(Compartment* p_from, Compartment* p_to) : p_from(p_from), p_to(p_to),
                                                     mRNA_hop_forward(this),
                                                     mRNA_hop_backward(this),
                                                     prot_hop_forward(this),
                                                     prot_hop_backward(this) {}
  
  virtual Type type() const = 0;

  virtual Junction& set_hopping_rate_constants() = 0;

  // virtual std::vector<double> hopping_rates() = 0;

  virtual ~Junction() = default;

  friend std::ostream& operator<<(std::ostream&, const Junction&);
};


#endif
