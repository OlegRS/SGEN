#ifndef __COMPARTMENT_HPP__
#define __COMPARTMENT_HPP__

#include <iostream>
#include <list>
#include <tuple>
#include "../engines/stochastic/events/Event.hpp"
#include <math.h>

#define PI 3.141592653589793

class Neuron;
class Spine;
class Dendritic_segment;
class Junction;
class Analytic_engine;
class Gillespie_engine;
class Morphologic_engine;

class Compartment {// Abstract class
  friend class Neuron;
  friend class Spine;
  friend class Dendritic_segment;
  friend class Junction;
  friend class Analytic_engine;
  friend class Gillespie_engine;
  friend class Morphologic_engine;
protected:
  
  struct Type {
#define SOMA 1 // Better replace these macros with enum
#define BASAL_DENDRITE 3
#define APICAL_DENDRITE 4
#define SPINE 7

    short id; // Type id consistent with .swc neuron morphology format
  
    Type(const int& type_id) : id(type_id) {}
    operator std::string() const;

    friend bool operator==(const Type& type, const short&& id) {return type.id == id;}
    friend bool operator!=(const Type& type, const short&& id) {return type.id != id;}
    friend std::ostream& operator<<(std::ostream& os, const Type& type) {
      return os << std::string(type);
    }
  };

  std::string name;

  double mRNA_decay_rate = 1.2e-5*3600,
    translation_rate = 0.021*3600,
    protein_decay_rate = 1.21e-6*3600,
    mRNA_diffusion_constant = 3.4e-3, // 3.4e-3 um2/s for CaMKII
    protein_diffusion_constant = .24, //0.24 um2/s for CaMKII
    mRNA_forward_trafficking_velocity = .5e-2, // 4e-2 um/s for CaMKII
    mRNA_backward_trafficking_velocity = .1e-2, // 4e-2 um/s for CaMKII
    protein_forward_trafficking_velocity = 0,
    protein_backward_trafficking_velocity = 0;

  double r=10, length = 200, //micrometers
    theta=0, phi=0, // Orientaiton in radians
    x, y, z; // Coordinates of the first end of the compartment

  std::string placement="end";

  std::list<Compartment*> p_descendants; // Descendants of the compartment in the tree

  std::list<std::list<Junction*>::iterator> it_p_out_junctions, // Outgoing junctions from compartment  
                                            it_p_in_junctions; // Incomming junctions to compartment
  
  Neuron *p_neuron = nullptr; // Pointer to its neuron
  std::list<Compartment*>::iterator iterator; // Its position in its neuron's container

  // For Analytic_engine
  size_t o1_index, mRNA_ind, prot_ind, id;

  // For stochastic engines
  struct Protein_creation : public Event {
    Protein_creation(Compartment* p_loc) : Event(p_loc) {}
    Event::Type type() const override {return Event::Type::PROTEIN_CREATION;}
    void operator()() override;
  } protein_creation;
  struct Protein_decay : public Event {
    Protein_decay(Compartment* p_loc) : Event(p_loc) {}
    Event::Type type() const override {return Event::Type::PROTEIN_DECAY;}
    void operator()() override;
  } protein_decay;
  struct MRNA_creation : public Event {
    MRNA_creation(Compartment* p_loc) : Event(p_loc) {}
    Event::Type type() const override {return Event::Type::MRNA_CREATION;}
    void operator()() override;
  } mRNA_creation;
  struct MRNA_decay : public Event {
    MRNA_decay(Compartment* p_loc) : Event(p_loc) {}
    Event::Type type() const override {return Event::Type::MRNA_DECAY;}
    void operator()() override;
  } mRNA_decay;

  double n_mRNA_expectation=0, n_prot_expectation=0;
  size_t n_mRNAs=0, n_proteins=0;

  Compartment& clear_junctions() {
    it_p_out_junctions.clear();
    it_p_in_junctions.clear();
    return *this;
  }

public:
  Compartment(const std::string& name = "no_name") : name(name),
                                                     protein_creation(this),
                                                     protein_decay(this),
                                                     mRNA_creation(this),
                                                     mRNA_decay(this) {}
  
  Compartment(const double& length=200, const std::string& name = "no_name", const double& x=0, const double& y=0, const double& z=0, const double& radius=10, const double& theta=0, const double& phi=0, const double& mRNA_decay_rate=0.0432, const double& translation_rate=75.6, const double& protein_decay_rate=0.004356, const double& mRNA_diffusion_constant=3.4e-3, const double& protein_diffusion_constant=.24, const double& mRNA_forward_trafficking_velocity=.5e-2, const double& mRNA_backward_trafficking_velocity=.1e-2, const double& protein_forward_trafficking_velocity=0, const double& protein_backward_trafficking_velocity=0) : name(name), length(length), x(x), y(y), z(z), r(radius), theta(theta), phi(phi), mRNA_decay_rate(mRNA_decay_rate), translation_rate(translation_rate), protein_decay_rate(protein_decay_rate), mRNA_diffusion_constant(mRNA_diffusion_constant), protein_diffusion_constant(protein_diffusion_constant), mRNA_forward_trafficking_velocity(mRNA_forward_trafficking_velocity), mRNA_backward_trafficking_velocity(mRNA_backward_trafficking_velocity), protein_forward_trafficking_velocity(protein_forward_trafficking_velocity), protein_backward_trafficking_velocity(protein_backward_trafficking_velocity), protein_creation(this), protein_decay(this), mRNA_creation(this), mRNA_decay(this) {}

  Compartment(const Compartment& parent, const double& x, const double& y, const double& z, const double& radius=10, const std::string& name = "no_name", const double& length=200, const double& d_theta=0, const double& d_phi=0, const double& mRNA_decay_rate=0.0432, const double& translation_rate=75.6, const double& protein_decay_rate=0.004356, const double& mRNA_diffusion_constant=3.4e-3, const double& protein_diffusion_constant=.24, const double& mRNA_forward_trafficking_velocity=.5e-2, const double& mRNA_backward_trafficking_velocity=.1e-2, const double& protein_forward_trafficking_velocity=0, const double& protein_backward_trafficking_velocity=0) : name(name),length(length),protein_creation(this),protein_decay(this),mRNA_creation(this),mRNA_decay(this), x(x), y(y), z(z), r(radius), theta(parent.theta+d_theta), phi(parent.phi+d_phi), mRNA_decay_rate(mRNA_decay_rate), translation_rate(translation_rate), protein_decay_rate(protein_decay_rate), mRNA_diffusion_constant(mRNA_diffusion_constant), protein_diffusion_constant(protein_diffusion_constant), mRNA_forward_trafficking_velocity(mRNA_forward_trafficking_velocity), mRNA_backward_trafficking_velocity(mRNA_backward_trafficking_velocity), protein_forward_trafficking_velocity(protein_forward_trafficking_velocity), protein_backward_trafficking_velocity(protein_backward_trafficking_velocity) {}

  Compartment(const Compartment& parent, const std::string& name = "no_name", const double& length=200, const double& radius=10, const double& d_theta=0, const double& d_phi=0, const std::string& placement="end", const double& mRNA_decay_rate=0.0432, const double& translation_rate=75.6, const double& protein_decay_rate=0.004356, const double& mRNA_diffusion_constant=3.4e-3, const double& protein_diffusion_constant=.24, const double& mRNA_forward_trafficking_velocity=.5e-2, const double& mRNA_backward_trafficking_velocity=.1e-2, const double& protein_forward_trafficking_velocity=0, const double& protein_backward_trafficking_velocity=0) : theta(parent.theta+d_theta), phi(parent.phi+d_phi), name(name),length(length), x(parent.x + length*sin(theta)*cos(phi)), y(parent.y + length*sin(theta)*sin(phi)), z(parent.z + length*cos(theta)), r(radius), placement(placement), mRNA_decay_rate(mRNA_decay_rate), translation_rate(translation_rate), protein_decay_rate(protein_decay_rate), mRNA_diffusion_constant(mRNA_diffusion_constant), protein_diffusion_constant(protein_diffusion_constant), mRNA_forward_trafficking_velocity(mRNA_forward_trafficking_velocity), mRNA_backward_trafficking_velocity(mRNA_backward_trafficking_velocity), protein_forward_trafficking_velocity(protein_forward_trafficking_velocity), protein_backward_trafficking_velocity(protein_backward_trafficking_velocity), protein_creation(this),protein_decay(this),mRNA_creation(this),mRNA_decay(this) {
    if (placement == "end") // Most of the time
      return;

    double par_l,
      par_theta = parent.get_theta(),
      par_phi = parent.get_phi();

    if(placement ==  "middle") 
      par_l = parent.get_length()/2;
    else if (placement == "random")
      par_l = parent.get_length()/2;
    else
      std::cerr << "--- ERROR: Unknown compartment placement method\n"
                << "- Should be either \"end\", \"middle\" or \"random\".\n";
    
    x -= par_l*sin(par_theta)*cos(par_phi);
    y -= par_l*sin(par_theta)*sin(par_phi);
    z -= par_l*cos(par_theta);
  }
  
  virtual Type type() const = 0;
  std::string get_name() const {return name;}

  Neuron* which_neuron() {return p_neuron;}

  double cross_section() {return PI*r*r;}
  
  Compartment& connect_to(Compartment&); // Linking another compartment
  Compartment& disconnect_from(Compartment&); // Unlinking another compartment

  Compartment& set_mRNA_decay_rate(const double& rate) {mRNA_decay_rate=rate; return *this;}
  Compartment& set_translation_rate(const double& rate) {translation_rate=rate; return *this;}
  Compartment& set_protein_decay_rate(const double& rate) {protein_decay_rate=rate; return *this;}
  
  size_t mRNA_count() const {return n_mRNAs;}
  size_t protein_count() const {return n_proteins;}
  double radius() const {return r;}
  Compartment& set_radius(const double& radius) {r=radius; return *this;}
  double get_length() const {return length;}
  double get_theta() const {return theta;}
  double get_phi() const {return phi;}
  std::tuple<double,double,double> position() const {return std::make_tuple(x,y,z);}
  std::tuple<double,double> orientation() const {return std::make_tuple(theta, phi);}
  
  virtual ~Compartment() = default;

  friend std::ostream& operator<<(std::ostream&, const Compartment&);
  friend std::ostream& operator<<(std::ostream&, const Junction&);
};

#endif
