#ifndef __SPINE_HPP__
#define __SPINE_HPP__

#include "../Neuron.hpp"

class Spine : public Compartment {
  friend class Neuron::Den_syn_junction;
  friend class Analytic_engine;
  
  //Parameters
  double protein_binding_rate = .6; // per hour
  double protein_unbinding_rate = 6;

  double length = 200/60,//2.5, //um
    protein_diffusion_constant = .24, //0.24 um2/s for CaMKII
    protein_forward_trafficking_velocity = 0, 
    protein_backward_trafficking_velocity = 0;
  
public:
  Spine();
  Spine(Compartment &parent, const std::string& name = "no_name", const double& length=2, const double& radius=1, const double& binding_rate=.6, const double& unbinding_rate=6, const double& d_theta=PI/2, const double& d_phi=0, const std::string& placement="end");
  Spine(Compartment &parent, const std::string& name, const double &protein_binding_rate, const double &protein_unbinding_rate, const double &protein_decay_rate=0, const unsigned int &protein_number = 0);

  const double& get_protein_binding_rate() const {return protein_binding_rate;}
  const double& get_protein_unbinding_rate() const {return protein_unbinding_rate;}

  Spine& set_protein_binding_rate(const double& rate) {protein_binding_rate=rate; return *this;}
  Spine& set_protein_unbinding_rate(const double& rate) {protein_unbinding_rate=rate; return *this;}

  Compartment::Type type() const override {return SPINE;}

  friend std::ostream& operator<<(std::ostream&, const Spine&);
};

#endif
