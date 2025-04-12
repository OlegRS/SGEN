#ifndef __DENDRITIC_SEGMENT_HPP__
#define __DENDRITIC_SEGMENT_HPP__

#include "../Neuron.hpp"

class Dendritic_segment : public Compartment {
  friend class Neuron::Den_den_junction;
  friend class Neuron::Som_den_junction;
  friend class Analytic_engine;

  size_t n_descending_DS = 0;// Number of decstnding dendritic segments 
  
  //Parameters
  
public:
  Dendritic_segment(const double &mRNA_diffusion_rate, const double &mRNA_forward_trafficking_rate, const double &protein_diffusion_rate, const double &protein_trafficking_rate, const double &translation_rate, const unsigned int &mRNA_numer = 0, const unsigned int &protein_number = 0, const std::string& name = "no_name");

  Dendritic_segment(Compartment& parent, const double& x, const double& y, const double& z, const double& radius=10, const std::string& name = "no_name", const double& length=200, const double& d_theta=0, const double& d_phi=0, const double& mRNA_decay_rate=0.0432, const double& translation_rate=75.6, const double& protein_decay_rate=0.004356, const double& mRNA_diffusion_constant=3.4e-3, const double& protein_diffusion_constant=.24, const double& mRNA_forward_trafficking_velocity=.5e-2, const double& mRNA_backward_trafficking_velocity=.1e-2, const double& protein_forward_trafficking_velocity=0, const double& protein_backward_trafficking_velocity=0);
  
  Dendritic_segment(Compartment &parent, const std::string& name = "no_name", const double& length=200, const double& radius=10, const double& d_theta=0, const double& d_phi=0, const double& mRNA_decay_rate=0.0432, const double& translation_rate=75.6, const double& protein_decay_rate=0.004356, const double& mRNA_diffusion_constant=3.4e-3, const double& protein_diffusion_constant=.24, const double& mRNA_forward_trafficking_velocity=.5e-2, const double& mRNA_backward_trafficking_velocity=.1e-2, const double& protein_forward_trafficking_velocity=0, const double& protein_backward_trafficking_velocity=0, const bool& MIDDLE_PLACEMENT=false);

  Dendritic_segment& set_translation_rate(const double& tr=75.6) {
    Compartment::translation_rate = tr;
    return *this;
  }

  double translation_rate() {return Compartment::translation_rate;}

  Compartment::Type type() const {return APICAL_DENDRITE;}

  friend std::ostream& operator<<(std::ostream&, const Dendritic_segment&);
};

#endif
