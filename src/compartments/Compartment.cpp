#include "../../include/compartments/Compartment.hpp"
#include "../../include/compartments/Soma.hpp"
#include "../../include/compartments/Dendritic_segment.hpp"
#include "../../include/compartments/Spine.hpp"
#include "../../include/randomisation/PRNG.hpp"

Compartment::Compartment(const Compartment& parent, const std::string& name, const double& length, const double& radius, const double& d_theta, const double& d_phi, const std::string& placement, const double& mRNA_decay_rate, const double& translation_rate, const double& protein_decay_rate, const double& mRNA_diffusion_constant, const double& protein_diffusion_constant, const double& mRNA_forward_trafficking_velocity, const double& mRNA_backward_trafficking_velocity, const double& protein_forward_trafficking_velocity, const double& protein_backward_trafficking_velocity) : theta(parent.theta+d_theta), phi(parent.phi+d_phi), name(name),length(length), x(parent.x + length*std::sin(theta)*std::cos(phi)), y(parent.y + length*std::sin(theta)*std::sin(phi)), z(parent.z + length*std::cos(theta)), r(radius), placement(placement), mRNA_decay_rate(mRNA_decay_rate), translation_rate(translation_rate), protein_decay_rate(protein_decay_rate), mRNA_diffusion_constant(mRNA_diffusion_constant), protein_diffusion_constant(protein_diffusion_constant), mRNA_forward_trafficking_velocity(mRNA_forward_trafficking_velocity), mRNA_backward_trafficking_velocity(mRNA_backward_trafficking_velocity), protein_forward_trafficking_velocity(protein_forward_trafficking_velocity), protein_backward_trafficking_velocity(protein_backward_trafficking_velocity), protein_creation(this),protein_decay(this),mRNA_creation(this),mRNA_decay(this) {
    if (placement == "end") // Most of the time
      return;

    if(placement ==  "middle") {
      offset = coordinate_offset(parent);
    }
    else if (placement == "random") {
      if (PRNG::instance().rand_max() != 1) {
        std::cerr << "--- WARNING: Changing rand_max from " << PRNG::instance().rand_max() << " to " << "1 for random compartment placement\n";
        PRNG::instance().set_max(1);
      }
      offset = coordinate_offset(parent, PRNG::instance()());
    }
    else
      std::cerr << "--- ERROR: Unknown compartment placement method\n"
                << "- Should be either \"end\", \"middle\" or \"random\".\n";

    x += offset[0];
    y += offset[1];
    z += offset[2];
}


std::vector<double> Compartment::coordinate_offset(const Compartment& parent, const double& length_fraction) {
     double
       par_l = parent.get_length()*length_fraction,
       par_theta = parent.get_theta(),
       par_phi = parent.get_phi();

     std::vector<double>
       orient = {std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta)},
       par_orient = {std::sin(par_theta)*std::cos(par_phi), std::sin(par_theta)*std::sin(par_phi), std::cos(par_theta)};

     double dot_prod = orient[0]*par_orient[0] + orient[1]*par_orient[1] + orient[2]*par_orient[2];
     std::vector<double> dr_vec(3);
     dr_vec[0] = orient[0] - dot_prod*par_orient[0];
     dr_vec[1] = orient[1] - dot_prod*par_orient[1];
     dr_vec[2] = orient[2] - dot_prod*par_orient[2];
     double dr_vec_norm = std::sqrt(dr_vec[0]*dr_vec[0] + dr_vec[1]*dr_vec[1] + dr_vec[2]*dr_vec[2]);
     dr_vec[0] *= parent.radius()/dr_vec_norm;
     dr_vec[1] *= parent.radius()/dr_vec_norm;
     dr_vec[2] *= parent.radius()/dr_vec_norm;
    
     return {dr_vec[0]-par_l*par_orient[0], dr_vec[1]-par_l*par_orient[1], dr_vec[2]-par_l*par_orient[2]};
}

std::ostream& operator<<(std::ostream &os, const Compartment &comp) {
  if(comp.type() == SPINE)
    os << static_cast<const Spine&>(comp);
  else if(comp.type() == APICAL_DENDRITE || comp.type() == BASAL_DENDRITE)
    os << static_cast<const Dendritic_segment&>(comp);
  else if(comp.type() == SOMA)
    os << static_cast<const Soma&>(comp);
  else {
    os << "\n!!!!! ERROR: UNKNOWN TYPE COMPARTMENT DETECTED !!!!!\n";
    exit(1);
  }

  os << " | n_mRNAs = " << comp.n_mRNAs << ", n_proteins = " << comp.n_proteins;
      
  return os;
}

// Compartment::Compartment(Compartment *parent, const std::string& name) : name(name), p_neuron(parent.p_neuron) {
//   // Construct the compartment and link it to its parent in the tree
//   p_attached_compartments.push_back(parent);
//   parent -> p_attached_compartments.push_back(*this);
//   p_neuron = parent->p_neuron;
//   if(p_neuron == nullptr)
//     return;
//   else
//     add_compartment(*this);
// }
