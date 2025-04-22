#include "../../include/junctions/Den_den_junction.hpp"
#include "../../include/compartments/Dendritic_segment.hpp"

Junction& Neuron::Den_den_junction::set_hopping_rate_constants() {

  double
    diffusion_scaling_factor = ( 1/((p_from->length)*(p_from->length)) + 1/((p_to->length)*(p_to->length)) )/2.,
    trafficking_scaling_factor = ( 1/(p_from->length) + 1/(p_to->length) )/2.; //p_from->cross_section()

  double
    mRNA_fwd_diff_rate = p_from->mRNA_diffusion_constant * diffusion_scaling_factor,
    mRNA_bkwd_diff_rate = p_to->mRNA_diffusion_constant * diffusion_scaling_factor,
    mRNA_fwd_traff_rate = p_from->mRNA_forward_trafficking_velocity * trafficking_scaling_factor,
    mRNA_bkwd_traff_rate = p_to->mRNA_backward_trafficking_velocity * trafficking_scaling_factor,
    prot_fwd_diff_rate = p_from->protein_diffusion_constant * diffusion_scaling_factor,
    prot_bkwd_diff_rate = p_to->protein_diffusion_constant * diffusion_scaling_factor,
    prot_fwd_traff_rate = p_from->protein_forward_trafficking_velocity * trafficking_scaling_factor,
    prot_bkwd_traff_rate = p_to->protein_backward_trafficking_velocity * trafficking_scaling_factor;
    
  fwd_mRNA_hop_rate = (mRNA_fwd_diff_rate + mRNA_fwd_traff_rate);///static_cast<Dendritic_segment*>(p_from)->n_descending_DS;
  bkwd_mRNA_hop_rate = mRNA_bkwd_diff_rate + mRNA_bkwd_traff_rate;

  fwd_prot_hop_rate = (prot_fwd_diff_rate + prot_fwd_traff_rate);///static_cast<Dendritic_segment*>(p_from)->n_descending_DS;
  bkwd_prot_hop_rate = prot_bkwd_diff_rate + prot_bkwd_traff_rate;

  return *this;
}
