#include "../../include/junctions/Som_den_junction.hpp"
#include "../../include/compartments/Dendritic_segment.hpp"
#include "../../include/compartments/Soma.hpp"

Junction& Neuron::Som_den_junction::set_hopping_rate_constants() {

  double min_cs = p_to->cross_section(); //Contact area between two compartments
  if(p_from->cross_section() < min_cs)
    min_cs = p_from->cross_section();

  double
    fwd_diff_factor = min_cs/((p_from->length)*(p_from->length)*p_from->cross_section()),
    bkwd_diff_factor = min_cs/((p_to->length)*(p_to->length)*p_to->cross_section()),
    
    mRNA_fwd_diff_rate = p_from->mRNA_diffusion_constant * fwd_diff_factor,
    mRNA_bkwd_diff_rate = p_to->mRNA_diffusion_constant * bkwd_diff_factor,    
    mRNA_fwd_traff_rate = p_from->mRNA_forward_trafficking_velocity / p_from->length,
    mRNA_bkwd_traff_rate = p_to->mRNA_backward_trafficking_velocity / p_to->length,    
    prot_fwd_diff_rate = p_from->protein_diffusion_constant * fwd_diff_factor,
    prot_bkwd_diff_rate = p_to->protein_diffusion_constant * bkwd_diff_factor,    
    prot_fwd_traff_rate = p_from->protein_forward_trafficking_velocity / p_from->length,
    prot_bkwd_traff_rate = p_to->protein_backward_trafficking_velocity / p_to->length;
    
  fwd_mRNA_hop_rate = mRNA_fwd_diff_rate + mRNA_fwd_traff_rate;
  bkwd_mRNA_hop_rate = mRNA_bkwd_diff_rate + mRNA_bkwd_traff_rate;

  fwd_prot_hop_rate = prot_fwd_diff_rate + prot_fwd_traff_rate;
  bkwd_prot_hop_rate = prot_bkwd_diff_rate + prot_bkwd_traff_rate;
  
  return *this;
}
