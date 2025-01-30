#include "../../../../include/compartments/Soma.hpp"

void Soma::Gene_activation::operator()() {
  auto& location = *((Soma*)p_location);
  location.n_active_genes++;

  rate -= location.gene_activation_rate;
  location.p_neuron->total_rate -= location.gene_activation_rate;
  
  location.gene_deactivation.rate += location.gene_deactivation_rate;
  location.p_neuron->total_rate += location.gene_deactivation_rate;

  location.mRNA_creation.rate += location.transcription_rate;
  location.p_neuron->total_rate += location.transcription_rate;
}
