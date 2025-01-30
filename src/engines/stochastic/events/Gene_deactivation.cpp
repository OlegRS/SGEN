#include "../../../../include/compartments/Soma.hpp"

void Soma::Gene_deactivation::operator()() {
  auto& location = *((Soma*)p_location);
  location.n_active_genes--;
  
  location.gene_activation.rate += location.gene_activation_rate;
  location.p_neuron->total_rate += location.gene_activation_rate;

  rate -= location.gene_deactivation_rate;
  location.p_neuron->total_rate -= location.gene_deactivation_rate;

  location.mRNA_creation.rate -= location.transcription_rate;
  location.p_neuron->total_rate -= location.transcription_rate;
}
