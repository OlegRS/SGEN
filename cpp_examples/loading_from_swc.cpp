#include "../Neuron.hpp"
#include "../engines/analytic/Analytic_engine.hpp"

using namespace std;

int main() {

  // Neuron neuron("../../data/morphologies/10_2REDO-850-GM18-Ctl-Ctl-Chow-BNL16A-CA1_Finished2g.swc", "laded_neuron");
  Neuron neuron("../../data/morphologies/10_2REDO-850-GM18-Ctl-Ctl-Chow-BNL16A-CA1_Finished2e.swc", "laded_neuron");
  // Neuron neuron("../../data/morphologies/test.swc", "laded_neuron");

  Analytic_engine ae(neuron);

  ae.mRNA_stationary_expectations().protein_stationary_expectations();
  
  ae.gene_mRNA_stationary_covariances().gene_protein_stationary_covariances().mRNA_mRNA_stationary_covariances().mRNA_protein_stationary_covariances().protein_protein_stationary_covariances();
  
  return 0;
}
