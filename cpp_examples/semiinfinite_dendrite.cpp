#include <vector>
#include "../Neuron.hpp"
#include "../compartments/Soma.hpp"
#include "../compartments/Dendritic_segment.hpp"
#include "../compartments/Spine.hpp"
#include "../engines/analytic/Analytic_engine.hpp"

using namespace std;

#define N_DS 250
#define LENGTH 5000.

int main() {

  Soma soma("soma", LENGTH/N_DS);
  
  // Constructing a neuron
  Dendritic_segment* p_ds = new Dendritic_segment(soma, "d_1-1", LENGTH/N_DS);
  
  for(size_t i=0; i<N_DS-1; ++i)
    p_ds = new Dendritic_segment(*p_ds, "d_1-" + to_string(i+2), LENGTH/N_DS);

  Neuron neuron(soma, "Test_neuron");
  
  cout << neuron << endl;


  Analytic_engine ae(neuron);
  
  //  ae.sem_stationary_expectations().sem_stationary_pearson_correlations();


  ae.mRNA_stationary_expectations().protein_stationary_expectations();

  ae.gene_mRNA_stationary_covariances().mRNA_mRNA_stationary_covariances();//.gene_protein_stationary_covariances().mRNA_protein_stationary_covariances().protein_protein_stationary_covariances();

  
  return 0;
}
