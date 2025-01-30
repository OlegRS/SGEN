#include <vector>
#include "../Neuron.hpp"
#include "../compartments/Soma.hpp"
#include "../compartments/Dendritic_segment.hpp"
#include "../compartments/Spine.hpp"
#include "../engines/analytic/Analytic_engine.hpp"

using namespace std;

#define N_DS 30
#define DS_LENGTH 200.

int main() {

  Soma soma("soma", DS_LENGTH/N_DS);
  
  // Constructing a neuron
  Dendritic_segment* p_ds = new Dendritic_segment(soma, "d_1-1", DS_LENGTH/N_DS);
  
  for(unsigned int i=0; i<N_DS-1; ++i) {
    p_ds = new Dendritic_segment(*p_ds, "d_1-" + to_string(i+2), DS_LENGTH/N_DS);
    if(i==N_DS/3)
      new Spine(*p_ds, "s_1_1");
    if(i==2*N_DS/3)
      new Spine(*p_ds, "s_1_2");
  }

  Dendritic_segment* p_ds_branching = p_ds;// = new Dendritic_segment(*p_ds, "d_BP", DS_LENGTH/10000.);

  for(unsigned int i=0; i<N_DS; ++i) {
    p_ds = new Dendritic_segment(*p_ds, "d_2-" + to_string(i+1), DS_LENGTH/N_DS);
    if(i==N_DS/3)
      new Spine(*p_ds, "s_1-1_1");
    if(i==2*N_DS/3)
      new Spine(*p_ds, "s_1-1_2");
  }

  p_ds = p_ds_branching;
  
  for(unsigned int i=0; i<N_DS; ++i) {
    p_ds = new Dendritic_segment(*p_ds, "d_3-" + to_string(i+1), DS_LENGTH/N_DS);
    if(i==N_DS/3)
      new Spine(*p_ds, "s_1-2_1");
    if(i==2*N_DS/3)
      new Spine(*p_ds, "s_1-2_2");
  }



  Neuron neuron(soma, "Test_neuron");
  
  cout << neuron << endl;


  Analytic_engine ae(neuron);
  
  //  ae.sem_stationary_expectations().sem_stationary_pearson_correlations();


  ae.mRNA_stationary_expectations().protein_stationary_expectations();

  ae.gene_mRNA_stationary_covariances().gene_protein_stationary_covariances().mRNA_mRNA_stationary_covariances().mRNA_protein_stationary_covariances().protein_protein_stationary_covariances();

  
  return 0;
}
