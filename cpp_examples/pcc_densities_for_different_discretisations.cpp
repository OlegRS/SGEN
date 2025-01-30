#include <vector>
#include "../Neuron.hpp"
#include "../compartments/Soma.hpp"
#include "../compartments/Dendritic_segment.hpp"
#include "../compartments/Spine.hpp"
#include "../engines/analytic/Analytic_engine.hpp"

using namespace std;

#define N_DS_MAX 200
#define N_POINTS 50
#define DS_LENGTH 500.

int main() {

  for(size_t n_ds=1; n_ds <= N_DS_MAX+1; n_ds+=N_DS_MAX/N_POINTS) {

    ofstream ofs("../../data/0p1_percent_translation_rate/DS_LENGTH" + to_string(int(DS_LENGTH)) + "N_DS_" + to_string(int(n_ds)));

    Soma soma("soma", DS_LENGTH/n_ds);

    cout << "Transcription_rate_constant = " << soma.get_transcription_rate() << endl;
  
    // Constructing a neuron
    Dendritic_segment* p_ds = new Dendritic_segment(soma, "d_1-1", DS_LENGTH/n_ds);
  
    for(unsigned int i=0; i<n_ds-1; ++i)
      p_ds = new Dendritic_segment(*p_ds, "d_1-" + to_string(i+2), DS_LENGTH/n_ds);

    Neuron neuron(soma, "Test_neuron");
  
    cout << neuron << endl;


    Analytic_engine ae(neuron);
  
    //  ae.sem_stationary_expectations().sem_stationary_pearson_correlations();


    ae.mRNA_stationary_expectations().protein_stationary_expectations();
    
    ae.gene_mRNA_stationary_covariances().gene_protein_stationary_covariances().mRNA_mRNA_stationary_covariances().mRNA_protein_stationary_covariances().protein_protein_stationary_covariances(ofs);

    ofs.close();
  }

  
  return 0;
}
