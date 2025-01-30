#include <vector>
#include "../include/Neuron.hpp"
#include "../include/compartments/Soma.hpp"
#include "../include/compartments/Dendritic_segment.hpp"
#include "../include/compartments/Spine.hpp"
#include "../include/engines/analytic/Analytic_engine.hpp"

using namespace std;

#define N_DS 10
#define LENGTH 5000.

int main() {

  Soma soma("soma", LENGTH/N_DS);
  
  // Constructing a neuron
  Dendritic_segment* p_ds = new Dendritic_segment(soma, "d_1-1", LENGTH/N_DS);

  // Spine* s = new Spine(*p_ds);
  
  for(unsigned int i=0; i<N_DS-1; ++i)
    p_ds = new Dendritic_segment(*p_ds, "d_1-" + to_string(i+2), LENGTH/N_DS);

  Neuron neuron(soma, "Test_neuron");
  
  cout << neuron << endl;


  Analytic_engine ae(neuron);

  ae.stationary_expectations_and_correlations();
  
  return 0;
}
