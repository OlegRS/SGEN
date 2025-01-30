#include <vector>
#include "../Neuron.hpp"
#include "../compartments/Soma.hpp"
#include "../compartments/Dendritic_segment.hpp"
#include "../compartments/Spine.hpp"
#include "../engines/analytic/Analytic_engine.hpp"

using namespace std;

#define N_DS 400
#define LENGTH 5000.

int main() {

  Soma soma("soma", LENGTH/N_DS);
  
  // Constructing a neuron
  Dendritic_segment* p_ds = new Dendritic_segment(soma, "d_1-1", LENGTH/N_DS);
  
  for(unsigned int i=0; i<N_DS-1; ++i)
    p_ds = new Dendritic_segment(*p_ds, "d_1-" + to_string(i+2), LENGTH/N_DS);

  Neuron neuron(soma, "Test_neuron");
  
  cout << neuron << endl;

  Analytic_engine ae(neuron);

  // Computing stationary moments by matrix inversion
  // ae.sem_stationary_expectations().sem_stationary_covariances();


#define t_start 0
#define t1 500
  double dt = .01;


#define dim ((N_DS + 1)*2 + 1)

  arma::vec expectations(dim);
  arma::mat covariances(dim,dim);

  ofstream ofs_mRNA_expectations("sem_mRNA_expectations_400_always_active_gene"),
    ofs_gene_mRNA_covariances("sem_gene_mRNA_covariances_400_always_active_gene"),
    ofs_mRNA_covariances("sem_mRNA_covariances_400_always_active_gene"),
    ofs_prot_expectations("sem_protein_expectations_400_always_active_gene"),
    ofs_prot_covariances("sem_protein_covariances_400_always_active_gene");

  for(double t=t_start; t<t1; t+=dt) {
    std::cout << "t=" << t << std::endl;
    // expectations = ae.get_expectations();
    // ofs_expectations << t;
    // for(size_t i=0; i<dim; ++i)
    //   ofs_expectations  << ',' << expectations(i);
    // ofs_expectations << endl;
    
    if(t==0)
      ae.sem_nonstationary_expectations_direct_ODE_solver_step(dt,true).sem_nonstationary_covariances_direct_ODE_solver_step(dt);
    else
      ae.sem_nonstationary_expectations_direct_ODE_solver_step(dt).sem_nonstationary_covariances_direct_ODE_solver_step(dt);
  }

  // Writing output
  expectations = ae.get_expectations();
  ofs_mRNA_expectations << expectations(1);
  for(size_t i=2; i<1+N_DS+1; ++i)
    ofs_mRNA_expectations  << ',' << expectations(i);
  ofs_mRNA_expectations << endl;

  ofs_prot_expectations << expectations(1+N_DS+1);
  for(size_t i=1+N_DS+2; i<dim; ++i)
    ofs_prot_expectations  << ',' << expectations(i);
  ofs_prot_expectations << endl;

  covariances = ae.get_covariances();
  ofs_mRNA_covariances << covariances(1,1) - expectations(1);
  for(size_t i=2; i<1+N_DS+1; ++i)
    ofs_mRNA_covariances  << ',' << covariances(i,i) - expectations(i);
  ofs_mRNA_covariances << endl;

  ofs_gene_mRNA_covariances << covariances(0,1);
  for(size_t i=2; i<1+N_DS+1; ++i)
    ofs_gene_mRNA_covariances  << ',' << covariances(0,i);
  ofs_gene_mRNA_covariances << endl;

  
  // ofs_mRNA_correlations << "t=" << t1 << std::endl;
  // for(size_t i=0; i<dim; ++i) {
  //   for(size_t j=0; j<dim; ++j)
  //     ofs_mRNA_correlations << mRNA_covariances(i,j) - mRNA_expectations(i)*mRNA_expectations(j) << ',';
  //   ofs_mRNA_correlations << std::endl;
  // }

  
  return 0;
}
