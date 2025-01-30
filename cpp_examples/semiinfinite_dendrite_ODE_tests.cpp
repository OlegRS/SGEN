#include <vector>
#include "../Neuron.hpp"
#include "../compartments/Soma.hpp"
#include "../compartments/Dendritic_segment.hpp"
#include "../compartments/Spine.hpp"
#include "../engines/analytic/Analytic_engine.hpp"

using namespace std;

#define N_DS 10
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
  ae.mRNA_stationary_expectations().protein_stationary_expectations();
  ae.gene_mRNA_stationary_covariances().mRNA_mRNA_stationary_covariances().gene_protein_stationary_covariances().mRNA_protein_stationary_covariances().protein_protein_stationary_covariances();


#define t_start 0
#define t1 1000
  double dt = .1;


#define dim N_DS + 1 //((N_DS + 1)*2 + 1)
  
  arma::mat mRNA_covariances(dim,dim), mRNA_prot_covariances(dim,dim), prot_prot_covariances(dim,dim);
  arma::vec mRNA_expectations(dim), mRNA_variances(dim),
    gene_mRNA_covariances(dim),
    gene_prot_covariances(dim),
    protein_expectations(dim);

  ofstream ofs_mRNA_expectations("mRNA_expectations"),
    ofs_protein_expectations("protein_expectations"),
    ofs_gene_mRNA_covariances("gene_mRNA_covariances"),
    ofs_gene_prot_covariances("gene_prot_covariances"),
    ofs_mRNA_covariances("mRNA_covariances"),
    ofs_mRNA_variances("mRNA_variances"),
    ofs_mRNA_correlations("mRNA_correlations"),
    ofs_mRNA_prot_covariances("mRNA_prot_covariances"),
    ofs_prot_prot_covariances("prot_prot_covariances");


  for(double t=t_start; t<t1; t+=dt) {
    std::cout << "t=" << t << std::endl;
    // expectations = ae.get_expectations();
    // ofs_expectations << t;
    // for(size_t i=0; i<dim; ++i)
    //   ofs_expectations  << ',' << expectations(i);
    // ofs_expectations << endl;
    
    if(t==0) {
      ae.nonstationary_active_genes_expectations_direct_ODE_solver_step(dt).nonstationary_mRNA_expectations_direct_ODE_solver_step(dt, true).nonstationary_gene_gene_covariances_direct_ODE_solver_step(dt).nonstationary_gene_mRNA_covariances_direct_ODE_solver_step(dt, true).nonstationary_mRNA_mRNA_covariances_direct_ODE_solver_step(dt, true);
      // ae.sem_nonstationary_expectations_direct_ODE_solver_step(dt,true).sem_nonstationary_covariances_direct_ODE_solver_step(dt);
    }
    else {
      ae.nonstationary_active_genes_expectations_direct_ODE_solver_step(dt).nonstationary_mRNA_expectations_direct_ODE_solver_step(dt).nonstationary_gene_gene_covariances_direct_ODE_solver_step(dt).nonstationary_gene_mRNA_covariances_direct_ODE_solver_step(dt).nonstationary_mRNA_mRNA_covariances_direct_ODE_solver_step(dt);
      // ae.sem_nonstationary_expectations_direct_ODE_solver_step(dt,true).sem_nonstationary_covariances_direct_ODE_solver_step(dt);
    }
  }

  ae.internalise_mRNA_expectations();
  
  for(double t=t_start; t<t1; t+=dt) {
    std::cout << "t=" << t << std::endl;    
    if(t==0)
      ae.nonstationary_protein_expectations_direct_ODE_solver_step(dt,true);
    else
      ae.nonstationary_protein_expectations_direct_ODE_solver_step(dt);
  }

  for(double t=t_start; t<t1; t+=dt) {
    std::cout << "t=" << t << std::endl;
    if(t==0)
      ae.nonstationary_gene_prot_covariances_direct_ODE_solver_step(dt, true);
    else
      ae.nonstationary_gene_prot_covariances_direct_ODE_solver_step(dt);
  }

  for(double t=t_start; t<t1; t+=dt) {
    std::cout << "t=" << t << std::endl;    
    if(t==0)
      ae.nonstationary_mRNA_prot_covariances_direct_ODE_solver_step(dt, true);
    else
      ae.nonstationary_mRNA_prot_covariances_direct_ODE_solver_step(dt);
  }
  
  for(double t=t_start; t<t1; t+=dt) {
    std::cout << "t=" << t << std::endl;    
    if(t==0)
      ae.nonstationary_prot_prot_covariances_direct_ODE_solver_step(dt, true);
    else
      ae.nonstationary_prot_prot_covariances_direct_ODE_solver_step(dt);
  }


  // Writing output
  mRNA_expectations = ae.get_mRNA_expectations();
  ofs_mRNA_expectations << mRNA_expectations(0);
  for(size_t i=1; i<dim; ++i)
    ofs_mRNA_expectations  << ',' << mRNA_expectations(i);
  ofs_mRNA_expectations << endl;

  protein_expectations = ae.get_protein_expectations();
  ofs_protein_expectations << protein_expectations(0);
  for(size_t i=1; i<dim; ++i)
    ofs_protein_expectations  << ',' << protein_expectations(i);
  ofs_protein_expectations << endl;

  gene_mRNA_covariances = ae.get_gene_mRNA_covariances();
  ofs_gene_mRNA_covariances << gene_mRNA_covariances(0);
  for(size_t i=1; i<dim; ++i)
    ofs_gene_mRNA_covariances  << ',' << gene_mRNA_covariances(i);
  ofs_gene_mRNA_covariances << endl;

  gene_prot_covariances = ae.get_gene_prot_covariances();
  ofs_gene_prot_covariances << gene_prot_covariances(0);
  for(size_t i=1; i<dim; ++i)
    ofs_gene_prot_covariances  << ',' << gene_prot_covariances(i);
  ofs_gene_prot_covariances << endl;

  mRNA_covariances = ae.get_mRNA_mRNA_covariances();
  ofs_mRNA_covariances << mRNA_covariances(0,0) - mRNA_expectations(0);
  for(size_t i=1; i<dim; ++i)
    ofs_mRNA_covariances  << ',' << mRNA_covariances(i,i) - mRNA_expectations(i);
  ofs_mRNA_covariances << endl;

  mRNA_prot_covariances = ae.get_mRNA_prot_covariances();
  ofs_mRNA_prot_covariances << mRNA_prot_covariances;

  prot_prot_covariances = ae.get_prot_prot_covariances();
  ofs_prot_prot_covariances << prot_prot_covariances;

  arma::vec mRNA_std(dim);
  arma::mat mRNA_PCCs(dim,dim);
    
  for(size_t i=0; i<dim; ++i) {   
    mRNA_std(i) = sqrt(mRNA_covariances(i,i) - mRNA_expectations(i)*mRNA_expectations(i));
    mRNA_PCCs(i,i) = 1;
  }

  for(size_t i=0; i<dim; ++i)
    for(size_t j=0; j<i; ++j)
      mRNA_PCCs(i,j) = mRNA_PCCs(j,i) = (mRNA_covariances(i,j) - mRNA_expectations(i)*mRNA_expectations(j))/(mRNA_std(i)*mRNA_std(j));
  
  ofs_mRNA_correlations << mRNA_PCCs; // Writing the whole correlation matrix

    
  // ofs_mRNA_correlations << "t=" << t1 << std::endl;
  // for(size_t i=0; i<dim; ++i) {
  //   for(size_t j=0; j<dim; ++j)
  //     ofs_mRNA_correlations << mRNA_covariances(i,j) - mRNA_expectations(i)*mRNA_expectations(j) << ',';
  //   ofs_mRNA_correlations << std::endl;
  // }

  
  return 0;
}
