#ifndef __ANALYTIC_ENGINE_HPP__
#define __ANALYTIC_ENGINE_HPP__

#include <armadillo>

#include "../../Neuron.hpp"
#include "../../compartments/Soma.hpp"
#include "../../compartments/Dendritic_segment.hpp"
#include "../../compartments/Spine.hpp"
#include "../../junctions/Som_den_junction.hpp"
#include "../../junctions/Den_syn_junction.hpp"
#include "../../junctions/Den_den_junction.hpp"

class Analytic_engine {
  // Parameters
  Neuron* p_neuron = nullptr;
  unsigned int o1_dim, o2_dim;
  arma::mat *p_Ap, *p_Am, *p_H,
    *p_mRNA_Ap, *p_mRNA_Am, *p_mRNA_H,
    *p_prot_Ap, *p_prot_Am, *p_prot_H,
    *p_PM, //This is the diagonal of protein-mRNA Am (which is somewhat diagonal, but not square) for separate ODEs
    *p_o1_mat, *p_o2_mat,
    o1_mRNA_matrix, o1_prot_matrix,
    *o2_gene_mRNA_mat, *o2_gene_prot_mat,
    *o2_mRNA_mRNA_mat, *o2_mRNA_prot_mat, *o2_prot_prot_mat,
    *o2_nonstationary_RHS_mat,
    *p_cov_mat, *p_mRNA_mRNA_cov_mat,
    *p_mRNA_prot_cov_mat, *p_prot_prot_cov_mat;
  arma::vec *p_b,
    *p_o1_RHS, *p_o2_RHS,
    expectations, *p_covariances,
    o1_mRNA_RHS, o1_prot_RHS,
    mRNA_expectations, protein_expectations,
    *o2_gene_mRNA, *o2_gene_mRNA_RHS,
    *o2_gene_prot, *o2_gene_prot_RHS,
    *o2_mRNA_mRNA, *o2_mRNA_mRNA_RHS,
    *o2_mRNA_prot, *o2_mRNA_prot_RHS,
    *o2_prot_prot, *o2_prot_prot_RHS;

  // Variables
  std::vector<std::string> o1_var_names, o1_mRNA_names, o1_prot_names,  *p_o2_var_names;
  std::vector<double*> p_o1_vars,
    p_o1_mRNA_expectations, p_o1_prot_expectations; // Pointers to variables within the compartments
  // /// These are needed to fill o1_matrix
  // unsigned int parent_start_ind = 0;
  // unsigned int desc_start_ind = 3;

  Analytic_engine& initialise_As_and_bs();
  Analytic_engine& initialise_mRNA_mRNA_cov_mat();
  Analytic_engine& initialise_mRNA_As();
  Analytic_engine& initialise_prot_As();
  Analytic_engine& initialise_o1_mat_and_RHS();
  // Sets 1st order matrix starting from the given compartment
  const Compartment* set_As_and_bs_soma();
  void set_As(const Compartment&);

  const Compartment* sem_set_As_and_bs_soma();
  void sem_set_As(const Compartment&);

  const Compartment* set_mRNA_As_soma();
  void set_mRNA_As(const Compartment&);
  const Compartment* set_prot_As_soma();
  void set_prot_As(const Compartment&);
  const Compartment* set_PM_soma();
  void set_PM(const Compartment&);  

  const Compartment* update_prot_source_soma();
  void update_prot_source(const Compartment&);
  
  Analytic_engine& initialise_hopping_rate_matrix();
  void set_hopping_rate_matrix(const Compartment&);
  void sem_set_hopping_rate_matrix(const Compartment&);

  Analytic_engine& initialise_mRNA_hopping_rate_matrix();
  void set_mRNA_hopping_rate_matrix(const Compartment&);
  Analytic_engine& initialise_prot_hopping_rate_matrix();
  void set_prot_hopping_rate_matrix(const Compartment&);

  const Compartment* set_o1_soma();
  void set_o1_matrix(const Compartment&);
  const Compartment* sem_set_o1_soma();
  void sem_set_o1_matrix(const Compartment&);

  // 1st order with separations of gene-mRNA-protein dynamics
  const Compartment* set_o1_mRNA_soma(); //Computes expectation of active genes and sets RHS for mRNA eqns
  const Compartment* set_o1_prot_soma();
  void set_o1_mRNA_matrix(const Compartment&);
  void set_o1_prot_matrix(const Compartment&);

  // 2nd order computation
  void initialise_o2();
  void set_prot_index_from(Compartment& compartment); // Needed for sem_o2
  void sem_initialise_o2();
  const Compartment* sem_set_soma();
  void sem_set_expectations(const Compartment& parent);
  void set_o2_soma();
  void set_o2_matrix();
  void set_o2_nonstationary_RHS_soma();
  void sem_set_o2_nonstationary_RHS_soma();
  void set_o2_nonstationary_RHS_mat();
  void sem_set_o2_nonstationary_RHS_mat();

  void sem_set_o2_soma();
  void sem_set_o2_matrix();

  double o2_gene_gene;
  const Compartment* set_o2_gene_mRNA_soma();
  void set_o2_gene_mRNA_matrix(const Compartment&);
  const Compartment* set_o2_gene_prot_soma();
  void set_o2_gene_prot_matrix(const Compartment&);
  const Compartment* set_o2_mRNA_mRNA_soma();
  void set_o2_mRNA_mRNA_matrix(const Compartment&);
  const Compartment* set_o2_mRNA_prot_soma();
  void set_o2_mRNA_prot_matrix(const Compartment&);
  const Compartment* set_o2_prot_prot_soma();
  void set_o2_prot_prot_matrix(const Compartment&);
   
  void set_o2_RHS();
  void sem_set_o2_RHS();
  Analytic_engine& clear_o1_mat_and_RHS();
  Analytic_engine& clear_As_and_bs();
  Analytic_engine& clear_mRNA_As();
  Analytic_engine& clear_prot_As();
  Analytic_engine& clear_hopping_rate_matrix();
  Analytic_engine& clear_mRNA_hopping_rate_matrix();
  Analytic_engine& clear_prot_hopping_rate_matrix();
  Analytic_engine& clear_o1();
  Analytic_engine& clear_o2_mat_and_RHS();
  Analytic_engine& clear_o2();
    
public:

  size_t o2_ind(const size_t &i, const size_t &j, const size_t &o1_dim) const; // o2 index conversion
  size_t o2_ind(const size_t &i, const size_t &j) const {return o2_ind(i,j,o1_dim);}
  size_t o2_ind_asym(const size_t &i, const size_t &j, const size_t &dim_x) const {return dim_x*i+j;}
  inline size_t sem_o2_ind(const size_t &i, const size_t &j) const; // semantic o2 index conversion

  
  Analytic_engine(Neuron& neuron, bool cov_mat_init=false) :
    p_neuron(&neuron),
    o1_dim(3 + 2*neuron.p_dend_segments.size() + neuron.p_synapses.size()),
    o2_dim(o1_dim*(o1_dim+1)/2),
    p_Ap(nullptr), p_Am(nullptr), p_H(nullptr), p_b(nullptr),
    p_mRNA_Ap(nullptr), p_mRNA_Am(nullptr), p_mRNA_H(nullptr),
    p_prot_Ap(nullptr), p_prot_Am(nullptr), p_prot_H(nullptr),
    p_PM(nullptr),
    p_mRNA_mRNA_cov_mat(nullptr), p_prot_prot_cov_mat(nullptr), p_mRNA_prot_cov_mat(nullptr),
    p_o1_mat(nullptr), p_o1_RHS(nullptr),
    o1_mRNA_matrix(1+neuron.p_dend_segments.size(), 1+neuron.p_dend_segments.size()),
    o1_prot_matrix(1+neuron.p_dend_segments.size()+neuron.p_synapses.size(), 1+neuron.p_dend_segments.size()+neuron.p_synapses.size()),
    o1_mRNA_RHS(1+neuron.p_dend_segments.size()),
    mRNA_expectations(1+neuron.p_dend_segments.size()),
    o1_mRNA_names(1+neuron.p_dend_segments.size()),
    o1_prot_RHS(1+neuron.p_dend_segments.size()+neuron.p_synapses.size()),
    protein_expectations(1+neuron.p_dend_segments.size()+neuron.p_synapses.size()),
    o1_prot_names(1+neuron.p_dend_segments.size()+neuron.p_synapses.size()),
    p_o2_mat(nullptr), p_o2_RHS(nullptr), p_o2_var_names(nullptr), p_covariances(nullptr),
    expectations(o1_dim),
    o1_var_names(o1_dim),
    p_o1_vars(o1_dim),
    p_o1_mRNA_expectations(1+neuron.p_dend_segments.size()),
    p_o1_prot_expectations(1+neuron.p_dend_segments.size()+neuron.p_synapses.size()),
    o2_gene_gene(neuron.p_soma->gene_activation_rate*(neuron.p_soma->number_of_gene_copies-1)/(neuron.p_soma->gene_activation_rate + neuron.p_soma->gene_deactivation_rate)*neuron.p_soma->n_active_genes_expectation),
    o2_gene_mRNA(nullptr), o2_gene_mRNA_RHS(nullptr), o2_gene_mRNA_mat(nullptr),
    o2_gene_prot(nullptr), o2_gene_prot_RHS(nullptr), o2_gene_prot_mat(nullptr),
    o2_mRNA_mRNA(nullptr), o2_mRNA_mRNA_RHS(nullptr), o2_mRNA_mRNA_mat(nullptr),
    o2_mRNA_prot(nullptr), o2_mRNA_prot_RHS(nullptr), o2_mRNA_prot_mat(nullptr),
    o2_prot_prot(nullptr), o2_prot_prot_RHS(nullptr), o2_prot_prot_mat(nullptr),
    o2_nonstationary_RHS_mat(nullptr)
  {
    neuron.prot_ind = neuron.p_dend_segments.size()+1;
    p_cov_mat = cov_mat_init ? new arma::mat(o1_dim, o1_dim) : nullptr;
  }
  
  Analytic_engine& initialize();

  size_t o2_dimension() {return o2_dim;};
  size_t o1_dimension() {return o1_dim;};
  
  Analytic_engine& stationary_expectations();//const Neuron& neur = *p_neuron);
  Analytic_engine& sem_stationary_expectations();//const Neuron& neur = *p_neuron);
  Analytic_engine& mRNA_stationary_expectations();
  std::vector<double> stationary_mRNA_expectations(); // For Python bindings
  Analytic_engine& mRNA_o1_eigen_decomposition();
  std::vector<double> mRNA_o1_eigenvalues();
  Analytic_engine& protein_stationary_expectations();
  std::vector<double> stationary_protein_expectations(); // For Python bindings
  Analytic_engine& protein_o1_eigen_decomposition();
  std::vector<double> protein_o1_eigenvalues();
  Analytic_engine& nonstationary_expectations(const std::list<double>& times);//const Neuron& neur = *p_neuron);
  Analytic_engine& nonstationary_expectations(const double& time, const bool& reset_matrices = false, const bool& internalise=false);
  Analytic_engine& sem_nonstationary_expectations(const double& time, const bool& reset_matrices = false, const bool& internalise=false);
  Analytic_engine& sem_nonstationary_expectations(const std::list<double>& times);
  
  Analytic_engine& stationary_covariances(bool write_covariance_matrix=false);//const Neuron& neur = *p_neuron);
  Analytic_engine& sem_stationary_covariances(); // Semantically grouped stationary covariances
  Analytic_engine& sem_stationary_pearson_correlations();

  double stationary_gene_gene_covariance() {return o2_gene_gene + p_neuron->p_soma->n_active_genes_expectation;}
  Analytic_engine& gene_mRNA_stationary_covariances();
  std::vector<double> stationary_gene_mRNA_covariances();  // For Python bindings
  Analytic_engine& gene_protein_stationary_covariances();
  std::vector<double> stationary_gene_protein_covariances();  // For Python bindings
  Analytic_engine& mRNA_mRNA_stationary_covariances();
  std::vector<std::vector<double>> stationary_mRNA_mRNA_covariances();
  Analytic_engine& mRNA_protein_stationary_covariances();
  std::vector<std::vector<double>> stationary_mRNA_protein_covariances();
  Analytic_engine& protein_protein_stationary_covariances();
  std::vector<std::vector<double>> stationary_protein_protein_covariances();
  Analytic_engine& protein_protein_stationary_covariances(std::ofstream &ofs);
  Analytic_engine& nonstationary_covariances(const std::list<double>& times, arma::vec* initial_G1, arma::vec* initial_G2);
  Analytic_engine& sem_nonstationary_covariances(const std::list<double>& times, arma::vec* initial_G1, arma::vec* initial_G2);
  Analytic_engine& sem_nonstationary_covariances_using_integral(const std::list<double>& times, arma::vec* initial_G1, arma::vec* initial_G2);
  Analytic_engine& sem_nonstationary_covariances_direct_ODE_solver(const std::list<double>& times, arma::vec* initial_G1, arma::vec* initial_G2);
  Analytic_engine& sem_nonstationary_covariances_direct_ODE_solver_no_D_matrix(const std::list<double>& times, arma::vec* initial_G1, arma::vec* initial_G2);

  Analytic_engine& nonstationary_active_genes_expectations_direct_ODE_solver_step(const double& dt);
  Analytic_engine& nonstationary_mRNA_expectations_direct_ODE_solver_step(const double& dt, const bool& reset_matrices=false, const bool& internalise=false);
  Analytic_engine& nonstationary_protein_expectations_direct_ODE_solver_step(const double& dt, const bool& reset_matrices=false, const bool& internalise=false, const bool& source_update=false);

  Analytic_engine& nonstationary_gene_gene_covariances_direct_ODE_solver_step(const double& dt, const bool& reset_matrices=false);
  Analytic_engine& nonstationary_gene_mRNA_covariances_direct_ODE_solver_step(const double& dt, const bool& reset_matrices=false);
  Analytic_engine& nonstationary_mRNA_mRNA_covariances_direct_ODE_solver_step(const double& dt, const bool& reset_matrices=false);
  Analytic_engine& nonstationary_gene_prot_covariances_direct_ODE_solver_step(const double& dt, const bool& reset_matrices=false);
  Analytic_engine& nonstationary_mRNA_prot_covariances_direct_ODE_solver_step(const double& dt, const bool& reset_matrices=false);
  Analytic_engine& nonstationary_prot_prot_covariances_direct_ODE_solver_step(const double& dt, const bool& reset_matrices=false);

  std::list<std::list<std::vector<std::vector<double>>>> nonstationary_moments_direct_ODE_solver(const std::vector<double>& write_times, const double& dt, const std::string& file_name);

  Analytic_engine& stationary_expectations_and_correlations(const std::string& active_gene_expectation_file_name="active_gene_expectation.csv", const std::string& active_gene_variance_file_name="active_gene_variance.csv", const std::string& mRNA_expectations_file_name="mRNA_expectations.csv", const std::string& protein_expectations_file_name="protein_expectations.csv", const std::string& gene_mRNA_covariances_file_name="gene_mRNA_covariances.csv", const std::string& gene_prot_covariances_file_name="gene_prot_covariances.csv", const std::string& mRNA_mRNA_covariances_file_name="mRNA_mRNA_covariances.csv", const std::string& mRNA_prot_covariances_file_name="mRNA_prot_covariances.csv", const std::string& prot_prot_covariances_file_name="prot_prot_covariances.csv");
  
  Analytic_engine& stationary_expectations_and_correlations(const double& dt_mRNA, const double& dt_prot, const double& t_fin_mRNA, const double& t_fin_prot, const std::string& active_gene_expectation_file_name="active_gene_expectation.csv", const std::string& active_gene_variance_file_name="active_gene_variance.csv", const std::string& mRNA_expectations_file_name="mRNA_expectations.csv", const std::string& protein_expectations_file_name="protein_expectations.csv", const std::string& gene_mRNA_covariances_file_name="gene_mRNA_covariances.csv", const std::string& gene_prot_covariances_file_name="gene_prot_covariances.csv", const std::string& mRNA_mRNA_covariances_file_name="mRNA_mRNA_covariances.csv", const std::string& mRNA_prot_covariances_file_name="mRNA_prot_covariances.csv", const std::string& prot_prot_covariances_file_name="prot_prot_covariances.csv");

  Analytic_engine& nonstationary_expectations_direct_ODE_solver_step(const double& dt, const bool& reset_matrices=false, const bool& internalise=false);
  Analytic_engine& nonstationary_covariances_direct_ODE_solver_step(const double& dt, const bool& reset_matrices=false);
  Analytic_engine& sem_nonstationary_expectations_direct_ODE_solver_step(const double& dt, const bool& reset_matrices=false, const bool& internalise=false);
  Analytic_engine& sem_nonstationary_covariances_direct_ODE_solver_step(const double& dt, const bool& reset_matrices=false);

  Analytic_engine& internalise_expectations(); // Writes expectations into compartments
  Analytic_engine& internalise_mRNA_expectations();
  Analytic_engine& internalise_prot_expectations();

  Analytic_engine& sem_stationary_time_correlations(const std::list<double>& times);

  std::vector<std::string>& o1_mRNA_var_names() {return o1_mRNA_names;}
  std::vector<std::string>& o1_prot_var_names() {return o1_prot_names;}
  
  double active_genes_expectation();
  double mRNA_expectation(const Compartment&);
  double protein_expectation(const Compartment&);
  double gene_mRNA_correlation(const Compartment&);
  double mRNA_mRNA_correlation(const Compartment&, const Compartment&);
  double gene_protein_correlation(const Compartment&);
  double mRNA_protein_correlation(const Compartment&, const Compartment&);
  double protein_protein_correlation(const Compartment&, const Compartment&);

  arma::vec* G1() {return &expectations;}
  std::vector<std::string>* o1_variable_names() {return &o1_var_names;}
  arma::vec* G2() {return p_covariances;}

  size_t dim_o1() const {return o1_dim;}

  const arma::vec& get_expectations() {return expectations;}
  const arma::mat& get_covariances() {return *p_cov_mat;}
  const arma::vec& get_mRNA_expectations() {return mRNA_expectations;}
  const arma::vec& get_protein_expectations() {return protein_expectations;}

  const arma::vec& get_gene_mRNA_covariances() {return *o2_gene_mRNA;}
  const arma::vec& get_gene_prot_covariances() {return *o2_gene_prot;}
  const arma::mat& get_mRNA_mRNA_covariances() {return *p_mRNA_mRNA_cov_mat;}
  const arma::mat& get_mRNA_prot_covariances() {return *p_mRNA_prot_cov_mat;}
  const arma::mat& get_prot_prot_covariances() {return *p_prot_prot_cov_mat;}

  const arma::mat& get_o1_matrix() {return *p_o1_mat;}

  Analytic_engine& set_neuron(Neuron *p_n) {p_neuron=p_n; return *this;}
  
  // std::vector<std::pair<std::string,double>> expectations();
  // std::vector<std::vector<std::pair<std::string,double>>> Analytic_engine& correlations();

  Analytic_engine& clear_all();
    
  ~Analytic_engine();
};

#endif
