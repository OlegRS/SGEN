#include "../../../include/engines/analytic/Analytic_engine.hpp"
#include <math.h>

const Compartment* Analytic_engine::set_o1_mRNA_soma() {
  auto& soma = *p_neuron->p_soma;
  o1_mRNA_matrix(0,0) -= soma.mRNA_decay_rate;
  o1_mRNA_RHS(0) = -soma.n_active_genes_expectation * soma.transcription_rate;

  o1_mRNA_names[0] = soma.name;
  p_o1_mRNA_expectations[0] = &soma.n_mRNA_expectation;

  return p_neuron->p_soma;
}

void Analytic_engine::set_o1_mRNA_matrix(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty())
    return;

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;
    
    size_t parent_start_ind = p_junc -> p_from -> mRNA_ind-1;
    size_t desc_start_ind = p_junc -> p_to -> mRNA_ind-1;
    
    if(p_junc->type() != DEN_SYN) {

      o1_mRNA_names[desc_start_ind] = p_junc->p_to->name;
      
      p_o1_mRNA_expectations[desc_start_ind] = &p_junc->p_to->n_mRNA_expectation;

      o1_mRNA_matrix(desc_start_ind, desc_start_ind) -= p_junc->p_to->mRNA_decay_rate;

      o1_mRNA_matrix(parent_start_ind, parent_start_ind) -= p_junc->fwd_mRNA_hop_rate;
      o1_mRNA_matrix(parent_start_ind, desc_start_ind) += p_junc->bkwd_mRNA_hop_rate;
      o1_mRNA_matrix(desc_start_ind, parent_start_ind) += p_junc->fwd_mRNA_hop_rate;
      o1_mRNA_matrix(desc_start_ind, desc_start_ind) -= p_junc->bkwd_mRNA_hop_rate;
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_o1_mRNA_matrix(*((*it_p_junc)->p_to));
}

Analytic_engine& Analytic_engine::mRNA_stationary_expectations() {

  set_o1_mRNA_matrix(*set_o1_mRNA_soma());
  mRNA_expectations = o1_mRNA_matrix.i()*o1_mRNA_RHS;
  
  return internalise_mRNA_expectations();
}

std::vector<double> Analytic_engine::stationary_mRNA_expectations() {

  mRNA_stationary_expectations();

  size_t dim = p_neuron->p_dend_segments.size() + p_neuron->p_synapses.size() + 1;
  std::vector<double> expectations(dim);

  expectations[0] = mRNA_expectations(0);
  for(auto& p_junc : p_neuron -> p_junctions)
    expectations[p_junc->p_to->id] = p_junc->p_to->n_mRNA_expectation;
  
  return expectations;
}

Analytic_engine& Analytic_engine::mRNA_o1_eigen_decomposition() {
  std::cerr << "Computing eigen decomposition...\n";
  arma::cx_vec eigval_c;
  arma::vec eigval(1+p_neuron->p_dend_segments.size());
  arma::cx_mat eigvec_c;
  arma::mat trans(1+p_neuron->p_dend_segments.size(), 1+p_neuron->p_dend_segments.size()); // Transition matrix
  arma::eig_gen(eigval_c, eigvec_c, o1_mRNA_matrix);

  std::cout << "mRNA_eigval:\n";
  for(size_t i=0; i<1+p_neuron->p_dend_segments.size(); ++i) {
    std::cout << (eigval(i) = -eigval_c(i).real()) << ',';
    for(size_t j=0; j<1+p_neuron->p_dend_segments.size(); ++j)
      trans(i,j) = -eigvec_c(i,j).real();
  }
  std::cout << "\nmRNA_eigvec:\n" << trans << std::endl;

  std::cout << "Compartment names for mRNA:\n";
  for(auto& name : o1_mRNA_names)
    std::cout << name + "__mRNA" << std::endl;

  return *this;
}

Analytic_engine& Analytic_engine::internalise_mRNA_expectations() {
  size_t mRNA_size = 1 + p_neuron->p_dend_segments.size();
  for(size_t i=0; i<mRNA_size; ++i)
    *p_o1_mRNA_expectations[i] = mRNA_expectations(i);
  return *this;
}

const Compartment* Analytic_engine::set_o1_prot_soma() {
  auto& soma = *p_neuron->p_soma;
  o1_prot_matrix(0,0) -= soma.protein_decay_rate;
  o1_prot_RHS(0) = -mRNA_expectations(0)*soma.translation_rate;

  o1_prot_names[0] = soma.name;
  p_o1_prot_expectations[0] = &soma.n_prot_expectation;

  return &soma;
}

void Analytic_engine::set_o1_prot_matrix(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty())
    return;

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;
    
    size_t parent_start_ind = p_junc -> p_from -> id;
    size_t desc_start_ind = p_junc -> p_to -> id;
    
    o1_prot_names[desc_start_ind] = p_junc->p_to->name;      
    p_o1_prot_expectations[desc_start_ind] = &p_junc->p_to->n_prot_expectation;

    o1_prot_matrix(parent_start_ind, parent_start_ind) -= p_junc->fwd_prot_hop_rate;
    o1_prot_matrix(parent_start_ind, desc_start_ind) += p_junc->bkwd_prot_hop_rate;
    o1_prot_matrix(desc_start_ind, parent_start_ind) += p_junc->fwd_prot_hop_rate;
    o1_prot_matrix(desc_start_ind, desc_start_ind) -= p_junc->bkwd_prot_hop_rate;


    if(p_junc->type() != DEN_SYN) {// Protein decay everywhere apart from synapses
      o1_prot_matrix(desc_start_ind, desc_start_ind) -= p_junc->p_to->protein_decay_rate;
      o1_prot_RHS(desc_start_ind) = -(p_junc->p_to->n_mRNA_expectation)*(p_junc->p_to->translation_rate);
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_o1_prot_matrix(*((*it_p_junc)->p_to));
}

Analytic_engine& Analytic_engine::protein_stationary_expectations() {

  set_o1_prot_matrix(*set_o1_prot_soma());
  protein_expectations = o1_prot_matrix.i()*o1_prot_RHS;
  
  return internalise_prot_expectations();
}

std::vector<double> Analytic_engine::stationary_protein_expectations() {

  protein_stationary_expectations();

  size_t dim = p_neuron->p_dend_segments.size() + p_neuron->p_synapses.size() + 1;
  
  std::vector<double> expectations(dim);
  for(size_t i=0; i<dim; ++i)
    expectations[i] = protein_expectations(i);
  
  return expectations;
}

Analytic_engine& Analytic_engine::protein_o1_eigen_decomposition() {
  std::cerr << "Computing eigen decomposition...\n";
  size_t sz = 1+p_neuron->p_dend_segments.size() + p_neuron->p_synapses.size();
  arma::cx_vec eigval_c;
  arma::vec eigval(sz);
  arma::cx_mat eigvec_c;
  arma::mat trans(sz, sz); // Transition matrix
  arma::eig_gen(eigval_c, eigvec_c, o1_prot_matrix);

  std::cout << "Protein_eigval:\n";
  for(size_t i=0; i<sz; ++i) {
    std::cout << (eigval(i) = -eigval_c(i).real()) << ',';
    for(size_t j=0; j<sz; ++j)
      trans(i,j) = -eigvec_c(i,j).real();
  }
  std::cout << "\nProtein_eigvec:\n" << trans << std::endl;

  std::cout << "Compartment names for proteins:\n";
  for(auto& name : o1_prot_names)
    std::cout << name + "__Prot" << std::endl;

  return *this;
}

Analytic_engine& Analytic_engine::internalise_prot_expectations() {
  size_t prot_size = 1 + p_neuron->p_dend_segments.size() + p_neuron->p_synapses.size();
  for(size_t i=0; i<prot_size; ++i)
    *p_o1_prot_expectations[i] = protein_expectations(i);
  return *this;
}

Analytic_engine& Analytic_engine::initialise_o1_mat_and_RHS() {
  clear_o1_mat_and_RHS();
  p_o1_mat = new arma::mat(o1_dim, o1_dim);
  p_o1_RHS = new arma::vec(o1_dim);
  
  return *this;
}

Analytic_engine& Analytic_engine::initialise_As_and_bs() {
  clear_As_and_bs();

  p_Ap = new arma::mat(o1_dim, o1_dim);
  p_Am = new arma::mat(o1_dim, o1_dim);
  
  p_b = new arma::vec(o1_dim);
  
  return *this;
}

Analytic_engine& Analytic_engine::initialise_mRNA_mRNA_cov_mat() {
  p_mRNA_mRNA_cov_mat = new arma::mat(o1_dim, o1_dim);
  return *this;
}

Analytic_engine& Analytic_engine::initialise_mRNA_As() {
  clear_mRNA_As();

  size_t mRNA_dim = 1+p_neuron->p_dend_segments.size();
  p_mRNA_Ap = new arma::mat(mRNA_dim, mRNA_dim);
  p_mRNA_Am = new arma::mat(mRNA_dim, mRNA_dim);
  
  return *this;
}

Analytic_engine& Analytic_engine::initialise_prot_As() {
  clear_prot_As();

  size_t prot_dim = 1+p_neuron->p_dend_segments.size()+p_neuron->p_synapses.size();

  p_prot_Ap = new arma::mat(prot_dim, prot_dim);
  p_prot_Am = new arma::mat(prot_dim, prot_dim);
  
  return *this;
}

Analytic_engine& Analytic_engine::initialise_hopping_rate_matrix() {
  clear_hopping_rate_matrix();

  p_H = new arma::mat(o1_dim, o1_dim);
  
  return *this;
}

Analytic_engine& Analytic_engine::initialise_mRNA_hopping_rate_matrix() {
  clear_mRNA_hopping_rate_matrix();

  size_t mRNA_dim = 1+p_neuron->p_dend_segments.size();
  p_mRNA_H = new arma::mat(mRNA_dim, mRNA_dim);
  
  return *this;
}

Analytic_engine& Analytic_engine::initialise_prot_hopping_rate_matrix() {
  clear_prot_hopping_rate_matrix();

  size_t prot_dim = 1+p_neuron->p_dend_segments.size()+p_neuron->p_synapses.size();
  p_prot_H = new arma::mat(prot_dim, prot_dim);
  
  return *this;
}

const Compartment* Analytic_engine::set_As_and_bs_soma() {
  initialise_As_and_bs();
  
  auto& soma = *p_neuron->p_soma;

  (*p_Am)(0,0) = -soma.gene_activation_rate - soma.gene_deactivation_rate;
  (*p_Am)(1,0) = soma.transcription_rate;
  (*p_Am)(1,1) = -soma.mRNA_decay_rate;
  (*p_Am)(2,1) = soma.translation_rate;
  (*p_Am)(2,2) = -soma.protein_decay_rate;

  (*p_Ap)(0,0) = -soma.gene_activation_rate + soma.gene_deactivation_rate;
  (*p_Ap)(1,0) = soma.transcription_rate;
  (*p_Ap)(1,1) = soma.mRNA_decay_rate;
  (*p_Ap)(2,1) = soma.translation_rate;
  (*p_Ap)(2,2) = soma.protein_decay_rate;

  o1_var_names[0] = soma.name + "__Gene";
  o1_var_names[1] = soma.name + "__mRNA";
  o1_var_names[2] = soma.name + "__Prot";

  p_o1_vars[0] = &soma.n_active_genes_expectation;
  p_o1_vars[1] = &soma.n_mRNA_expectation;
  p_o1_vars[2] = &soma.n_prot_expectation;

  (*p_b)(0) = (p_neuron->p_soma->gene_activation_rate) * (p_neuron->p_soma->number_of_gene_copies);
  
  return p_neuron->p_soma;
}

const Compartment* Analytic_engine::sem_set_As_and_bs_soma() {
  initialise_As_and_bs();
  auto& soma = *p_neuron->p_soma;

  set_prot_index_from(soma);

  (*p_Am)(0,0) = -soma.gene_activation_rate - soma.gene_deactivation_rate;
  (*p_Am)(1,0) = soma.transcription_rate;
  (*p_Am)(1,1) = -soma.mRNA_decay_rate;
  (*p_Am)(soma.prot_ind,1) = soma.translation_rate;
  (*p_Am)(soma.prot_ind,soma.prot_ind) = -soma.protein_decay_rate;

  (*p_Ap)(0,0) = -soma.gene_activation_rate + soma.gene_deactivation_rate;
  (*p_Ap)(1,0) = soma.transcription_rate;
  (*p_Ap)(1,1) = soma.mRNA_decay_rate;
  (*p_Ap)(soma.prot_ind,1) = soma.translation_rate;
  (*p_Ap)(soma.prot_ind,soma.prot_ind) = soma.protein_decay_rate;

  o1_var_names[0] = soma.name + "__Gene";
  o1_var_names[1] = soma.name + "__mRNA";
  o1_var_names[soma.prot_ind] = soma.name + "__Prot";

  p_o1_vars[0] = &soma.n_active_genes_expectation;
  p_o1_vars[1] = &soma.n_mRNA_expectation;
  p_o1_vars[soma.prot_ind] = &soma.n_prot_expectation;

  (*p_b)(0) = (p_neuron->p_soma->gene_activation_rate) * (p_neuron->p_soma->number_of_gene_copies);
  
  return p_neuron->p_soma;
}


void Analytic_engine::set_As(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty())
    return;

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;
    
    size_t& parent_start_ind = p_junc -> p_from -> o1_index;
    size_t& desc_start_ind = p_junc -> p_to -> o1_index;
    
    if(p_junc->type() == DEN_SYN) {
      o1_var_names[p_junc->p_to->o1_index] = p_junc->p_to->name + "__Prot";
      p_o1_vars[p_junc->p_to->o1_index] = &p_junc->p_to->n_prot_expectation;
      
      (*p_Am)(parent_start_ind+1, parent_start_ind+1) -= p_junc->fwd_prot_hop_rate;
      (*p_Am)(parent_start_ind+1, desc_start_ind) += p_junc->bkwd_prot_hop_rate;
      (*p_Am)(desc_start_ind, parent_start_ind+1) += p_junc->fwd_prot_hop_rate;
      (*p_Am)(desc_start_ind, desc_start_ind) -= p_junc->bkwd_prot_hop_rate;      
      (*p_Am)(desc_start_ind, desc_start_ind) -= p_junc->p_to->protein_decay_rate;

      (*p_Ap)(parent_start_ind+1, parent_start_ind+1) += p_junc->fwd_prot_hop_rate;
      (*p_Ap)(parent_start_ind+1, desc_start_ind) += p_junc->bkwd_prot_hop_rate;
      (*p_Ap)(desc_start_ind, parent_start_ind+1) += p_junc->fwd_prot_hop_rate;
      (*p_Ap)(desc_start_ind, desc_start_ind) += p_junc->bkwd_prot_hop_rate;
      (*p_Ap)(desc_start_ind, desc_start_ind) += p_junc->p_to->protein_decay_rate;
    }
    else if(p_junc->type() == DEN_DEN) {
      o1_var_names[p_junc->p_to->o1_index] = p_junc->p_to->name + "__mRNA";
      o1_var_names[p_junc->p_to->o1_index+1] = p_junc->p_to->name + "__Prot";
      p_o1_vars[p_junc->p_to->o1_index] = &p_junc->p_to->n_mRNA_expectation;
      p_o1_vars[p_junc->p_to->o1_index+1] = &p_junc->p_to->n_prot_expectation;

      (*p_Am)(desc_start_ind, desc_start_ind) -= p_junc->p_to->mRNA_decay_rate;
      (*p_Am)(desc_start_ind+1, desc_start_ind) += p_junc->p_to->translation_rate;
      (*p_Am)(desc_start_ind+1, desc_start_ind+1) -= p_junc->p_to->protein_decay_rate;
      (*p_Am)(parent_start_ind, parent_start_ind) -= p_junc->fwd_mRNA_hop_rate;
      (*p_Am)(parent_start_ind, desc_start_ind) += p_junc->bkwd_mRNA_hop_rate;
      (*p_Am)(desc_start_ind, parent_start_ind) += p_junc->fwd_mRNA_hop_rate;
      (*p_Am)(desc_start_ind, desc_start_ind) -= p_junc->bkwd_mRNA_hop_rate;
      (*p_Am)(parent_start_ind+1, parent_start_ind+1) -= p_junc->fwd_prot_hop_rate;
      (*p_Am)(parent_start_ind+1, desc_start_ind+1) += p_junc->bkwd_prot_hop_rate;
      (*p_Am)(desc_start_ind+1, parent_start_ind+1) += p_junc->fwd_prot_hop_rate;
      (*p_Am)(desc_start_ind+1, desc_start_ind+1) -= p_junc->bkwd_prot_hop_rate;

      (*p_Ap)(desc_start_ind, desc_start_ind) += p_junc->p_to->mRNA_decay_rate;
      (*p_Ap)(desc_start_ind+1, desc_start_ind) += p_junc->p_to->translation_rate;
      (*p_Ap)(desc_start_ind+1, desc_start_ind+1) += p_junc->p_to->protein_decay_rate;
      (*p_Ap)(parent_start_ind, parent_start_ind) += p_junc->fwd_mRNA_hop_rate;
      (*p_Ap)(parent_start_ind, desc_start_ind) += p_junc->bkwd_mRNA_hop_rate;
      (*p_Ap)(desc_start_ind, parent_start_ind) += p_junc->fwd_mRNA_hop_rate;
      (*p_Ap)(desc_start_ind, desc_start_ind) += p_junc->bkwd_mRNA_hop_rate;
      (*p_Ap)(parent_start_ind+1, parent_start_ind+1) += p_junc->fwd_prot_hop_rate;
      (*p_Ap)(parent_start_ind+1, desc_start_ind+1) += p_junc->bkwd_prot_hop_rate;
      (*p_Ap)(desc_start_ind+1, parent_start_ind+1) += p_junc->fwd_prot_hop_rate;
      (*p_Ap)(desc_start_ind+1, desc_start_ind+1) += p_junc->bkwd_prot_hop_rate;
    }
    else if(p_junc->type() == SOM_DEN) {
      o1_var_names[p_junc->p_to->o1_index] = p_junc->p_to->name + "__mRNA";
      o1_var_names[p_junc->p_to->o1_index+1] = p_junc->p_to->name + "__Prot";
      p_o1_vars[p_junc->p_to->o1_index] = &p_junc->p_to->n_mRNA_expectation;
      p_o1_vars[p_junc->p_to->o1_index+1] = &p_junc->p_to->n_prot_expectation;
      
      (*p_Am)(desc_start_ind, desc_start_ind) -= p_junc->p_to->mRNA_decay_rate;
      (*p_Am)(desc_start_ind+1, desc_start_ind) += p_junc->p_to->translation_rate;
      (*p_Am)(desc_start_ind+1, desc_start_ind+1) -= p_junc->p_to->protein_decay_rate;
      (*p_Am)(parent_start_ind+1, parent_start_ind+1) -= p_junc->fwd_mRNA_hop_rate;
      (*p_Am)(parent_start_ind+1, desc_start_ind) += p_junc->bkwd_mRNA_hop_rate;
      (*p_Am)(desc_start_ind, parent_start_ind+1) += p_junc->fwd_mRNA_hop_rate;
      (*p_Am)(desc_start_ind, desc_start_ind) -= p_junc->bkwd_mRNA_hop_rate;
      (*p_Am)(parent_start_ind+2, parent_start_ind+2) -= p_junc->fwd_prot_hop_rate;
      (*p_Am)(parent_start_ind+2, desc_start_ind+1) += p_junc->bkwd_prot_hop_rate;
      (*p_Am)(desc_start_ind+1, parent_start_ind+2) += p_junc->fwd_prot_hop_rate;
      (*p_Am)(desc_start_ind+1, desc_start_ind+1) -= p_junc->bkwd_prot_hop_rate;

      (*p_Ap)(desc_start_ind, desc_start_ind) += p_junc->p_to->mRNA_decay_rate;
      (*p_Ap)(desc_start_ind+1, desc_start_ind) += p_junc->p_to->translation_rate;
      (*p_Ap)(desc_start_ind+1, desc_start_ind+1) += p_junc->p_to->protein_decay_rate;
      (*p_Ap)(parent_start_ind+1, parent_start_ind+1) += p_junc->fwd_mRNA_hop_rate;
      (*p_Ap)(parent_start_ind+1, desc_start_ind) += p_junc->bkwd_mRNA_hop_rate;
      (*p_Ap)(desc_start_ind, parent_start_ind+1) += p_junc->fwd_mRNA_hop_rate;
      (*p_Ap)(desc_start_ind, desc_start_ind) += p_junc->bkwd_mRNA_hop_rate;
      (*p_Ap)(parent_start_ind+2, parent_start_ind+2) += p_junc->fwd_prot_hop_rate;
      (*p_Ap)(parent_start_ind+2, desc_start_ind+1) += p_junc->bkwd_prot_hop_rate;
      (*p_Ap)(desc_start_ind+1, parent_start_ind+2) += p_junc->fwd_prot_hop_rate;
      (*p_Ap)(desc_start_ind+1, desc_start_ind+1) += p_junc->bkwd_prot_hop_rate;

    }
    else {
      std::cerr << "-----------------------------------------\n"
                << "ERROR: UNKNOWN TYPE JUNCTION FOUND\n"
                << "-----------------------------------------\n";
      exit(1);
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_As(*((*it_p_junc)->p_to));
}

void Analytic_engine::sem_set_As(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty())
    return;

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;

    size_t& parent_prot_ind = p_junc -> p_from -> prot_ind;
    size_t& desc_prot_ind = p_junc -> p_to -> prot_ind;
    size_t& parent_mRNA_ind = p_junc -> p_from -> mRNA_ind;
    size_t& desc_mRNA_ind = p_junc -> p_to -> mRNA_ind;
        
    if(p_junc->type() == DEN_SYN) {
      o1_var_names[desc_prot_ind] = p_junc->p_to->name + "__Prot";
      p_o1_vars[desc_prot_ind] = &p_junc->p_to->n_prot_expectation;

      (*p_Am)(parent_prot_ind, parent_prot_ind) -= p_junc->fwd_prot_hop_rate;
      (*p_Am)(parent_prot_ind, desc_prot_ind) += p_junc->bkwd_prot_hop_rate;
      (*p_Am)(desc_prot_ind, parent_prot_ind) += p_junc->fwd_prot_hop_rate;
      (*p_Am)(desc_prot_ind, desc_prot_ind) -= p_junc->bkwd_prot_hop_rate;      
      (*p_Am)(desc_prot_ind, desc_prot_ind) -= p_junc->p_to->protein_decay_rate;

      (*p_Ap)(parent_prot_ind, parent_prot_ind) += p_junc->fwd_prot_hop_rate;
      (*p_Ap)(parent_prot_ind, desc_prot_ind) += p_junc->bkwd_prot_hop_rate;
      (*p_Ap)(desc_prot_ind, parent_prot_ind) += p_junc->fwd_prot_hop_rate;
      (*p_Ap)(desc_prot_ind, desc_prot_ind) += p_junc->bkwd_prot_hop_rate;
      (*p_Ap)(desc_prot_ind, desc_prot_ind) += p_junc->p_to->protein_decay_rate;
    }
    else if(p_junc->type() == DEN_DEN) {
      o1_var_names[desc_mRNA_ind] = p_junc->p_to->name + "__mRNA";
      o1_var_names[desc_prot_ind] = p_junc->p_to->name + "__Prot";
      p_o1_vars[desc_mRNA_ind] = &p_junc->p_to->n_mRNA_expectation;
      p_o1_vars[desc_prot_ind] = &p_junc->p_to->n_prot_expectation;

      (*p_Am)(desc_mRNA_ind, desc_mRNA_ind) -= p_junc->p_to->mRNA_decay_rate;
      (*p_Am)(desc_prot_ind, desc_mRNA_ind) += p_junc->p_to->translation_rate;
      (*p_Am)(desc_prot_ind, desc_prot_ind) -= p_junc->p_to->protein_decay_rate;
      (*p_Am)(parent_mRNA_ind, parent_mRNA_ind) -= p_junc->fwd_mRNA_hop_rate;
      (*p_Am)(parent_mRNA_ind, desc_mRNA_ind) += p_junc->bkwd_mRNA_hop_rate;
      (*p_Am)(desc_mRNA_ind, parent_mRNA_ind) += p_junc->fwd_mRNA_hop_rate;
      (*p_Am)(desc_mRNA_ind, desc_mRNA_ind) -= p_junc->bkwd_mRNA_hop_rate;
      (*p_Am)(parent_prot_ind, parent_prot_ind) -= p_junc->fwd_prot_hop_rate;
      (*p_Am)(parent_prot_ind, desc_prot_ind) += p_junc->bkwd_prot_hop_rate;
      (*p_Am)(desc_prot_ind, parent_prot_ind) += p_junc->fwd_prot_hop_rate;
      (*p_Am)(desc_prot_ind, desc_prot_ind) -= p_junc->bkwd_prot_hop_rate;

      (*p_Ap)(desc_mRNA_ind, desc_mRNA_ind) += p_junc->p_to->mRNA_decay_rate;
      (*p_Ap)(desc_prot_ind, desc_mRNA_ind) += p_junc->p_to->translation_rate;
      (*p_Ap)(desc_prot_ind, desc_prot_ind) += p_junc->p_to->protein_decay_rate;
      (*p_Ap)(parent_mRNA_ind, parent_mRNA_ind) += p_junc->fwd_mRNA_hop_rate;
      (*p_Ap)(parent_mRNA_ind, desc_mRNA_ind) += p_junc->bkwd_mRNA_hop_rate;
      (*p_Ap)(desc_mRNA_ind, parent_mRNA_ind) += p_junc->fwd_mRNA_hop_rate;
      (*p_Ap)(desc_mRNA_ind, desc_mRNA_ind) += p_junc->bkwd_mRNA_hop_rate;
      (*p_Ap)(parent_prot_ind, parent_prot_ind) += p_junc->fwd_prot_hop_rate;
      (*p_Ap)(parent_prot_ind, desc_prot_ind) += p_junc->bkwd_prot_hop_rate;
      (*p_Ap)(desc_prot_ind, parent_prot_ind) += p_junc->fwd_prot_hop_rate;
      (*p_Ap)(desc_prot_ind, desc_prot_ind) += p_junc->bkwd_prot_hop_rate;

    }
    else if(p_junc->type() == SOM_DEN) {
      o1_var_names[desc_mRNA_ind] = p_junc->p_to->name + "__mRNA";
      o1_var_names[desc_prot_ind] = p_junc->p_to->name + "__Prot";
      p_o1_vars[desc_mRNA_ind] = &p_junc->p_to->n_mRNA_expectation;
      p_o1_vars[desc_prot_ind] = &p_junc->p_to->n_prot_expectation;
      
      (*p_Am)(desc_mRNA_ind, desc_mRNA_ind) -= p_junc->p_to->mRNA_decay_rate;
      (*p_Am)(desc_prot_ind, desc_mRNA_ind) += p_junc->p_to->translation_rate;
      (*p_Am)(desc_prot_ind, desc_prot_ind) -= p_junc->p_to->protein_decay_rate;
      (*p_Am)(parent_mRNA_ind, parent_mRNA_ind) -= p_junc->fwd_mRNA_hop_rate;
      (*p_Am)(parent_mRNA_ind, desc_mRNA_ind) += p_junc->bkwd_mRNA_hop_rate;
      (*p_Am)(desc_mRNA_ind, parent_mRNA_ind) += p_junc->fwd_mRNA_hop_rate;
      (*p_Am)(desc_mRNA_ind, desc_mRNA_ind) -= p_junc->bkwd_mRNA_hop_rate;
      (*p_Am)(parent_prot_ind, parent_prot_ind) -= p_junc->fwd_prot_hop_rate;
      (*p_Am)(parent_prot_ind, desc_prot_ind) += p_junc->bkwd_prot_hop_rate;
      (*p_Am)(desc_prot_ind, parent_prot_ind) += p_junc->fwd_prot_hop_rate;
      (*p_Am)(desc_prot_ind, desc_prot_ind) -= p_junc->bkwd_prot_hop_rate;


      (*p_Ap)(desc_mRNA_ind, desc_mRNA_ind) += p_junc->p_to->mRNA_decay_rate;
      (*p_Ap)(desc_prot_ind, desc_mRNA_ind) += p_junc->p_to->translation_rate;
      (*p_Ap)(desc_prot_ind, desc_prot_ind) += p_junc->p_to->protein_decay_rate;
      (*p_Ap)(parent_mRNA_ind, parent_mRNA_ind) += p_junc->fwd_mRNA_hop_rate;
      (*p_Ap)(parent_mRNA_ind, desc_mRNA_ind) += p_junc->bkwd_mRNA_hop_rate;
      (*p_Ap)(desc_mRNA_ind, parent_mRNA_ind) += p_junc->fwd_mRNA_hop_rate;
      (*p_Ap)(desc_mRNA_ind, desc_mRNA_ind) += p_junc->bkwd_mRNA_hop_rate;
      (*p_Ap)(parent_prot_ind, parent_prot_ind) += p_junc->fwd_prot_hop_rate;
      (*p_Ap)(parent_prot_ind, desc_prot_ind) += p_junc->bkwd_prot_hop_rate;
      (*p_Ap)(desc_prot_ind, parent_prot_ind) += p_junc->fwd_prot_hop_rate;
      (*p_Ap)(desc_prot_ind, desc_prot_ind) += p_junc->bkwd_prot_hop_rate;
    }
    else {
      std::cerr << "-----------------------------------------\n"
                << "ERROR: UNKNOWN TYPE JUNCTION FOUND\n"
                << "-----------------------------------------\n";
      exit(1);
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    sem_set_As(*((*it_p_junc)->p_to));
}

const Compartment* Analytic_engine::set_mRNA_As_soma() {
  initialise_mRNA_As();

  auto& soma = *p_neuron->p_soma;
  (*p_mRNA_Am)(0,0) = -soma.mRNA_decay_rate;
  (*p_mRNA_Ap)(0,0) = soma.mRNA_decay_rate;
  
  o1_mRNA_names[0] = soma.name;
  p_o1_mRNA_expectations[0] = &soma.n_mRNA_expectation;
   
  return &soma;
}

void Analytic_engine::set_mRNA_As(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty())
    return;

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;
    
    size_t parent_start_ind = p_junc -> p_from -> mRNA_ind-1;
    size_t desc_start_ind = p_junc -> p_to -> mRNA_ind-1;
    
    if(p_junc->type() != DEN_SYN) {

      o1_mRNA_names[desc_start_ind] = p_junc->p_to->name;
      
      p_o1_mRNA_expectations[desc_start_ind] = &p_junc->p_to->n_mRNA_expectation;

      (*p_mRNA_Am)(desc_start_ind, desc_start_ind) -= p_junc->p_to->mRNA_decay_rate;
      (*p_mRNA_Am)(parent_start_ind, parent_start_ind) -= p_junc->fwd_mRNA_hop_rate;
      (*p_mRNA_Am)(parent_start_ind, desc_start_ind) += p_junc->bkwd_mRNA_hop_rate;
      (*p_mRNA_Am)(desc_start_ind, parent_start_ind) += p_junc->fwd_mRNA_hop_rate;
      (*p_mRNA_Am)(desc_start_ind, desc_start_ind) -= p_junc->bkwd_mRNA_hop_rate;

      (*p_mRNA_Ap)(desc_start_ind, desc_start_ind) += p_junc->p_to->mRNA_decay_rate;
      (*p_mRNA_Ap)(parent_start_ind, parent_start_ind) += p_junc->fwd_mRNA_hop_rate;
      (*p_mRNA_Ap)(parent_start_ind, desc_start_ind) += p_junc->bkwd_mRNA_hop_rate;
      (*p_mRNA_Ap)(desc_start_ind, parent_start_ind) += p_junc->fwd_mRNA_hop_rate;
      (*p_mRNA_Ap)(desc_start_ind, desc_start_ind) += p_junc->bkwd_mRNA_hop_rate;
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_mRNA_As(*((*it_p_junc)->p_to));
}

void Analytic_engine::set_mRNA_hopping_rate_matrix(const Compartment& parent) {

  if(parent.it_p_out_junctions.empty())
    return;

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;

    size_t parent_start_ind = p_junc -> p_from -> mRNA_ind-1;
    size_t desc_start_ind = p_junc -> p_to -> mRNA_ind-1;
    
    if(p_junc->type() != DEN_SYN) {
      (*p_mRNA_H)(parent_start_ind, desc_start_ind) = p_junc->fwd_mRNA_hop_rate;
      (*p_mRNA_H)(desc_start_ind, parent_start_ind) = p_junc->bkwd_mRNA_hop_rate;
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_mRNA_hopping_rate_matrix(*((*it_p_junc)->p_to));
}

void Analytic_engine::set_prot_hopping_rate_matrix(const Compartment& parent) {

  if(parent.it_p_out_junctions.empty())
    return;

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;

    size_t& parent_ind = p_junc -> p_from -> id;
    size_t& desc_ind = p_junc -> p_to -> id;
    
    (*p_prot_H)(parent_ind, desc_ind) = p_junc->fwd_prot_hop_rate;
    (*p_prot_H)(desc_ind, parent_ind) = p_junc->bkwd_prot_hop_rate;
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_prot_hopping_rate_matrix(*((*it_p_junc)->p_to));
}


const Compartment* Analytic_engine::set_prot_As_soma() {
  initialise_prot_As();

  auto& soma = *p_neuron->p_soma;
  (*p_prot_Am)(0,0) = -soma.protein_decay_rate;
  (*p_prot_Ap)(0,0) = soma.protein_decay_rate;

  o1_prot_RHS(0) = -mRNA_expectations(0)*soma.translation_rate;
  
  o1_prot_names[0] = soma.name;
  p_o1_prot_expectations[0] = &soma.n_prot_expectation;
   
  return &soma;
}

void Analytic_engine::set_prot_As(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty())
    return;

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;
    
    size_t parent_start_ind = p_junc -> p_from -> id;
    size_t desc_start_ind = p_junc -> p_to -> id;
    
    o1_prot_names[desc_start_ind] = p_junc->p_to->name;      
    p_o1_prot_expectations[desc_start_ind] = &p_junc->p_to->n_prot_expectation;

    (*p_prot_Am)(desc_start_ind, desc_start_ind) -= p_junc->p_to->protein_decay_rate;
    (*p_prot_Am)(parent_start_ind, parent_start_ind) -= p_junc->fwd_prot_hop_rate;
    (*p_prot_Am)(parent_start_ind, desc_start_ind) += p_junc->bkwd_prot_hop_rate;
    (*p_prot_Am)(desc_start_ind, parent_start_ind) += p_junc->fwd_prot_hop_rate;
    (*p_prot_Am)(desc_start_ind, desc_start_ind) -= p_junc->bkwd_prot_hop_rate;

    (*p_prot_Ap)(desc_start_ind, desc_start_ind) += p_junc->p_to->protein_decay_rate;
    (*p_prot_Ap)(parent_start_ind, parent_start_ind) += p_junc->fwd_prot_hop_rate;
    (*p_prot_Ap)(parent_start_ind, desc_start_ind) += p_junc->bkwd_prot_hop_rate;
    (*p_prot_Ap)(desc_start_ind, parent_start_ind) += p_junc->fwd_prot_hop_rate;
    (*p_prot_Ap)(desc_start_ind, desc_start_ind) += p_junc->bkwd_prot_hop_rate;

    o1_prot_RHS(desc_start_ind) = -(p_junc->p_to->n_mRNA_expectation)*(p_junc->p_to->translation_rate);
  }
  
  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_prot_As(*((*it_p_junc)->p_to));
}

const Compartment* Analytic_engine::set_PM_soma() {
  size_t mRNA_dim = 1 + p_neuron->p_dend_segments.size(),
    prot_dim = mRNA_dim + p_neuron->p_synapses.size();

  if(p_PM)
    delete p_PM;

  p_PM = new arma::mat(prot_dim, mRNA_dim);

  auto& soma = *p_neuron->p_soma;
  (*p_PM)(0) = soma.translation_rate;

  return &soma;
}

void Analytic_engine::set_PM(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty())
    return;

  for(auto& it_p_junc : parent.it_p_out_junctions) 
    if((*it_p_junc)->type() != DEN_SYN) {
      size_t prot_ind = (*it_p_junc)->p_to->id,
        mRNA_ind = (*it_p_junc)->p_to->mRNA_ind-1;
      (*p_PM)(prot_ind, mRNA_ind) = (*it_p_junc)->p_to->translation_rate;
    }
  
  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_PM(*((*it_p_junc)->p_to));
}


const Compartment* Analytic_engine::update_prot_source_soma() {
  auto& soma = *p_neuron->p_soma;
  o1_prot_RHS(0) = -mRNA_expectations(0)*soma.translation_rate;
  return &soma;
}

void Analytic_engine::update_prot_source(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty())
    return;

  for(auto& it_p_junc : parent.it_p_out_junctions)
    o1_prot_RHS((*it_p_junc)->p_to->id) = -((*it_p_junc)->p_to->n_mRNA_expectation)*((*it_p_junc)->p_to->translation_rate);
  
  for(auto& it_p_junc : parent.it_p_out_junctions)
    update_prot_source(*((*it_p_junc)->p_to));
}


void Analytic_engine::set_hopping_rate_matrix(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty())
    return;

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;
    
    size_t& parent_start_ind = p_junc -> p_from -> o1_index;
    size_t& desc_start_ind = p_junc -> p_to -> o1_index;
    
    if(p_junc->type() == DEN_SYN) {
      (*p_H)(parent_start_ind+1, desc_start_ind) = p_junc->fwd_prot_hop_rate;
      (*p_H)(desc_start_ind, parent_start_ind+1) = p_junc->bkwd_prot_hop_rate;
    }
    else if(p_junc->type() == DEN_DEN) {
      (*p_H)(parent_start_ind, desc_start_ind) = p_junc->fwd_mRNA_hop_rate;
      (*p_H)(desc_start_ind, parent_start_ind) = p_junc->bkwd_mRNA_hop_rate;
      
      (*p_H)(parent_start_ind+1, desc_start_ind+1) = p_junc->fwd_prot_hop_rate;
      (*p_H)(desc_start_ind+1, parent_start_ind+1) = p_junc->bkwd_prot_hop_rate;
    }
    else if(p_junc->type() == SOM_DEN) {
      (*p_H)(parent_start_ind+1, desc_start_ind) = p_junc->fwd_mRNA_hop_rate;
      (*p_H)(desc_start_ind, parent_start_ind+1) = p_junc->bkwd_mRNA_hop_rate;

      (*p_H)(parent_start_ind+2, desc_start_ind+1) = p_junc->fwd_prot_hop_rate;
      (*p_H)(desc_start_ind+1, parent_start_ind+2) = p_junc->bkwd_prot_hop_rate;
    }
    else {
      std::cerr << "-----------------------------------------\n"
                << "ERROR: UNKNOWN TYPE JUNCTION FOUND\n"
                << "-----------------------------------------\n";
      exit(1);
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_hopping_rate_matrix(*((*it_p_junc)->p_to));
}

void Analytic_engine::sem_set_hopping_rate_matrix(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty())
    return;

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;

    size_t& parent_prot_ind = p_junc -> p_from -> prot_ind;
    size_t& desc_prot_ind = p_junc -> p_to -> prot_ind;
    size_t& parent_mRNA_ind = p_junc -> p_from -> mRNA_ind;
    size_t& desc_mRNA_ind = p_junc -> p_to -> mRNA_ind;
    
    if(p_junc->type() == DEN_SYN) {
      (*p_H)(parent_prot_ind, desc_prot_ind) = p_junc->fwd_prot_hop_rate;
      (*p_H)(desc_prot_ind, parent_prot_ind) = p_junc->bkwd_prot_hop_rate;
    }
    else if(p_junc->type() == DEN_DEN) {
      (*p_H)(parent_mRNA_ind, desc_mRNA_ind) = p_junc->fwd_mRNA_hop_rate;
      (*p_H)(desc_mRNA_ind, parent_mRNA_ind) = p_junc->bkwd_mRNA_hop_rate;
      
      (*p_H)(parent_prot_ind, desc_prot_ind) = p_junc->fwd_prot_hop_rate;
      (*p_H)(desc_prot_ind, parent_prot_ind) = p_junc->bkwd_prot_hop_rate;
    }
    else if(p_junc->type() == SOM_DEN) {
      (*p_H)(parent_mRNA_ind, desc_mRNA_ind) = p_junc->fwd_mRNA_hop_rate;
      (*p_H)(desc_mRNA_ind, parent_mRNA_ind) = p_junc->bkwd_mRNA_hop_rate;

      (*p_H)(parent_prot_ind, desc_prot_ind) = p_junc->fwd_prot_hop_rate;
      (*p_H)(desc_prot_ind, parent_prot_ind) = p_junc->bkwd_prot_hop_rate;
    }
    else {
      std::cerr << "-----------------------------------------\n"
                << "ERROR: UNKNOWN TYPE JUNCTION FOUND\n"
                << "-----------------------------------------\n";
      exit(1);
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    sem_set_hopping_rate_matrix(*((*it_p_junc)->p_to));
}

const Compartment* Analytic_engine::set_o1_soma() {
  initialise_o1_mat_and_RHS();
  
  auto& soma = *p_neuron->p_soma;
  auto& o1_mat = *p_o1_mat;

  o1_mat(0,0) = -soma.gene_activation_rate - soma.gene_deactivation_rate;

  o1_mat(1,0) = soma.transcription_rate;
  o1_mat(1,1) = -soma.mRNA_decay_rate;

  o1_mat(2,1) = soma.translation_rate;
  o1_mat(2,2) = -soma.protein_decay_rate;

  o1_var_names[0] = soma.name + "__Gene";
  o1_var_names[1] = soma.name + "__mRNA";
  o1_var_names[2] = soma.name + "__Prot";

  p_o1_vars[0] = &soma.n_active_genes_expectation;
  p_o1_vars[1] = &soma.n_mRNA_expectation;
  p_o1_vars[2] = &soma.n_prot_expectation;

  (*p_o1_RHS)(0) = -(p_neuron->p_soma->gene_activation_rate) * (p_neuron->p_soma->number_of_gene_copies);
  
  return p_neuron->p_soma;
}

void Analytic_engine::set_o1_matrix(const Compartment& parent) {

  auto& o1_mat = *p_o1_mat;
  
  if(parent.it_p_out_junctions.empty())
    return;

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;
    
    size_t& parent_start_ind = p_junc -> p_from -> o1_index;
    size_t& desc_start_ind = p_junc -> p_to -> o1_index;
    
    if(p_junc->type() == DEN_SYN) {
      o1_var_names[p_junc->p_to->o1_index] = p_junc->p_to->name + "__Prot";
      p_o1_vars[p_junc->p_to->o1_index] = &p_junc->p_to->n_prot_expectation;

      o1_mat(desc_start_ind, desc_start_ind) -= p_junc->p_to->protein_decay_rate;
      
      o1_mat(parent_start_ind+1, parent_start_ind+1) -= p_junc->fwd_prot_hop_rate;
      o1_mat(parent_start_ind+1, desc_start_ind) += p_junc->bkwd_prot_hop_rate;
      o1_mat(desc_start_ind, parent_start_ind+1) += p_junc->fwd_prot_hop_rate;
      o1_mat(desc_start_ind, desc_start_ind) -= p_junc->bkwd_prot_hop_rate;
    }
    else if(p_junc->type() == DEN_DEN) {
      o1_var_names[p_junc->p_to->o1_index] = p_junc->p_to->name + "__mRNA";
      o1_var_names[p_junc->p_to->o1_index+1] = p_junc->p_to->name + "__Prot";
      p_o1_vars[p_junc->p_to->o1_index] = &p_junc->p_to->n_mRNA_expectation;
      p_o1_vars[p_junc->p_to->o1_index+1] = &p_junc->p_to->n_prot_expectation;

      o1_mat(desc_start_ind, desc_start_ind) -= p_junc->p_to->mRNA_decay_rate;
      o1_mat(desc_start_ind+1, desc_start_ind) += p_junc->p_to->translation_rate;
      o1_mat(desc_start_ind+1, desc_start_ind+1) -= p_junc->p_to->protein_decay_rate;

      o1_mat(parent_start_ind, parent_start_ind) -= p_junc->fwd_mRNA_hop_rate;
      o1_mat(parent_start_ind, desc_start_ind) += p_junc->bkwd_mRNA_hop_rate;
      o1_mat(desc_start_ind, parent_start_ind) += p_junc->fwd_mRNA_hop_rate;
      o1_mat(desc_start_ind, desc_start_ind) -= p_junc->bkwd_mRNA_hop_rate;

      o1_mat(parent_start_ind+1, parent_start_ind+1) -= p_junc->fwd_prot_hop_rate;
      o1_mat(parent_start_ind+1, desc_start_ind+1) += p_junc->bkwd_prot_hop_rate;
      o1_mat(desc_start_ind+1, parent_start_ind+1) += p_junc->fwd_prot_hop_rate;
      o1_mat(desc_start_ind+1, desc_start_ind+1) -= p_junc->bkwd_prot_hop_rate;
    }
    else if(p_junc->type() == SOM_DEN) {
      o1_var_names[p_junc->p_to->o1_index] = p_junc->p_to->name + "__mRNA";
      o1_var_names[p_junc->p_to->o1_index+1] = p_junc->p_to->name + "__Prot";
      p_o1_vars[p_junc->p_to->o1_index] = &p_junc->p_to->n_mRNA_expectation;
      p_o1_vars[p_junc->p_to->o1_index+1] = &p_junc->p_to->n_prot_expectation;
      
      o1_mat(desc_start_ind, desc_start_ind) -= p_junc->p_to->mRNA_decay_rate;
      o1_mat(desc_start_ind+1, desc_start_ind) += p_junc->p_to->translation_rate;
      o1_mat(desc_start_ind+1, desc_start_ind+1) -= p_junc->p_to->protein_decay_rate;

      o1_mat(parent_start_ind+1, parent_start_ind+1) -= p_junc->fwd_mRNA_hop_rate;
      o1_mat(parent_start_ind+1, desc_start_ind) += p_junc->bkwd_mRNA_hop_rate;
      o1_mat(desc_start_ind, parent_start_ind+1) += p_junc->fwd_mRNA_hop_rate;
      o1_mat(desc_start_ind, desc_start_ind) -= p_junc->bkwd_mRNA_hop_rate;

      o1_mat(parent_start_ind+2, parent_start_ind+2) -= p_junc->fwd_prot_hop_rate;
      o1_mat(parent_start_ind+2, desc_start_ind+1) += p_junc->bkwd_prot_hop_rate;
      o1_mat(desc_start_ind+1, parent_start_ind+2) += p_junc->fwd_prot_hop_rate;
      o1_mat(desc_start_ind+1, desc_start_ind+1) -= p_junc->bkwd_prot_hop_rate;
    }
    else {
      std::cerr << "-----------------------------------------\n"
                << "ERROR: UNKNOWN TYPE JUNCTION FOUND\n"
                << "-----------------------------------------\n";
      exit(1);
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_o1_matrix(*((*it_p_junc)->p_to));
}

const Compartment* Analytic_engine::sem_set_o1_soma() {
  initialise_o1_mat_and_RHS();
  
  auto& soma = *p_neuron->p_soma;
  auto& o1_mat = *p_o1_mat;

  set_prot_index_from(soma);

  o1_mat(0,0) = -soma.gene_activation_rate - soma.gene_deactivation_rate;

  o1_mat(1,0) = soma.transcription_rate; //soma.mRNA_ind==1
  o1_mat(1,1) = -soma.mRNA_decay_rate;

  o1_mat(soma.prot_ind,1) = soma.translation_rate;
  o1_mat(soma.prot_ind,soma.prot_ind) = -soma.protein_decay_rate;

  o1_var_names[0] = soma.name + "__Gene";
  o1_var_names[1] = soma.name + "__mRNA";
  o1_var_names[soma.prot_ind] = soma.name + "__Prot";

  p_o1_vars[0] = &soma.n_active_genes_expectation;
  p_o1_vars[1] = &soma.n_mRNA_expectation;
  p_o1_vars[soma.prot_ind] = &soma.n_prot_expectation;

  (*p_o1_RHS)(0) = -(p_neuron->p_soma->gene_activation_rate) * (p_neuron->p_soma->number_of_gene_copies);
  
  return p_neuron->p_soma;
}

void Analytic_engine::sem_set_o1_matrix(const Compartment& parent) {
  auto& o1_mat = *p_o1_mat;
  
  if(parent.it_p_out_junctions.empty())
    return;

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;
    
    size_t& parent_prot_ind = p_junc -> p_from -> prot_ind;
    size_t& desc_prot_ind = p_junc -> p_to -> prot_ind;
    size_t& parent_mRNA_ind = p_junc -> p_from -> mRNA_ind;
    size_t& desc_mRNA_ind = p_junc -> p_to -> mRNA_ind;

    if(p_junc->type() == DEN_SYN) {
      o1_var_names[desc_prot_ind] = p_junc->p_to->name + "__Prot";
      p_o1_vars[desc_prot_ind] = &p_junc->p_to->n_prot_expectation;

      o1_mat(desc_prot_ind, desc_prot_ind) -= p_junc->p_to->protein_decay_rate;
      
      o1_mat(parent_prot_ind, parent_prot_ind) -= p_junc->fwd_prot_hop_rate;
      o1_mat(parent_prot_ind, desc_prot_ind) += p_junc->bkwd_prot_hop_rate;
      o1_mat(desc_prot_ind, parent_prot_ind) += p_junc->fwd_prot_hop_rate;
      o1_mat(desc_prot_ind, desc_prot_ind) -= p_junc->bkwd_prot_hop_rate;
    }
    else if(p_junc->type() == DEN_DEN) {
      o1_var_names[desc_mRNA_ind] = p_junc->p_to->name + "__mRNA";
      o1_var_names[desc_prot_ind] = p_junc->p_to->name + "__Prot";
      p_o1_vars[desc_mRNA_ind] = &p_junc->p_to->n_mRNA_expectation;
      p_o1_vars[desc_prot_ind] = &p_junc->p_to->n_prot_expectation;

      o1_mat(desc_mRNA_ind, desc_mRNA_ind) -= p_junc->p_to->mRNA_decay_rate;
      o1_mat(desc_prot_ind, desc_mRNA_ind) += p_junc->p_to->translation_rate;
      o1_mat(desc_prot_ind, desc_prot_ind) -= p_junc->p_to->protein_decay_rate;

      o1_mat(parent_mRNA_ind, parent_mRNA_ind) -= p_junc->fwd_mRNA_hop_rate;
      o1_mat(parent_mRNA_ind, desc_mRNA_ind) += p_junc->bkwd_mRNA_hop_rate;
      o1_mat(desc_mRNA_ind, parent_mRNA_ind) += p_junc->fwd_mRNA_hop_rate;
      o1_mat(desc_mRNA_ind, desc_mRNA_ind) -= p_junc->bkwd_mRNA_hop_rate;

      o1_mat(parent_prot_ind, parent_prot_ind) -= p_junc->fwd_prot_hop_rate;
      o1_mat(parent_prot_ind, desc_prot_ind) += p_junc->bkwd_prot_hop_rate;
      o1_mat(desc_prot_ind, parent_prot_ind) += p_junc->fwd_prot_hop_rate;
      o1_mat(desc_prot_ind, desc_prot_ind) -= p_junc->bkwd_prot_hop_rate;
    }
    else if(p_junc->type() == SOM_DEN) {
      o1_var_names[desc_mRNA_ind] = p_junc->p_to->name + "__mRNA";
      o1_var_names[desc_prot_ind] = p_junc->p_to->name + "__Prot";
      p_o1_vars[desc_mRNA_ind] = &p_junc->p_to->n_mRNA_expectation;
      p_o1_vars[desc_prot_ind] = &p_junc->p_to->n_prot_expectation;
      
      o1_mat(desc_mRNA_ind, desc_mRNA_ind) -= p_junc->p_to->mRNA_decay_rate;
      o1_mat(desc_prot_ind, desc_mRNA_ind) += p_junc->p_to->translation_rate;
      o1_mat(desc_prot_ind, desc_prot_ind) -= p_junc->p_to->protein_decay_rate;

      o1_mat(parent_mRNA_ind, parent_mRNA_ind) -= p_junc->fwd_mRNA_hop_rate;
      o1_mat(parent_mRNA_ind, desc_mRNA_ind) += p_junc->bkwd_mRNA_hop_rate;
      o1_mat(desc_mRNA_ind, parent_mRNA_ind) += p_junc->fwd_mRNA_hop_rate;
      o1_mat(desc_mRNA_ind, desc_mRNA_ind) -= p_junc->bkwd_mRNA_hop_rate;

      o1_mat(parent_prot_ind, parent_prot_ind) -= p_junc->fwd_prot_hop_rate;
      o1_mat(parent_prot_ind, desc_prot_ind) += p_junc->bkwd_prot_hop_rate;
      o1_mat(desc_prot_ind, parent_prot_ind) += p_junc->fwd_prot_hop_rate;
      o1_mat(desc_prot_ind, desc_prot_ind) -= p_junc->bkwd_prot_hop_rate;
    }
    else {
      std::cerr << "-----------------------------------------\n"
                << "ERROR: UNKNOWN TYPE JUNCTION FOUND\n"
                << "-----------------------------------------\n";
      exit(1);
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    sem_set_o1_matrix(*((*it_p_junc)->p_to));
}

Analytic_engine& Analytic_engine::internalise_expectations() {
  for(size_t i=0; i<o1_dim; ++i)
    *p_o1_vars[i] = expectations(i);
  return *this;
}

Analytic_engine& Analytic_engine::stationary_expectations() {
  // std::cerr << "Setting o1_matrix...\n";
  set_o1_matrix(*set_o1_soma());
  // std::cerr << "Inverting o1_matrix...\n";
  arma::mat inv_o1_matrix = (*p_o1_mat).i();
  // std::cerr << "Done with o1_matrix inversion\n";

  expectations = inv_o1_matrix * (*p_o1_RHS);

  // size_t i=0;
  // for(auto& o1_var_name : o1_var_names)
  //   std::cout << o1_var_name << "=\t" << expectations[i++] << std::endl;
  
  return internalise_expectations();
}

Analytic_engine& Analytic_engine::sem_stationary_expectations() {
  std::cerr << "Setting o1_matrix...\n";
  sem_set_o1_matrix(*sem_set_o1_soma());

  std::cerr << "Inverting o1_matrix...\n";
  arma::mat inv_o1_matrix = (*p_o1_mat).i();
  std::cerr << "Done with o1_matrix inversion\n";

  expectations = inv_o1_matrix * (*p_o1_RHS);

  size_t i=0;
  for(auto& o1_var_name : o1_var_names)
    std::cout << o1_var_name << "=\t" << expectations[i++] << std::endl;
  
  return internalise_expectations();
}

Analytic_engine& Analytic_engine::nonstationary_expectations(const std::list<double>& times) { //const Neuron& neur) {

  std::cerr << "Setting o1_matrix...\n";
  set_o1_matrix(*set_o1_soma());
  auto& o1_mat = *p_o1_mat;
  o1_mat = -o1_mat;
  
  std::cerr << "Computing eigen decomposition...\n";
  arma::cx_vec eigval_c;
  arma::vec eigval(o1_dim);
  arma::cx_mat eigvec_c;
  arma::mat tm(o1_dim, o1_dim); // Transition matrix
  arma::eig_gen(eigval_c, eigvec_c, o1_mat);
  
  for(size_t i=0; i<o1_dim; ++i) {
    eigval(i) = eigval_c(i).real();
    for(size_t j=0; j<o1_dim; ++j)
      tm(i,j) = eigvec_c(i,j).real();
  }

  std::cerr << "Inverting transition matrix...\n";
  arma::mat inv_tm = tm.i(); // Transition matrix
  
  std::cerr << "Eigenvalues:\n"
            << eigval << std::endl;
  std::cerr << "Eigenvectors:\n"
            << tm << std::endl;

  for(auto& o1_var_name : o1_var_names)
    std::cout << o1_var_name << ',';
  std::cout << std::endl;
  
  arma::vec stationary_expectations(o1_dim);
  // Computing stationary part
  for(size_t i=0; i<o1_dim; ++i) {
    double sum = 0;
    for(size_t k=0; k<o1_dim; ++k)
      sum += 1/eigval[k]*tm(i,k)*inv_tm(k,0);
    std::cout << "stationary_expectations[" << i << "]=" << (stationary_expectations(i) = (p_neuron->p_soma->gene_activation_rate)*(p_neuron->p_soma->number_of_gene_copies)*sum)  << std::endl;
  }

  // Setting integrals of motion c = G(0) - A*b, where G(0)=0 for now
  arma::vec c = 0 - stationary_expectations;

  arma::vec l_sum(o1_dim); // Precomputed l_sum
  for(size_t k=0; k<o1_dim; ++k)
    for(size_t l=0; l<o1_dim; ++l)
      l_sum(k) += inv_tm(k,l)*c(l);

  for(auto& t : times) {
    std::cout << t << ',';
    for(size_t i=0; i<o1_dim; ++i) {
      double k_sum = 0;
      for(size_t k=0; k<o1_dim; ++k)
        k_sum += exp(-eigval(k)*t) * tm(i,k) * l_sum(k);
      expectations(i) = stationary_expectations(i) + k_sum;
      std::cout << expectations(i) << ',';
    }
    std::cout << std::endl;
  }

  return internalise_expectations();
}


Analytic_engine& Analytic_engine::nonstationary_expectations(const double& t, const bool& reset_matrices, const bool& internalise) {
  
  if(reset_matrices) { // Setting o1 matrix
    set_o1_matrix(*set_o1_soma());
    set_As(*set_As_and_bs_soma());
  }
  
  // Computing eigen decomposition
  arma::cx_vec eigval_c;
  arma::vec eigval(o1_dim);
  arma::cx_mat eigvec_c;
  arma::mat tm(o1_dim, o1_dim); // Transition matrix
  arma::eig_gen(eigval_c, eigvec_c, (-*p_o1_mat));
  
  for(size_t i=0; i<o1_dim; ++i) {
    eigval(i) = eigval_c(i).real();
    for(size_t j=0; j<o1_dim; ++j)
      tm(i,j) = eigvec_c(i,j).real();
  }

  // Inverting transition matrix
  arma::mat inv_tm = tm.i(); // Transition matrix
        
  arma::vec stationary_expectations(o1_dim);
  // Computing stationary part
  for(size_t i=0; i<o1_dim; ++i) {
    double sum = 0;
    for(size_t k=0; k<o1_dim; ++k)
      sum += 1/eigval[k]*tm(i,k)*inv_tm(k,0);
    stationary_expectations(i) = (p_neuron->p_soma->gene_activation_rate)*(p_neuron->p_soma->number_of_gene_copies)*sum;
  }

  // Setting integrals of motion c = G(0) - A*b, where G(0)=0 for now
  arma::vec c = 0 - stationary_expectations;

  arma::vec l_sum(o1_dim); // Precomputed l_sum
  for(size_t k=0; k<o1_dim; ++k)
    for(size_t l=0; l<o1_dim; ++l)
      l_sum(k) += inv_tm(k,l)*c(l);

  
  for(size_t i=0; i<o1_dim; ++i) {
    double k_sum = 0;
    for(size_t k=0; k<o1_dim; ++k)
      k_sum += exp(-eigval(k)*t) * tm(i,k) * l_sum(k);
    expectations(i) = stationary_expectations(i) + k_sum;
  }
  
  return internalise ? internalise_expectations() : *this;
}


Analytic_engine& Analytic_engine::sem_nonstationary_expectations(const double& t, const bool& reset_matrices, const bool& internalise) { 
  
  if(reset_matrices) // Setting o1_matrix
    sem_set_o1_matrix(*sem_set_o1_soma());
  
  // Computing eigen decomposition
  arma::cx_vec eigval_c;
  arma::vec eigval(o1_dim);
  arma::cx_mat eigvec_c;
  arma::mat tm(o1_dim, o1_dim); // Transition matrix
  arma::eig_gen(eigval_c, eigvec_c, (-*p_o1_mat));
  
  for(size_t i=0; i<o1_dim; ++i) {
    eigval(i) = eigval_c(i).real();
    for(size_t j=0; j<o1_dim; ++j)
      tm(i,j) = eigvec_c(i,j).real();
  }

  // Inverting transition matrix
  arma::mat inv_tm = tm.i(); // Transition matrix
        
  arma::vec stationary_expectations(o1_dim);
  // Computing stationary part
  for(size_t i=0; i<o1_dim; ++i) {
    double sum = 0;
    for(size_t k=0; k<o1_dim; ++k)
      sum += 1/eigval[k]*tm(i,k)*inv_tm(k,0);
    stationary_expectations(i) = (p_neuron->p_soma->gene_activation_rate)*(p_neuron->p_soma->number_of_gene_copies)*sum;
  }

  // Setting integrals of motion c = G(0) - A*b, where G(0)=0 for now
  arma::vec c = 0 - stationary_expectations;

  arma::vec l_sum(o1_dim); // Precomputed l_sum
  for(size_t k=0; k<o1_dim; ++k)
    for(size_t l=0; l<o1_dim; ++l)
      l_sum(k) += inv_tm(k,l)*c(l);

  
  for(size_t i=0; i<o1_dim; ++i) {
    double k_sum = 0;
    for(size_t k=0; k<o1_dim; ++k)
      k_sum += exp(-eigval(k)*t) * tm(i,k) * l_sum(k);
    expectations(i) = stationary_expectations(i) + k_sum;
  }
  
  return internalise ? internalise_expectations() : *this;
}

Analytic_engine& Analytic_engine::sem_nonstationary_expectations(const std::list<double>& times) { //const Neuron& neur) {
  
  std::cerr << "Setting o1_matrix...\n";
  sem_set_o1_matrix(*sem_set_o1_soma());
  auto& o1_mat = *p_o1_mat;
  o1_mat = -o1_mat;

  std::cerr << "Computing eigen decomposition...\n";
  arma::cx_vec eigval_c;
  arma::vec eigval(o1_dim);
  arma::cx_mat eigvec_c;
  arma::mat tm(o1_dim, o1_dim); // Transition matrix
  arma::eig_gen(eigval_c, eigvec_c, o1_mat);
  
  for(size_t i=0; i<o1_dim; ++i) {
    eigval(i) = eigval_c(i).real();
    for(size_t j=0; j<o1_dim; ++j)
      tm(i,j) = eigvec_c(i,j).real();
  }

  std::cerr << "Inverting transition matrix...\n";
  arma::mat inv_tm = tm.i(); // Transition matrix
  
  std::cerr << "Eigenvalues:\n"
            << eigval << std::endl;
  std::cerr << "Eigenvectors:\n"
            << tm << std::endl;
    
  for(auto& o1_var_name : o1_var_names)
    std::cout << o1_var_name << ',';
  std::cout << std::endl;
  
  arma::vec stationary_expectations(o1_dim);
  // Computing stationary part
  for(size_t i=0; i<o1_dim; ++i) {
    double sum = 0;
    for(size_t k=0; k<o1_dim; ++k)
      sum += 1/eigval[k]*tm(i,k)*inv_tm(k,0);
    std::cout << "stationary_expectations[" << i << "]=" << (stationary_expectations(i) = (p_neuron->p_soma->gene_activation_rate)*(p_neuron->p_soma->number_of_gene_copies)*sum)  << std::endl;
  }

  // Setting integrals of motion c = G(0) - A*b, where G(0)=0 for now
  arma::vec c = 0 - stationary_expectations;

  arma::vec l_sum(o1_dim); // Precomputed l_sum
  for(size_t k=0; k<o1_dim; ++k)
    for(size_t l=0; l<o1_dim; ++l)
      l_sum(k) += inv_tm(k,l)*c(l);

  for(auto& t : times) {
    std::cout << t << ',';
    for(size_t i=0; i<o1_dim; ++i) {
      double k_sum = 0;
      for(size_t k=0; k<o1_dim; ++k)
        k_sum += exp(-eigval(k)*t) * tm(i,k) * l_sum(k);
      expectations(i) = stationary_expectations(i) + k_sum;
      std::cout << expectations(i) << ',';
    }
    std::cout << std::endl;
  }

  return internalise_expectations();
}

Analytic_engine& Analytic_engine::nonstationary_active_genes_expectations_direct_ODE_solver_step(const double& dt) {
  auto& soma = *p_neuron->p_soma;
  
  soma.n_active_genes_expectation += ((soma.number_of_gene_copies - soma.n_active_genes_expectation)*soma.gene_activation_rate - soma.n_active_genes_expectation*soma.gene_deactivation_rate)*dt;
  
  return *this;
}

Analytic_engine& Analytic_engine::nonstationary_mRNA_expectations_direct_ODE_solver_step(const double& dt, const bool& reset_matrices, const bool& internalise) {
    
  auto& soma = *p_neuron->p_soma;

  if(reset_matrices) {// Setting matrices
    set_mRNA_As(*set_mRNA_As_soma());
    std::cerr << "mRNA_Ap:\n" << *p_mRNA_Ap << std::endl
              << "mRNA_Am:\n" << *p_mRNA_Am << std::endl
              << "mRNA_expectations:\n" << expectations << std::endl;
  }

  mRNA_expectations += (*p_mRNA_Am)*mRNA_expectations*dt;
  mRNA_expectations(0) += soma.n_active_genes_expectation * soma.transcription_rate*dt;

  // std::cerr << (((*p_Am)*expectations + (*p_b))*dt).t() << '\n';
  
  return internalise ? internalise_mRNA_expectations() : *this;
}

Analytic_engine& Analytic_engine::nonstationary_protein_expectations_direct_ODE_solver_step(const double& dt, const bool& reset_matrices, const bool& internalise, const bool& source_update) {
    
  auto& soma = *p_neuron->p_soma;
  
  if(reset_matrices) {// Setting matrices
    set_prot_As(*set_prot_As_soma());
    std::cerr << "prot_Ap:\n" << *p_prot_Ap << std::endl
              << "prot_Am:\n" << *p_prot_Am << std::endl
              << "prot_expectations:\n" << expectations << std::endl;
  }
  else if(source_update)
    update_prot_source(*update_prot_source_soma());

  protein_expectations += ((*p_prot_Am)*protein_expectations - o1_prot_RHS)*dt;
    
  return internalise ? internalise_prot_expectations() : *this;
}


Analytic_engine& Analytic_engine::nonstationary_expectations_direct_ODE_solver_step(const double& dt, const bool& reset_matrices, const bool& internalise) {
  if(reset_matrices) {// Setting matrices
    set_As(*set_As_and_bs_soma());
    initialise_hopping_rate_matrix();
    set_hopping_rate_matrix(*p_neuron->p_soma);
    
    std::cerr << "Ap:\n" << *p_Ap << std::endl
              << "Am:\n" << *p_Am << std::endl
              << "b:\n" << (*p_b).t() << std::endl
              << "H:\n" << (*p_H) << std::endl
              << "expectations:\n" << expectations << std::endl;
  }

  expectations += ((*p_Am)*expectations + (*p_b))*dt;
  // std::cerr << (((*p_Am)*expectations + (*p_b))*dt).t() << '\n';
  
  return internalise ? internalise_expectations() : *this;
}

Analytic_engine& Analytic_engine::nonstationary_gene_gene_covariances_direct_ODE_solver_step(const double& dt, const bool& reset_matrices) {
  auto& soma = *p_neuron->p_soma;

  if(reset_matrices) {// Setting matrices
    set_mRNA_As(*set_mRNA_As_soma());
    std::cerr << "mRNA_Ap:\n" << *p_mRNA_Ap << std::endl
              << "mRNA_Am:\n" << *p_mRNA_Am << std::endl
              << "mRNA_expectations:\n" << mRNA_expectations << std::endl;
  }

  soma.n_active_genes_variance += (-2*(soma.gene_activation_rate + soma.gene_deactivation_rate)*soma.n_active_genes_variance + ((2*soma.number_of_gene_copies-1)*soma.gene_activation_rate + soma.gene_deactivation_rate)*soma.n_active_genes_expectation + soma.number_of_gene_copies*soma.gene_activation_rate)*dt;
    
  return *this;
}

Analytic_engine& Analytic_engine::nonstationary_gene_mRNA_covariances_direct_ODE_solver_step(const double& dt, const bool& reset) {
  
  auto& soma = *p_neuron->p_soma;
  size_t mRNA_dim = 1+p_neuron->p_dend_segments.size();                                         
  
  if(reset)
    o2_gene_mRNA = new arma::vec(mRNA_dim);
  
  (*o2_gene_mRNA) += (-(soma.gene_activation_rate+soma.gene_deactivation_rate)*(*o2_gene_mRNA) + (*p_mRNA_Am)*(*o2_gene_mRNA) + soma.number_of_gene_copies*soma.gene_activation_rate * mRNA_expectations) * dt;

  (*o2_gene_mRNA)(0) += soma.transcription_rate*soma.n_active_genes_variance*dt;
  
  return *this;
}

Analytic_engine& Analytic_engine::nonstationary_mRNA_mRNA_covariances_direct_ODE_solver_step(const double& dt, const bool& reset) {
  auto& soma = *p_neuron->p_soma;  

  size_t mRNA_dim = 1+p_neuron->p_dend_segments.size();
    
  if(reset) {
    if(p_mRNA_mRNA_cov_mat) delete p_mRNA_mRNA_cov_mat;
    p_mRNA_mRNA_cov_mat = new arma::mat(mRNA_dim, mRNA_dim);
    set_mRNA_As(*set_mRNA_As_soma());
    initialise_mRNA_hopping_rate_matrix();
    set_mRNA_hopping_rate_matrix(soma);
  }
  
  arma::mat MM = (*p_mRNA_Am)*(*p_mRNA_mRNA_cov_mat);
  (*p_mRNA_mRNA_cov_mat) += (MM + MM.t())*dt;

  (*p_mRNA_mRNA_cov_mat)(0,0) += soma.transcription_rate*(soma.n_active_genes_expectation + (*o2_gene_mRNA)(0))*dt;

  arma::vec V = (*p_mRNA_Ap)*mRNA_expectations;
  for(size_t i=0; i<mRNA_dim; ++i) {
    (*p_mRNA_mRNA_cov_mat)(i,i) += V(i)*dt;
    (*p_mRNA_mRNA_cov_mat)(i,0) = ((*p_mRNA_mRNA_cov_mat)(0,i) += soma.transcription_rate*(*o2_gene_mRNA)(i)*dt);
  }
      
  for(size_t i=0; i<mRNA_dim; ++i)
    for(size_t j=0; j<i; ++j)
      (*p_mRNA_mRNA_cov_mat)(j,i) = ((*p_mRNA_mRNA_cov_mat)(i,j) -= ((*p_mRNA_H)(i,j)*mRNA_expectations(i) + (*p_mRNA_H)(j,i)*mRNA_expectations(j))*dt);

  return *this;
}

Analytic_engine& Analytic_engine::nonstationary_prot_prot_covariances_direct_ODE_solver_step(const double& dt, const bool& reset) {
  auto& soma = *p_neuron->p_soma;  

  size_t prot_dim = 1+p_neuron->p_dend_segments.size()+p_neuron->p_synapses.size();
    
  if(reset) {
    if(p_prot_prot_cov_mat) delete p_prot_prot_cov_mat;
    p_prot_prot_cov_mat = new arma::mat(prot_dim, prot_dim);
    //    std::cerr << "p_prot_prot_cov_mat INITIALISATION:\n" << *p_prot_prot_cov_mat << std::endl;
    initialise_prot_hopping_rate_matrix();
    set_prot_hopping_rate_matrix(soma);
  }
  
  arma::mat PP = (*p_prot_Am)*(*p_prot_prot_cov_mat);
  (*p_prot_prot_cov_mat) += (PP + PP.t())*dt;

  arma::mat PM = (*p_PM)*(*p_mRNA_prot_cov_mat);
  (*p_prot_prot_cov_mat) += (PM + PM.t())*dt;


  (*p_mRNA_mRNA_cov_mat)(0,0) += soma.transcription_rate*(soma.n_active_genes_expectation + (*o2_gene_mRNA)(0))*dt;

  arma::vec V = (*p_PM)*mRNA_expectations + (*p_prot_Ap)*protein_expectations;
  for(size_t i=0; i<prot_dim; ++i)
    (*p_prot_prot_cov_mat)(i,i) += V(i)*dt;
      
  for(size_t i=0; i<prot_dim; ++i)
    for(size_t j=0; j<i; ++j)
      (*p_prot_prot_cov_mat)(j,i) = ((*p_prot_prot_cov_mat)(i,j) -= ((*p_prot_H)(i,j)*protein_expectations(i) + (*p_prot_H)(j,i)*protein_expectations(j))*dt);

  return *this;
}

Analytic_engine& Analytic_engine::stationary_expectations_and_correlations(const std::string& active_gene_expectation_file_name, const std::string& active_gene_variance_file_name, const std::string& mRNA_expectations_file_name, const std::string& protein_expectations_file_name, const std::string& gene_mRNA_covariances_file_name, const std::string& gene_prot_covariances_file_name, const std::string& mRNA_mRNA_covariances_file_name, const std::string& mRNA_prot_covariances_file_name, const std::string& prot_prot_covariances_file_name) {
  std::ofstream
    ofs_active_gene_expectations(active_gene_expectation_file_name),
    ofs_active_gene_variance(active_gene_variance_file_name),
    ofs_mRNA_expectations(mRNA_expectations_file_name),
    ofs_protein_expectations(protein_expectations_file_name),
    ofs_gene_mRNA_covariances(gene_mRNA_covariances_file_name),
    ofs_gene_prot_covariances(gene_prot_covariances_file_name),
    ofs_mRNA_covariances(mRNA_mRNA_covariances_file_name),
    ofs_mRNA_prot_covariances(mRNA_prot_covariances_file_name),
    ofs_prot_prot_covariances(prot_prot_covariances_file_name);

  auto& soma = *p_neuron->p_soma;
  size_t mRNA_dim = 1 + p_neuron->p_dend_segments.size(),
    prot_dim = mRNA_dim + p_neuron->p_synapses.size();

  ofs_active_gene_expectations << soma.n_active_genes_expectation;
  ofs_active_gene_expectations.close();
  ofs_active_gene_variance << soma.n_active_genes_variance;
  ofs_active_gene_variance.close();

  // Computing stationary expectations using matrix inversion
  mRNA_stationary_expectations().protein_stationary_expectations();

  ofs_mRNA_expectations << mRNA_expectations(0);
  for(size_t i=1; i<mRNA_dim; ++i)
    ofs_mRNA_expectations << ',' << mRNA_expectations(i);  
  ofs_mRNA_expectations.close();

  ofs_protein_expectations << protein_expectations(0);
  for(size_t i=1; i<prot_dim; ++i)
    ofs_protein_expectations << ',' << protein_expectations(i);  
  ofs_protein_expectations.close();

  // Computing gene-mRNA and gene-Protein covariances using matrix inversion
  gene_mRNA_stationary_covariances().gene_protein_stationary_covariances();

  ofs_gene_mRNA_covariances << (*o2_gene_mRNA)(0);
  for(size_t i=1; i<mRNA_dim; ++i)
    ofs_gene_mRNA_covariances << ',' << (*o2_gene_mRNA)(i);  
  ofs_gene_mRNA_covariances.close();

  ofs_gene_prot_covariances << (*o2_gene_prot)(0);
  for(size_t i=1; i<prot_dim; ++i)
    ofs_gene_prot_covariances << ',' << (*o2_gene_prot)(i);  
  ofs_gene_prot_covariances.close();

  // Computing the rest by solving the ODEs to convergence
  //// Determining timescales through eigen decomposition
  if(p_mRNA_mRNA_cov_mat) delete p_mRNA_mRNA_cov_mat;
  p_mRNA_mRNA_cov_mat = new arma::mat(mRNA_dim, mRNA_dim);
  set_mRNA_As(*set_mRNA_As_soma());
  initialise_mRNA_hopping_rate_matrix();
  set_mRNA_hopping_rate_matrix(soma);

  arma::cx_vec eigval_c;
  arma::cx_mat eigvec_c;
  arma::vec eigval(mRNA_dim);
  arma::vec eigval_min(2);
  arma::vec eigval_max(2);
  arma::eig_gen(eigval_c, eigvec_c, *p_mRNA_Am);
  for(size_t i=0; i<mRNA_dim; ++i)
    eigval(i) = abs(eigval_c(i).real());
  eigval_min(0) = arma::min(eigval);
  eigval_max(0) = arma::max(eigval);
  
  arma::eig_gen(eigval_c, eigvec_c, *p_mRNA_Ap);
  for(size_t i=0; i<mRNA_dim; ++i)
    eigval(i) = abs(eigval_c(i).real());
  eigval_min(1) = arma::min(eigval);
  eigval_max(1) = arma::max(eigval);

  // arma::eig_gen(eigval_c, eigvec_c, *p_mRNA_H);
  // for(size_t i=0; i<mRNA_dim; ++i)
  //   eigval(i) = abs(eigval_c(i).real());
  // eigval_min(2) = arma::min(eigval);
  // eigval_max(2) = arma::max(eigval);

  double t_start=0,
    t_fin= 10/arma::min(eigval_min),// 3 x (slowest timescale)
    dt = 1/(5*arma::max(eigval_max)); // 1/5 x (fastest timescale)
  
  std::cout << "Computing mRNA-mRNA covariances...\n";
  std::cout << "mRNA_Am_min_eig = " << eigval_min(0)
            << ";  mRNA_Am_max_eig = " << eigval_max(0) << std::endl
            << "mRNA_Ap_min_eig = " << eigval_min(1)
            << ";  mRNA_Ap_max_eig = " << eigval_max(1) << std::endl;
            // << "mRNA_H_min_eig = " << eigval_min(2)
            // << ";  mRNA_H_max_eig = " << eigval_max(2) << std::endl;
  std::cout << "t_fin = " << t_fin << ";  dt = " << dt << std::endl;

  for(double t=t_start; t<t_fin; t+=dt)
      nonstationary_mRNA_mRNA_covariances_direct_ODE_solver_step(dt);

  for(size_t i=0; i<mRNA_dim; ++i) {
    ofs_mRNA_covariances << (*p_mRNA_mRNA_cov_mat)(i,0);
    for(size_t j=1; j<mRNA_dim; ++j)
      ofs_mRNA_covariances << ',' << (*p_mRNA_mRNA_cov_mat)(i,j);
    ofs_mRNA_covariances << std::endl;
  }  
  ofs_mRNA_covariances.close();

  set_prot_As(*set_prot_As_soma());
  initialise_prot_hopping_rate_matrix();
  set_prot_hopping_rate_matrix(soma);
  set_PM(*set_PM_soma());
  
  arma::cx_vec eigval_c_prot;
  arma::cx_mat eigvec_c_prot;
  arma::vec eigval_prot(prot_dim);
  arma::vec eigval_min_prot(2);
  arma::vec eigval_max_prot(2);
  arma::eig_gen(eigval_c_prot, eigvec_c_prot, *p_prot_Am);
  for(size_t i=0; i<prot_dim; ++i)
    eigval_prot(i) = abs(eigval_c_prot(i).real());
  eigval_min_prot(0) = arma::min(eigval_prot);
  eigval_max_prot(0) = arma::max(eigval_prot);
  
  arma::eig_gen(eigval_c_prot, eigvec_c_prot, *p_prot_Ap);
  for(size_t i=0; i<prot_dim; ++i)
    eigval_prot(i) = abs(eigval_c_prot(i).real());
  eigval_min_prot(1) = arma::min(eigval_prot);
  eigval_max_prot(1) = arma::max(eigval_prot);

  double t_fin_prot= 3/arma::min(eigval_min_prot),// 3 x (slowest timescale)
    dt_prot = 1/(5*arma::max(eigval_max_prot)); // 1/5 x (fastest timescale)

  t_fin = std::max(t_fin, t_fin_prot);
  dt = std::min(dt, dt_prot);
  
  std::cout << "Computing mRNA-Protein covariances...\n";
  std::cout << "prot_Am_min_eig = " << eigval_min_prot(0)
            << ";  prot_Am_max_eig = " << eigval_max_prot(0) << std::endl
            << "prot_Ap_min_eig = " << eigval_min_prot(1)
            << ";  prot_Ap_max_eig = " << eigval_max_prot(1) << std::endl;
  // << "prot_H_min_eig = " << eigval_min(2)
  // << ";  prot_H_max_eig = " << eigval_max(2) << std::endl;
  std::cout << "t_fin = " << t_fin << ";  dt = " << dt << std::endl;

  if(o2_mRNA_prot) delete o2_mRNA_prot;
  p_mRNA_prot_cov_mat = new arma::mat(mRNA_dim, prot_dim);
  
  for(double t=t_start; t<t_fin; t+=dt)
    nonstationary_mRNA_prot_covariances_direct_ODE_solver_step(dt);

  for(size_t i=0; i<mRNA_dim; ++i) {
    ofs_mRNA_prot_covariances << (*p_mRNA_prot_cov_mat)(i,0);
    for(size_t j=1; j<prot_dim; ++j)
      ofs_mRNA_prot_covariances << ',' << (*p_mRNA_prot_cov_mat)(i,j);
    ofs_mRNA_prot_covariances << std::endl;
  }  
  ofs_mRNA_prot_covariances.close();
  
  std::cout << "Computing Protein-Protein covariances...\n";
  std::cout << "t_fin = " << t_fin_prot << ";  dt = " << dt << std::endl;
  if(p_prot_prot_cov_mat) delete p_prot_prot_cov_mat;
  p_prot_prot_cov_mat = new arma::mat(prot_dim, prot_dim);

  for(double t=t_start; t<t_fin_prot; t+=dt)
    nonstationary_prot_prot_covariances_direct_ODE_solver_step(dt);

  for(size_t i=0; i<mRNA_dim; ++i) {
    ofs_prot_prot_covariances << (*p_prot_prot_cov_mat)(i,0);
    for(size_t j=1; j<prot_dim; ++j)
      ofs_prot_prot_covariances << ',' << (*p_prot_prot_cov_mat)(i,j);
    ofs_prot_prot_covariances << std::endl;
  }
  ofs_prot_prot_covariances.close();
  
  return *this;
}

Analytic_engine& Analytic_engine::stationary_expectations_and_correlations(const double& dt_mRNA, const double& dt_prot, const double& t_fin_mRNA, const double& t_fin_prot, const std::string& active_gene_expectation_file_name, const std::string& active_gene_variance_file_name, const std::string& mRNA_expectations_file_name, const std::string& protein_expectations_file_name, const std::string& gene_mRNA_covariances_file_name, const std::string& gene_prot_covariances_file_name, const std::string& mRNA_mRNA_covariances_file_name, const std::string& mRNA_prot_covariances_file_name, const std::string& prot_prot_covariances_file_name) {
  std::ofstream
    ofs_active_gene_expectations(active_gene_expectation_file_name),
    ofs_active_gene_variance(active_gene_variance_file_name),
    ofs_mRNA_expectations(mRNA_expectations_file_name),
    ofs_protein_expectations(protein_expectations_file_name),
    ofs_gene_mRNA_covariances(gene_mRNA_covariances_file_name),
    ofs_gene_prot_covariances(gene_prot_covariances_file_name),
    ofs_mRNA_covariances(mRNA_mRNA_covariances_file_name),
    ofs_mRNA_prot_covariances(mRNA_prot_covariances_file_name),
    ofs_prot_prot_covariances(prot_prot_covariances_file_name);

  auto& soma = *p_neuron->p_soma;
  size_t mRNA_dim = 1 + p_neuron->p_dend_segments.size(),
    prot_dim = mRNA_dim + p_neuron->p_synapses.size();

  ofs_active_gene_expectations << soma.n_active_genes_expectation;
  ofs_active_gene_expectations.close();
  ofs_active_gene_variance << soma.n_active_genes_variance;
  ofs_active_gene_variance.close();

  // Computing stationary expectations using matrix inversion
  mRNA_stationary_expectations().protein_stationary_expectations();

  ofs_mRNA_expectations << mRNA_expectations(0);
  for(size_t i=1; i<mRNA_dim; ++i)
    ofs_mRNA_expectations << ',' << mRNA_expectations(i);  
  ofs_mRNA_expectations.close();

  ofs_protein_expectations << protein_expectations(0);
  for(size_t i=1; i<prot_dim; ++i)
    ofs_protein_expectations << ',' << protein_expectations(i);  
  ofs_protein_expectations.close();

  // Computing gene-mRNA and gene-Protein covariances using matrix inversion
  gene_mRNA_stationary_covariances().gene_protein_stationary_covariances();

  ofs_gene_mRNA_covariances << (*o2_gene_mRNA)(0);
  for(size_t i=1; i<mRNA_dim; ++i)
    ofs_gene_mRNA_covariances << ',' << (*o2_gene_mRNA)(i);  
  ofs_gene_mRNA_covariances.close();

  ofs_gene_prot_covariances << (*o2_gene_prot)(0);
  for(size_t i=1; i<prot_dim; ++i)
    ofs_gene_prot_covariances << ',' << (*o2_gene_prot)(i);  
  ofs_gene_prot_covariances.close();

  // Computing the rest by solving the ODEs to convergence
  if(p_mRNA_mRNA_cov_mat) delete p_mRNA_mRNA_cov_mat;
  p_mRNA_mRNA_cov_mat = new arma::mat(mRNA_dim, mRNA_dim);
  set_mRNA_As(*set_mRNA_As_soma());
  initialise_mRNA_hopping_rate_matrix();
  set_mRNA_hopping_rate_matrix(soma);

  std::cout << "Computing mRNA-mRNA covariances...\n";
  
  for(double t=0; t<t_fin_mRNA; t+=dt_mRNA)
      nonstationary_mRNA_mRNA_covariances_direct_ODE_solver_step(dt_mRNA);

  for(size_t i=0; i<mRNA_dim; ++i) {
    ofs_mRNA_covariances << (*p_mRNA_mRNA_cov_mat)(i,0);
    for(size_t j=1; j<mRNA_dim; ++j)
      ofs_mRNA_covariances << ',' << (*p_mRNA_mRNA_cov_mat)(i,j);
    ofs_mRNA_covariances << std::endl;
  }  
  ofs_mRNA_covariances.close();

  set_prot_As(*set_prot_As_soma());
  initialise_prot_hopping_rate_matrix();
  set_prot_hopping_rate_matrix(soma);
  set_PM(*set_PM_soma());
  
  double t_fin = std::max(t_fin_mRNA, t_fin_prot);
  double dt = std::min(dt_mRNA, dt_prot);
  
  std::cout << "Computing mRNA-Protein covariances...\n";

  if(o2_mRNA_prot) delete o2_mRNA_prot;
  p_mRNA_prot_cov_mat = new arma::mat(mRNA_dim, prot_dim);
  
  for(double t=0; t<t_fin; t+=dt)
    nonstationary_mRNA_prot_covariances_direct_ODE_solver_step(dt);

  for(size_t i=0; i<mRNA_dim; ++i) {
    ofs_mRNA_prot_covariances << (*p_mRNA_prot_cov_mat)(i,0);
    for(size_t j=1; j<prot_dim; ++j)
      ofs_mRNA_prot_covariances << ',' << (*p_mRNA_prot_cov_mat)(i,j);
    ofs_mRNA_prot_covariances << std::endl;
  }  
  ofs_mRNA_prot_covariances.close();
  
  std::cout << "Computing Protein-Protein covariances...\n";
  if(p_prot_prot_cov_mat) delete p_prot_prot_cov_mat;
  p_prot_prot_cov_mat = new arma::mat(prot_dim, prot_dim);

  for(double t=0; t<t_fin_prot; t+=dt)
    nonstationary_prot_prot_covariances_direct_ODE_solver_step(dt);

  for(size_t i=0; i<mRNA_dim; ++i) {
    ofs_prot_prot_covariances << (*p_prot_prot_cov_mat)(i,0);
    for(size_t j=1; j<prot_dim; ++j)
      ofs_prot_prot_covariances << ',' << (*p_prot_prot_cov_mat)(i,j);
    ofs_prot_prot_covariances << std::endl;
  }
  ofs_prot_prot_covariances.close();
  
  return *this;
}

Analytic_engine& Analytic_engine::nonstationary_gene_prot_covariances_direct_ODE_solver_step(const double& dt, const bool& reset) {

  auto& soma = *p_neuron->p_soma;
  size_t prot_dim = 1+p_neuron->p_dend_segments.size()+p_neuron->p_synapses.size();

  if(reset) {
    if(o2_gene_prot) delete o2_gene_prot;
    o2_gene_prot = new arma::vec(prot_dim);
    set_PM(*set_PM_soma());
  }
  
  (*o2_gene_prot) += ( -(soma.gene_activation_rate+soma.gene_deactivation_rate)*(*o2_gene_prot)
                       + (*p_prot_Am)*(*o2_gene_prot)
                       + (*p_PM)*(*o2_gene_mRNA)
                       + soma.gene_activation_rate*soma.number_of_gene_copies*protein_expectations )*dt;
  
  return *this;
}

Analytic_engine& Analytic_engine::nonstationary_mRNA_prot_covariances_direct_ODE_solver_step(const double& dt, const bool& reset) {

  auto& soma = *p_neuron->p_soma;
  size_t mRNA_dim = 1 + p_neuron->p_dend_segments.size(),
    prot_dim = mRNA_dim + p_neuron->p_synapses.size();

  if(reset) {
    if(o2_mRNA_prot) delete o2_mRNA_prot;
    p_mRNA_prot_cov_mat = new arma::mat(mRNA_dim, prot_dim);
    set_prot_As(*set_prot_As_soma());
    set_PM(*set_PM_soma());
  }
  
  (*p_mRNA_prot_cov_mat) += ( (*p_mRNA_Am)*(*p_mRNA_prot_cov_mat)
                              + ((*p_PM)*(*p_mRNA_mRNA_cov_mat)).t()
                              + (*p_mRNA_prot_cov_mat)*(*p_prot_Am).t() )*dt;

  for(size_t j=0; j<prot_dim; ++j)
    (*p_mRNA_prot_cov_mat)(0,j) += (soma.transcription_rate*(*o2_gene_prot)(j))*dt;
  
  return *this;
}

Analytic_engine& Analytic_engine::nonstationary_covariances_direct_ODE_solver_step(const double& dt, const bool& reset_matrices) {
  if(reset_matrices) {// Setting matrices
    set_As(*set_As_and_bs_soma());
    set_hopping_rate_matrix(*p_neuron->p_soma);
    std::cerr << "Ap:\n" << *p_Ap << std::endl
              << "Am:\n" << *p_Am << std::endl
              << "b:\n" << (*p_b).t() << std::endl
              << "H:\n" << (*p_H) << std::endl;
  }
  
  arma::mat M = (*p_Am)*(*p_cov_mat);
  (*p_cov_mat) += (M + M.t())*dt;

  // Computing diffusion matrix
  arma::vec Aps = (*p_Ap)*expectations;  
  (*p_cov_mat)(0,0) += ( 2*(*p_b)[0]*expectations[0] +  Aps[0] + (*p_b)[0] )*dt;
  for(size_t i=1; i<o1_dim; ++i) {
    (*p_cov_mat)(i,0) = ((*p_cov_mat)(0,i) += (*p_b)[0]*expectations[i]*dt);
    (*p_cov_mat)(i,i) += Aps[i]*dt;
  }
  for(size_t i=1; i<o1_dim; ++i)
    for(size_t j=1; j<i; ++j)
      (*p_cov_mat)(j,i) = ((*p_cov_mat)(i,j) -= ((*p_H)(i,j)*expectations[i] + (*p_H)(j,i)*expectations[j])*dt);

  return *this;
}

Analytic_engine& Analytic_engine::sem_nonstationary_expectations_direct_ODE_solver_step(const double& dt, const bool& reset_matrices, const bool& internalise) {
  if(reset_matrices) {// Setting matrices
    sem_set_As(*sem_set_As_and_bs_soma());
    initialise_hopping_rate_matrix();
    sem_set_hopping_rate_matrix(*p_neuron->p_soma);
    std::cerr << "Ap:\n" << *p_Ap << std::endl
              << "Am:\n" << *p_Am << std::endl
              << "b:\n" << (*p_b).t() << std::endl
              << "H:\n" << (*p_H) << std::endl
              << "expectations:\n" << expectations << std::endl;
  }

  expectations += ((*p_Am)*expectations + (*p_b))*dt;
  // std::cerr << (((*p_Am)*expectations + (*p_b))*dt).t() << '\n';
  
  return internalise ? internalise_expectations() : *this;
}

Analytic_engine& Analytic_engine::sem_nonstationary_covariances_direct_ODE_solver_step(const double& dt, const bool& reset_matrices) {
  if(reset_matrices) {// Setting matrices
    set_As(*set_As_and_bs_soma());
    set_hopping_rate_matrix(*p_neuron->p_soma);
    std::cerr << "Ap:\n" << *p_Ap << std::endl
              << "Am:\n" << *p_Am << std::endl
              << "b:\n" << (*p_b).t() << std::endl
              << "H:\n" << (*p_H) << std::endl;
  }
  
  arma::mat M = (*p_Am)*(*p_cov_mat);
  (*p_cov_mat) += (M + M.t())*dt;

  // Computing diffusion matrix
  arma::vec Aps = (*p_Ap)*expectations;  
  (*p_cov_mat)(0,0) += ( 2*(*p_b)[0]*expectations[0] +  Aps[0] + (*p_b)[0] )*dt;
  for(size_t i=1; i<o1_dim; ++i) {
    (*p_cov_mat)(i,0) = ((*p_cov_mat)(0,i) += (*p_b)[0]*expectations[i]*dt);
    (*p_cov_mat)(i,i) += Aps[i]*dt;
  }
  for(size_t i=1; i<o1_dim; ++i)
    for(size_t j=1; j<i; ++j)
      (*p_cov_mat)(j,i) = ((*p_cov_mat)(i,j) -= ((*p_H)(i,j)*expectations[i] + (*p_H)(j,i)*expectations[j])*dt);

  return *this;
}


Analytic_engine& Analytic_engine::nonstationary_covariances(const std::list<double>& times, arma::vec* initial_G1, arma::vec* initial_G2) {
  std::cerr << "Setting o1_matrix...\n";
  set_o1_matrix(*set_o1_soma());
  auto& o1_mat = *p_o1_mat;
  o1_mat = -o1_mat;
  
  std::cerr << "Computing eigen decomposition...\n";
  arma::cx_vec eigval_c;
  arma::vec eigval(o1_dim);
  arma::cx_mat eigvec_c;
  arma::mat tm(o1_dim, o1_dim); // Transition matrix
  arma::eig_gen(eigval_c, eigvec_c, o1_mat);
  
  for(size_t i=0; i<o1_dim; ++i) {
    eigval(i) = eigval_c(i).real();
    for(size_t j=0; j<o1_dim; ++j)
      tm(i,j) = eigvec_c(i,j).real();
  }

  std::cerr << "Inverting transition matrix...\n";
  arma::mat inv_tm = tm.i(); // Transition matrix
  
  std::cerr << "Eigenvalues:\n"
            << eigval << std::endl;
  std::cerr << "Eigenvectors:\n"
            << tm << std::endl;

  for(auto& o1_var_name : o1_var_names)
    std::cout << o1_var_name << ',';
  std::cout << std::endl;
  
  arma::vec stationary_expectations(o1_dim);
  // Computing stationary part
  for(size_t i=0; i<o1_dim; ++i) {
    double sum = 0;
    for(size_t k=0; k<o1_dim; ++k)
      sum += 1/eigval[k]*tm(i,k)*inv_tm(k,0);
    stationary_expectations(i) = p_neuron->p_soma->gene_activation_rate * p_neuron->p_soma->number_of_gene_copies * sum;
  }

  for(size_t i=0; i<o1_dim; ++i)
    std::cout << o1_var_names[i] << ": " << stationary_expectations(i) << std::endl;

  // Setting integrals of motion c = G(0) - A*b, where G(0)=0 for now
  arma::vec c_vec = *initial_G1 - stationary_expectations;

  // Initialising second order
  initialise_o2();
  set_o2_matrix();
  set_o2_nonstationary_RHS_mat();

  std::vector<double> k_sum(o1_dim); // to store precomputed k-sum
  
  auto& covariances = *p_covariances;
  auto& o2_mat = *p_o2_mat;
  o2_mat *= -1;

  for(size_t i=0; i<o1_dim; ++i)
    for(size_t j=0; j<o1_dim; ++j)
      (*o2_nonstationary_RHS_mat)(o2_ind(i,i),j) *= 2;

  for(size_t i=0; i<o1_dim; ++i)
    for(size_t j=0; j<o2_dim; ++j)
      o2_mat(o2_ind(i,i),j) *= 2;
  
  std::cerr << "Computing o2 eigen decomposition...\n";
  arma::cx_vec o2_eigval_c;
  arma::vec o2_eigval(o2_dim);
  arma::cx_mat o2_eigvec_c;
  arma::mat o2_tm(o2_dim, o2_dim); // Transition matrix
  arma::eig_gen(o2_eigval_c, o2_eigvec_c, o2_mat);

  for(size_t i=0; i<o2_dim; ++i) {
    o2_eigval(i) = o2_eigval_c(i).real();
    for(size_t j=0; j<o2_dim; ++j)
      o2_tm(i,j) = o2_eigvec_c(i,j).real();
  }

  std::cerr << "o2_eigenvalues:\n"
            << o2_eigval << std::endl;
            // << "o2_eigenvectors:\n"
            // << o2_tm << std::endl;

  std::cerr << "Inverting o2 transition matrix...\n";
  arma::mat o2_inv_tm = o2_tm.i(); // Transition matrix
  std::cerr << "o2 transition matrix inverted\n";
  
  for(auto& o2_var_name : *p_o2_var_names)
    std::cout << o2_var_name << ',';
  std::cout << std::endl;

  // Precomputing o1 matrix products
  arma::vec o1_l_sum(o1_dim);
  for(size_t k=0; k<o1_dim; ++k)
    for(size_t l=0; l<o1_dim; ++l)
      o1_l_sum(k) += inv_tm(k,l)*c_vec(l); 
  // Precomputing o2 matrix products
  arma::vec o2_s_sum(o2_dim), o2_eta_sum(o2_dim), o2_nu_sum(o1_dim);
  arma::mat o2_zeta_sum(o2_dim, o1_dim);
  for(size_t j=0; j<o2_dim; ++j) {    
    for(size_t s=0; s<o2_dim; ++s)
      o2_s_sum(j) += o2_inv_tm(j,s) * (*initial_G2)(s);

    arma::vec o2_k_sum(o1_dim);
    for(size_t zeta=0; zeta<o1_dim; ++zeta)
      for(size_t k=0; k<o2_dim; ++k)
        o2_k_sum(zeta) += o2_inv_tm(j,k)*(*o2_nonstationary_RHS_mat)(k,zeta);
    
    for(size_t eta=0; eta<o1_dim; ++eta) // computing o2_eta_sum
      o2_eta_sum(j) += o2_k_sum(eta)*stationary_expectations(eta);
    o2_eta_sum(j) /= o2_eigval(j);

    for(size_t mu=0; mu<o1_dim; ++mu) // computing o2_zeta_sum
      for(size_t zeta=0; zeta<o1_dim; ++zeta)
        o2_zeta_sum(j,mu) += o2_k_sum(zeta)*tm(zeta,mu);
  }
  for(size_t mu=0; mu<o1_dim; ++mu) // computing o2_nu_sum
    for(size_t nu=0; nu<o1_dim; ++nu)
      o2_nu_sum(mu) += inv_tm(mu,nu)*c_vec(nu);

  // main loop
  std::cout << "Main loop...\n";
  for(auto& t : times) {
    // o1
    for(size_t i=0; i<o1_dim; ++i) {
      double k_sum = 0;
      for(size_t k=0; k<o1_dim; ++k)
        k_sum += exp(-eigval(k)*t) * tm(i,k) * o1_l_sum(k);     
      expectations(i) = stationary_expectations(i) + k_sum;
    }
    // o2
    arma::vec braces(o2_dim);
    for(size_t j=0; j<o2_dim; ++j) {
      double o2_mu_sum = 0;
      for(size_t mu=0; mu<o1_dim; ++mu)
        if(o2_eigval(j) != eigval(mu))
          o2_mu_sum += (exp(-eigval(mu)*t)-exp(-o2_eigval(j)*t))/(o2_eigval(j)-eigval(mu)) * o2_zeta_sum(j,mu) * o2_nu_sum(mu);
        else
          o2_mu_sum += t*exp(-o2_eigval(j)*t) * o2_zeta_sum(j,mu) * o2_nu_sum(mu);
      braces(j) = exp(-o2_eigval(j)*t)*o2_s_sum(j) + (1-exp(-o2_eigval(j)*t))*o2_eta_sum(j) + o2_mu_sum;
    }
    for(size_t i=0; i<o2_dim; ++i) {
      (*p_covariances)[i] = 0;
      for(size_t j=0; j<o2_dim; ++j)
        (*p_covariances)[i] += o2_tm(i,j)*braces(j);
    }

    //////////// COMPUTING VARIANCES //////////
    std::vector<double> rmss(o1_dim);
    std::cout << t;
    for(size_t i=0; i<o1_dim; ++i) {
      rmss[i] = sqrt((*p_covariances)(o2_ind(i,i)) - expectations(i)*(expectations(i)-1));
      std::cout  << ',' << expectations(i) << ',' << rmss[i];
    }
    std::cout << std::endl;
    //////////////////////////////////////////
  }

  return internalise_expectations();
}

Analytic_engine& Analytic_engine::sem_nonstationary_covariances(const std::list<double>& times, arma::vec* initial_G1, arma::vec* initial_G2) {
  std::cout << "sem_nonstationary_covariances\n";
  std::cerr << "Setting o1_matrix...\n";
  sem_set_o1_matrix(*sem_set_o1_soma());
  auto& o1_mat = *p_o1_mat;
  o1_mat = -o1_mat;
    
  std::cerr << "Computing eigen decomposition...\n";
  arma::cx_vec eigval_c;
  arma::vec eigval(o1_dim);
  arma::cx_mat eigvec_c;
  arma::mat tm(o1_dim, o1_dim); // Transition matrix
  arma::eig_gen(eigval_c, eigvec_c, o1_mat);
  
  for(size_t i=0; i<o1_dim; ++i) {
    eigval(i) = eigval_c(i).real();
    for(size_t j=0; j<o1_dim; ++j)
      tm(i,j) = eigvec_c(i,j).real();
  }

  std::cerr << "Inverting transition matrix...\n";
  arma::mat inv_tm = tm.i(); // Transition matrix
  
  std::cerr << "Eigenvalues:\n"
            << eigval << std::endl;
  std::cerr << "Eigenvectors:\n"
            << tm << std::endl;

  for(auto& o1_var_name : o1_var_names)
    std::cout << o1_var_name << ',';
  std::cout << std::endl;
  
  arma::vec stationary_expectations(o1_dim);
  // Computing stationary part
  for(size_t i=0; i<o1_dim; ++i) {
    double sum = 0;
    for(size_t k=0; k<o1_dim; ++k)
      sum += 1/eigval[k]*tm(i,k)*inv_tm(k,0);
    stationary_expectations(i) = p_neuron->p_soma->gene_activation_rate * p_neuron->p_soma->number_of_gene_copies * sum;
  }

  std::cout << "Stationary expectations from nonstationary covariances algorithm:\n";
  for(size_t i=0; i<o1_dim; ++i)
    std::cout << o1_var_names[i] << ": " << stationary_expectations(i) << std::endl;

  // Setting integrals of motion c = G(0) - A*b
  arma::vec c_vec = *initial_G1 - stationary_expectations;
  
  std::cerr << "*initial_G1:\n" << *initial_G1 << std::endl;
  std::cerr << "*initial_G2:\n" << *initial_G2 << std::endl;
  std::cerr << "c_vec" << c_vec << std::endl;

  // Initialising second order
  sem_initialise_o2();
  sem_set_o2_matrix();
  sem_set_o2_nonstationary_RHS_mat();
      
  auto& covariances = *p_covariances;
  auto& o2_mat = *p_o2_mat;
  o2_mat *= -1;

  for(size_t i=0; i<o1_dim; ++i)
    for(size_t j=0; j<o1_dim; ++j)
      (*o2_nonstationary_RHS_mat)(sem_o2_ind(i,i),j) *= 2;

  for(size_t i=0; i<o1_dim; ++i)
    for(size_t j=0; j<o2_dim; ++j)
      o2_mat(sem_o2_ind(i,i),j) *= 2;


  std::cerr << "Computing o2 eigen decomposition...\n";
  arma::cx_vec o2_eigval_c;
  arma::vec o2_eigval(o2_dim);
  arma::cx_mat o2_eigvec_c;
  arma::mat o2_tm(o2_dim, o2_dim); // Transition matrix
  arma::eig_gen(o2_eigval_c, o2_eigvec_c, o2_mat);

  std::cerr << "o2_eigvec_c:\n" << o2_eigvec_c << std::endl;
 
  for(size_t i=0; i<o2_dim; ++i) {
    o2_eigval(i) = o2_eigval_c(i).real();
    for(size_t j=0; j<o2_dim; ++j)
      o2_tm(i,j) = o2_eigvec_c(i,j).real();
  }

  std::cerr << "o2_eigenvalues:\n"
            << o2_eigval << std::endl;
            // << "o2_eigenvectors:\n"
            // << o2_tm << std::endl;

  std::cerr << "Inverting o2 transition matrix...\n";
  arma::mat o2_inv_tm = o2_tm.i(); // Transition matrix
  std::cerr << "o2 transition matrix inverted\n";
  std::cerr << "Computing o2_tm eigen decomposition...\n";
  arma::eig_gen(o2_eigval_c, o2_eigvec_c, o2_tm);
  std::cerr << "o2_eigval_c:\n" << o2_eigval_c << std::endl;
  
  for(auto& o2_var_name : *p_o2_var_names)
    std::cout << o2_var_name << ',';
  std::cout << std::endl;

  // Precomputing
  // o1
  arma::vec o1_l_sum(o1_dim);
  for(size_t k=0; k<o1_dim; ++k)
    for(size_t l=0; l<o1_dim; ++l)
      o1_l_sum(k) += inv_tm(k,l)*c_vec(l); 
  // o2
  arma::vec o2_s_sum(o2_dim), o2_eta_sum(o2_dim), o2_nu_sum(o1_dim);
  arma::mat o2_zeta_sum(o2_dim, o1_dim);
  for(size_t j=0; j<o2_dim; ++j) {    
    for(size_t s=0; s<o2_dim; ++s)
      o2_s_sum(j) += o2_inv_tm(j,s) * (*initial_G2)(s);

    arma::vec o2_k_sum(o1_dim);
    for(size_t zeta=0; zeta<o1_dim; ++zeta)
      for(size_t k=0; k<o2_dim; ++k)
        o2_k_sum(zeta) += o2_inv_tm(j,k)*(*o2_nonstationary_RHS_mat)(k,zeta);
    
    for(size_t eta=0; eta<o1_dim; ++eta) // computing o2_eta_sum
      o2_eta_sum(j) += o2_k_sum(eta)*stationary_expectations(eta);
    o2_eta_sum(j) /= o2_eigval(j);

    for(size_t mu=0; mu<o1_dim; ++mu) // computing o2_zeta_sum
      for(size_t zeta=0; zeta<o1_dim; ++zeta)
        o2_zeta_sum(j,mu) += o2_k_sum(zeta)*tm(zeta,mu);
  }
  for(size_t mu=0; mu<o1_dim; ++mu) // computing o2_nu_sum
    for(size_t nu=0; nu<o1_dim; ++nu)
      o2_nu_sum(mu) += inv_tm(mu,nu)*c_vec(nu);

  // main loop
  std::cout << "Main loop...\n";
  for(auto& t : times) {
    // std::cout << "t=" << t << '\n';
    // o1
    for(size_t i=0; i<o1_dim; ++i) {
      double k_sum = 0;
      for(size_t k=0; k<o1_dim; ++k)
        k_sum += exp(-eigval(k)*t) * tm(i,k) * o1_l_sum(k);     
      expectations(i) = stationary_expectations(i) + k_sum;
    }
    // o2
    arma::vec braces(o2_dim);
    for(size_t j=0; j<o2_dim; ++j) {
      double o2_mu_sum = 0;
      for(size_t mu=0; mu<o1_dim; ++mu)
        if(o2_eigval(j) != eigval(mu))
          o2_mu_sum += (exp(-eigval(mu)*t)-exp(-o2_eigval(j)*t))/(o2_eigval(j)-eigval(mu)) * o2_zeta_sum(j,mu) * o2_nu_sum(mu);
        else
          o2_mu_sum += t*exp(-o2_eigval(j)*t) * o2_zeta_sum(j,mu) * o2_nu_sum(mu);
      braces(j) = exp(-o2_eigval(j)*t)*o2_s_sum(j) + (1-exp(-o2_eigval(j)*t))*o2_eta_sum(j) + o2_mu_sum;
    }
    for(size_t i=0; i<o2_dim; ++i) {
      (*p_covariances)[i] = 0;
      for(size_t j=0; j<o2_dim; ++j)
        (*p_covariances)[i] += o2_tm(i,j)*braces(j);
    }
    
    //////////// COMPUTING VARIANCES //////////
    std::vector<double> rmss(o1_dim);
    std::cout << t;
    for(size_t i=0; i<o1_dim; ++i) {
      rmss[i] = sqrt((*p_covariances)(sem_o2_ind(i,i)) - expectations(i)*(expectations(i)-1));
      std::cout  << ',' << expectations(i) << ',' << rmss[i];
    }
    std::cout << std::endl;
    //////////////////////////////////////////
  }

  return internalise_expectations();
}

Analytic_engine& Analytic_engine::sem_nonstationary_covariances_using_integral(const std::list<double>& times, arma::vec* initial_G1, arma::vec* initial_G2) {
  std::cout << "sem_nonstationary_covariances_using_integral\n";
  std::cerr << "Setting o1_matrix...\n";
  sem_set_o1_matrix(*sem_set_o1_soma());
  auto& o1_mat = *p_o1_mat;
  o1_mat = -o1_mat;
    
  std::cerr << "Computing eigen decomposition...\n";
  arma::cx_vec eigval_c;
  arma::vec eigval(o1_dim);
  arma::cx_mat eigvec_c;
  arma::mat tm(o1_dim, o1_dim); // Transition matrix
  arma::eig_gen(eigval_c, eigvec_c, o1_mat);
  
  for(size_t i=0; i<o1_dim; ++i) {
    eigval(i) = eigval_c(i).real();
    for(size_t j=0; j<o1_dim; ++j)
      tm(i,j) = eigvec_c(i,j).real();
  }

  std::cerr << "Inverting transition matrix...\n";
  arma::mat inv_tm = tm.i(); // Transition matrix
  
  std::cerr << "Eigenvalues:\n"
            << eigval << std::endl;
  std::cerr << "Eigenvectors:\n"
            << tm << std::endl;

  for(auto& o1_var_name : o1_var_names)
    std::cout << o1_var_name << ',';
  std::cout << std::endl;
  
  arma::vec stationary_expectations(o1_dim);
  // Computing stationary part
  for(size_t i=0; i<o1_dim; ++i) {
    double sum = 0;
    for(size_t k=0; k<o1_dim; ++k)
      sum += 1/eigval[k]*tm(i,k)*inv_tm(k,0);
    stationary_expectations(i) = p_neuron->p_soma->gene_activation_rate * p_neuron->p_soma->number_of_gene_copies * sum;
  }

  std::cout << "Stationary expectations from nonstationary covariances algorithm:\n";
  for(size_t i=0; i<o1_dim; ++i)
    std::cout << o1_var_names[i] << ": " << stationary_expectations(i) << std::endl;

  // Setting integrals of motion c = G(0) - A*b
  arma::vec c_vec = *initial_G1 - stationary_expectations;

  // Initialising second order
  sem_initialise_o2();
  sem_set_o2_matrix();
  sem_set_o2_nonstationary_RHS_mat();
      
  auto& covariances = *p_covariances;
  auto& o2_mat = *p_o2_mat;
  o2_mat *= -1;

  for(size_t i=0; i<o1_dim; ++i)
    for(size_t j=0; j<o1_dim; ++j)
      (*o2_nonstationary_RHS_mat)(sem_o2_ind(i,i),j) *= 2;

  for(size_t i=0; i<o1_dim; ++i)
    for(size_t j=0; j<o2_dim; ++j)
      o2_mat(sem_o2_ind(i,i),j) *= 2;

  std::cerr << "Computing o2 eigen decomposition...\n";
  arma::cx_vec o2_eigval_c;
  arma::vec o2_eigval(o2_dim);
  arma::cx_mat o2_eigvec_c;
  arma::mat o2_tm(o2_dim, o2_dim); // Transition matrix
  arma::eig_gen(o2_eigval_c, o2_eigvec_c, o2_mat);

  std::cerr << "o2_eigvec_c:\n" << o2_eigvec_c << std::endl;
 
  for(size_t i=0; i<o2_dim; ++i) {
    o2_eigval(i) = o2_eigval_c(i).real();
    for(size_t j=0; j<o2_dim; ++j)
      o2_tm(i,j) = o2_eigvec_c(i,j).real();
  }

  std::cerr << "o2_eigenvalues:\n"
            << o2_eigval << std::endl;
            // << "o2_eigenvectors:\n"
            // << o2_tm << std::endl;

  std::cerr << "Inverting o2 transition matrix...\n";
  arma::mat o2_inv_tm = o2_tm.i(); // Transition matrix
  std::cerr << "o2 transition matrix inverted\n";
  std::cerr << "Computing o2_tm eigen decomposition...\n";
  arma::eig_gen(o2_eigval_c, o2_eigvec_c, o2_tm);
  std::cerr << "o2_eigval_c:\n" << o2_eigval_c << std::endl;
  
  for(auto& o2_var_name : *p_o2_var_names)
    std::cout << o2_var_name << ',';
  std::cout << std::endl;

  // Precomputing o1 matrix products
  arma::vec o1_l_sum(o1_dim);
  for(size_t k=0; k<o1_dim; ++k)
    for(size_t l=0; l<o1_dim; ++l)
      o1_l_sum(k) += inv_tm(k,l)*c_vec(l); 
  // Precomputing o2 matrix products
  arma::vec o2_s_sum(o2_dim);
  arma::mat o2_k_sum(o2_dim, o1_dim);
  for(size_t j=0; j<o2_dim; ++j) {
    for(size_t s=0; s<o2_dim; ++s)
      o2_s_sum(j) += o2_inv_tm(j,s) * (*initial_G2)(s);

    for(size_t eta=0; eta<o1_dim; ++eta)
      for(size_t k=0; k<o2_dim; ++k)
        o2_k_sum(j,eta) += o2_inv_tm(j,k)*(*o2_nonstationary_RHS_mat)(k,eta);
  }

  // main loop
  arma::mat integral(o2_dim, o1_dim); // to store the integral over time
  std::vector<double> braces(o2_dim);
  std::cout << "Main loop...\n";
  double t_prev = times.front();
  for(auto& t : times) {
    // std::cout << "t=" << t << '\n';
    // o1
    for(size_t i=0; i<o1_dim; ++i) {
      double k_sum = 0;
      for(size_t k=0; k<o1_dim; ++k)
        k_sum += exp(-eigval(k)*t) * tm(i,k) * o1_l_sum(k);
      expectations(i) = stationary_expectations(i) + k_sum;
    }
    // o2
    for(size_t j=0; j<o2_dim; ++j) {
      double eta_sum=0;
      for(size_t eta=0; eta<o1_dim; ++eta) { // Precomputing the integral for all eta
        integral(j,eta) += (expectations(eta) - o2_eigval(j)*integral(j,eta))*(t-t_prev);
        eta_sum += o2_k_sum(j,eta)*integral(j,eta);
      }
      braces[j] = exp(-o2_eigval(j)*t)*o2_s_sum(j) + eta_sum;
    }
    for(size_t i=0; i<o2_dim; ++i) {
      (*p_covariances)(i) = 0;
      for(size_t j=0; j<o2_dim; ++j)        
        (*p_covariances)(i) += o2_tm(i,j)*braces[j];
    }
    t_prev = t;

    //////////// COMPUTING VARIANCES //////////
    std::vector<double> rmss(o1_dim);
    std::cout << t;
    for(size_t i=0; i<o1_dim; ++i) {
      rmss[i] = sqrt((*p_covariances)(sem_o2_ind(i,i)) - expectations(i)*(expectations(i)-1));
      // if(rmss[i]>=0) {
      //   std::cout << o1_var_names[i] + ": " << expectations(i) << ", " << rmss[i] << std::endl;
      // }
      // else
      //   std::cout << o1_var_names[i] + ": " << expectations(i) << ", " << rmss[i] << " NEGATIVE!\n";
      std::cout  << ',' << expectations(i) << ',' << rmss[i];
    }
    std::cout << std::endl;
    //////////////////////////////////////////    
  }

  return internalise_expectations();
}


Analytic_engine& Analytic_engine::sem_nonstationary_covariances_direct_ODE_solver(const std::list<double>& times, arma::vec* initial_G1, arma::vec* initial_G2) {
  std::cout << "sem_nonstationary_covariances_using_direct_ODE_solver\n";
  std::cerr << "Setting o1_matrix...\n";
  sem_set_o1_matrix(*sem_set_o1_soma());
  auto& o1_mat = *p_o1_mat;
  o1_mat = -o1_mat;
    
  std::cerr << "Computing eigen decomposition...\n";
  arma::cx_vec eigval_c;
  arma::vec eigval(o1_dim);
  arma::cx_mat eigvec_c;
  arma::mat tm(o1_dim, o1_dim); // Transition matrix
  arma::eig_gen(eigval_c, eigvec_c, o1_mat);
  
  for(size_t i=0; i<o1_dim; ++i) {
    eigval(i) = eigval_c(i).real();
    for(size_t j=0; j<o1_dim; ++j)
      tm(i,j) = eigvec_c(i,j).real();
  }

  std::cerr << "Inverting transition matrix...\n";
  arma::mat inv_tm = tm.i(); // Transition matrix
  
  std::cerr << "Eigenvalues:\n"
            << eigval << std::endl;
  std::cerr << "Eigenvectors:\n"
            << tm << std::endl;

  for(auto& o1_var_name : o1_var_names)
    std::cout << o1_var_name << ',';
  std::cout << std::endl;
  
  arma::vec stationary_expectations(o1_dim);
  // Computing stationary part
  for(size_t i=0; i<o1_dim; ++i) {
    double sum = 0;
    for(size_t k=0; k<o1_dim; ++k)
      sum += 1/eigval[k]*tm(i,k)*inv_tm(k,0);
    stationary_expectations(i) = p_neuron->p_soma->gene_activation_rate * p_neuron->p_soma->number_of_gene_copies * sum;
  }

  std::cout << "Stationary expectations from nonstationary covariances algorithm:\n";
  for(size_t i=0; i<o1_dim; ++i)
    std::cout << o1_var_names[i] << ": " << stationary_expectations(i) << std::endl;

  // Setting integrals of motion c = G(0) - A*b
  arma::vec c_vec = *initial_G1 - stationary_expectations;

  // Initialising second order
  sem_initialise_o2();
  sem_set_o2_matrix();
  sem_set_o2_nonstationary_RHS_mat();
      
  auto& covariances = *p_covariances;
  auto& o2_mat = *p_o2_mat;
  o2_mat *= -1;

  for(size_t i=0; i<o1_dim; ++i)
    for(size_t j=0; j<o1_dim; ++j)
      (*o2_nonstationary_RHS_mat)(sem_o2_ind(i,i),j) *= 2;

  for(size_t i=0; i<o1_dim; ++i)
    for(size_t j=0; j<o2_dim; ++j)
      o2_mat(sem_o2_ind(i,i),j) *= 2;

  
  for(auto& o2_var_name : *p_o2_var_names)
    std::cout << o2_var_name << ',';
  std::cout << std::endl;

  // Precomputing o1 matrix products
  arma::vec o1_l_sum(o1_dim);
  for(size_t k=0; k<o1_dim; ++k)
    for(size_t l=0; l<o1_dim; ++l)
      o1_l_sum(k) += inv_tm(k,l)*c_vec(l); 

  // main loop
  std::cout << "Main loop...\n";
  double t_prev = times.front();
  for(auto& t : times) {
    // std::cout << "t=" << t << '\n';
    // o1
    for(size_t i=0; i<o1_dim; ++i) {
      double k_sum = 0;
      for(size_t k=0; k<o1_dim; ++k)
        k_sum += exp(-eigval(k)*t) * tm(i,k) * o1_l_sum(k);
      expectations(i) = stationary_expectations(i) + k_sum;
    }
    // o2

    (*p_covariances) += ((*o2_nonstationary_RHS_mat) * expectations - (*p_o2_mat)*(*p_covariances))*(t-t_prev);
    
    t_prev = t;

    //////////// COMPUTING VARIANCES //////////
    std::vector<double> rmss(o1_dim);
    std::cout << t;
    for(size_t i=0; i<o1_dim; ++i) {
      rmss[i] = sqrt((*p_covariances)(sem_o2_ind(i,i)) - expectations(i)*(expectations(i)-1));
      std::cout  << ',' << expectations(i) << ',' << rmss[i];
    }
    std::cout << std::endl;
    //////////////////////////////////////////
  }

  return internalise_expectations();
}

Analytic_engine& Analytic_engine::sem_nonstationary_covariances_direct_ODE_solver_no_D_matrix(const std::list<double>& times, arma::vec* initial_G1, arma::vec* initial_G2) {
  std::cout << "sem_nonstationary_covariances_using_direct_ODE_solver\n";
  std::cerr << "Setting o1_matrix...\n";
  sem_set_o1_matrix(*sem_set_o1_soma());
  auto& o1_mat = *p_o1_mat;
  o1_mat = -o1_mat;
    
  std::cerr << "Computing eigen decomposition...\n";
  arma::cx_vec eigval_c;
  arma::vec eigval(o1_dim);
  arma::cx_mat eigvec_c;
  arma::mat tm(o1_dim, o1_dim); // Transition matrix
  arma::eig_gen(eigval_c, eigvec_c, o1_mat);
  
  for(size_t i=0; i<o1_dim; ++i) {
    eigval(i) = eigval_c(i).real();
    for(size_t j=0; j<o1_dim; ++j)
      tm(i,j) = eigvec_c(i,j).real();
  }

  std::cerr << "Inverting transition matrix...\n";
  arma::mat inv_tm = tm.i(); // Transition matrix
  
  std::cerr << "Eigenvalues:\n"
            << eigval << std::endl;
  std::cerr << "Eigenvectors:\n"
            << tm << std::endl;

  for(auto& o1_var_name : o1_var_names)
    std::cout << o1_var_name << ',';
  std::cout << std::endl;
  
  arma::vec stationary_expectations(o1_dim);
  // Computing stationary part
  for(size_t i=0; i<o1_dim; ++i) {
    double sum = 0;
    for(size_t k=0; k<o1_dim; ++k)
      sum += 1/eigval[k]*tm(i,k)*inv_tm(k,0);
    stationary_expectations(i) = p_neuron->p_soma->gene_activation_rate * p_neuron->p_soma->number_of_gene_copies * sum;
  }

  std::cout << "Stationary expectations from nonstationary covariances algorithm:\n";
  for(size_t i=0; i<o1_dim; ++i)
    std::cout << o1_var_names[i] << ": " << stationary_expectations(i) << std::endl;

  // Setting integrals of motion c = G(0) - A*b
  arma::vec c_vec = *initial_G1 - stationary_expectations;

  // Initialising second order
  sem_initialise_o2();
  sem_set_o2_matrix();
      
  auto& covariances = *p_covariances;
  auto& o2_mat = *p_o2_mat;
  o2_mat *= -1;

  for(size_t i=0; i<o1_dim; ++i)
    for(size_t j=0; j<o2_dim; ++j)
      o2_mat(sem_o2_ind(i,i),j) *= 2;  
  
  for(auto& o2_var_name : *p_o2_var_names)
    std::cout << o2_var_name << ',';
  std::cout << std::endl;

  // Precomputing o1 matrix products
  arma::vec o1_l_sum(o1_dim);
  for(size_t k=0; k<o1_dim; ++k)
    for(size_t l=0; l<o1_dim; ++l)
      o1_l_sum(k) += inv_tm(k,l)*c_vec(l); 

  // main loop
  std::cout << "Main loop...\n";
  double t_prev = times.front();
  for(auto& t : times) {
    // std::cout << "t=" << t << '\n';
    // o1
    for(size_t i=0; i<o1_dim; ++i) {
      double k_sum = 0;
      for(size_t k=0; k<o1_dim; ++k)
        k_sum += exp(-eigval(k)*t) * tm(i,k) * o1_l_sum(k);
      expectations(i) = stationary_expectations(i) + k_sum;
    }
    // o2
    for(size_t i=0; i<o2_dim; ++i)
      (*p_o2_RHS)[i] = 0;
    sem_set_o2_RHS();
    for(size_t i=0; i<o1_dim; ++i)
      (*p_o2_RHS)(sem_o2_ind(i,i)) *= 2;
    (*p_covariances) += (-(*p_o2_RHS) - o2_mat*(*p_covariances))*(t-t_prev);
    
    t_prev = t;

    //////////// COMPUTING VARIANCES //////////
    std::vector<double> rmss(o1_dim);
    std::cout << t;
    for(size_t i=0; i<o1_dim; ++i) {
      rmss[i] = sqrt((*p_covariances)(o2_ind(i,i)) - expectations(i)*(expectations(i)-1));
      std::cout  << ',' << expectations(i) << ',' << rmss[i];
    }
    std::cout << std::endl;
    //////////////////////////////////////////
  }

  return internalise_expectations();
}

size_t Analytic_engine::o2_ind(const size_t &i, const size_t &j, const size_t &dim) const {
  if(i<=j)
    return (2*dim-i-1)*i/2+j;
  else
    return o2_ind(j, i, dim);
}

size_t Analytic_engine::sem_o2_ind(const size_t &i, const size_t &j) const {
  if(i<=j) {
    size_t n_dend = p_neuron->p_dend_segments.size(),
      n_p = 1 + n_dend + p_neuron->p_synapses.size();
    if(i==0)
      return j;
    else if(i>0 && i<n_dend+2 && j<n_dend+2)
      return o1_dim + (2+2*n_dend-i)*(i-1)/2 + j-1;
    else if(i>0 && i<n_dend+2 && j>n_dend+1)
      return (n_dend+2)*(n_dend+1)/2 + n_p*i + j;
    else
      return n_p*(n_dend+2) + (n_dend+2)*(n_dend+1)/2 + (2*n_p+n_dend-i+1)*(i-n_dend-2)/2 + j;
    }
  else
    return sem_o2_ind(j, i);
}

void Analytic_engine::initialise_o2() {
  p_o2_mat = new arma::mat(o2_dim, o2_dim);
  p_o2_RHS =  new arma::vec(o2_dim);
  p_covariances = new arma::vec(o2_dim);
  
  p_o2_var_names = new std::vector<std::string>(o2_dim);
  std::vector<std::string> &o2_var_names = *p_o2_var_names;
  for(size_t i=0; i<o1_dim; ++i)
    for(size_t j=i; j<o1_dim; ++j)
      o2_var_names[o2_ind(i,j)] = o1_var_names[i] + '-' + o1_var_names[j];
}

void Analytic_engine::set_prot_index_from(Compartment& compartment) {
  if(compartment.type() == SOMA)
    compartment.prot_ind = (p_neuron->prot_ind = p_neuron->p_dend_segments.size()+2)++;
  else
    compartment.prot_ind = p_neuron->prot_ind++;

  for (auto& p_d_comp : compartment.p_descendants)  
    set_prot_index_from(*p_d_comp);
}

void Analytic_engine::sem_initialise_o2() {
  set_prot_index_from(*p_neuron->p_soma);
  
  if(p_o2_mat) clear_o2_mat_and_RHS();
  p_o2_mat = new arma::mat(o2_dim, o2_dim);
  p_o2_RHS =  new arma::vec(o2_dim);
  if(!p_covariances)
    p_covariances = new arma::vec(o2_dim);
  if(!p_o2_var_names)
    p_o2_var_names = new std::vector<std::string>(o2_dim);

  // sem_set_expectations(*sem_set_soma());

  std::vector<std::string> &o2_var_names = *p_o2_var_names;
  for(size_t i=0; i<o1_dim; ++i)
    for(size_t j=i; j<o1_dim; ++j)
      o2_var_names[sem_o2_ind(i,j)] = o1_var_names[i] + '-' + o1_var_names[j];
}

const Compartment* Analytic_engine::sem_set_soma() {
  auto& soma = *p_neuron->p_soma;

  o1_var_names[0] = soma.name + "__Gene";
  o1_var_names[1] = soma.name + "__mRNA";
  o1_var_names[soma.prot_ind] = soma.name + "__Prot";

  expectations(0) = soma.n_active_genes_expectation;
  expectations(1) = soma.n_mRNA_expectation;
  expectations(soma.prot_ind) = soma.n_prot_expectation;
  
  return &soma;
}

void Analytic_engine::sem_set_expectations(const Compartment& parent) {

  if(parent.it_p_out_junctions.empty())
    return;

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& desc = *(*it_p_junc)->p_to;
        
    if(desc.type() == SPINE) {
      o1_var_names[desc.prot_ind] = desc.name + "__Prot";

      expectations(desc.prot_ind) = desc.n_prot_expectation;
    }
    else if(desc.type() == APICAL_DENDRITE || desc.type() == BASAL_DENDRITE) {
      o1_var_names[desc.mRNA_ind] = desc.name + "__mRNA";
      o1_var_names[desc.prot_ind] = desc.name + "__Prot";

      expectations(desc.mRNA_ind) = desc.n_mRNA_expectation;
      expectations(desc.prot_ind) = desc.n_prot_expectation;
    }
    else {
      std::cerr << "-----------------------------------------\n"
                << "ERROR: UNKNOWN TYPE JUNCTION FOUND\n"
                << "-----------------------------------------\n";
      exit(1);
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    sem_set_expectations(*((*it_p_junc)->p_to));
}

void Analytic_engine::set_o2_nonstationary_RHS_soma() {

  o2_nonstationary_RHS_mat = new arma::mat(o2_dim, o1_dim);
  
  const auto& soma = *p_neuron->p_soma;

  (*o2_nonstationary_RHS_mat)(o2_ind(0,0), 0) -= soma.gene_activation_rate;
  (*o2_nonstationary_RHS_mat)(o2_ind(0,1), 0) += soma.transcription_rate;
  (*o2_nonstationary_RHS_mat)(o2_ind(1,2), 1) += soma.translation_rate;
  
  for(size_t i=0; i<o1_dim; ++i) 
    (*o2_nonstationary_RHS_mat)(o2_ind(0,i), i) += soma.gene_activation_rate*soma.number_of_gene_copies;
}

void Analytic_engine::sem_set_o2_nonstationary_RHS_soma() {

  o2_nonstationary_RHS_mat = new arma::mat(o2_dim, o1_dim);
  
  const auto& soma = *p_neuron->p_soma;

  (*o2_nonstationary_RHS_mat)(sem_o2_ind(0,0), 0) -= soma.gene_activation_rate;
  (*o2_nonstationary_RHS_mat)(sem_o2_ind(0,1), 0) += soma.transcription_rate;
  (*o2_nonstationary_RHS_mat)(sem_o2_ind(1,soma.prot_ind), 1) += soma.translation_rate;
  
  for(size_t i=0; i<o1_dim; ++i) 
    (*o2_nonstationary_RHS_mat)(sem_o2_ind(0,i), i) += soma.gene_activation_rate*soma.number_of_gene_copies;
}

void Analytic_engine::set_o2_soma() {
  const auto& soma = *p_neuron->p_soma;
  auto& o2_mat = *p_o2_mat;
  
  for(size_t i=0; i<o1_dim; ++i) {
    o2_mat(o2_ind(0,i), o2_ind(0,i)) -= (soma.gene_activation_rate + soma.gene_deactivation_rate);
    o2_mat(o2_ind(1,i), o2_ind(0,i)) += soma.transcription_rate;
    o2_mat(o2_ind(1,i), o2_ind(1,i)) -= soma.mRNA_decay_rate;
    o2_mat(o2_ind(2,i), o2_ind(1,i)) += soma.translation_rate;
    o2_mat(o2_ind(2,i), o2_ind(2,i)) -= soma.protein_decay_rate;
  }
}

void Analytic_engine::set_o2_nonstationary_RHS_mat() {
  set_o2_nonstationary_RHS_soma();

  for(auto& p_junc : p_neuron->p_junctions) { // Looping on junctions, setting descendants

    size_t& desc_start_ind = p_junc -> p_to -> o1_index;

    if(p_junc->type() != DEN_SYN)
      (*o2_nonstationary_RHS_mat)(o2_ind(desc_start_ind, desc_start_ind+1), desc_start_ind) += p_junc->p_to->translation_rate;
  }
}

void Analytic_engine::sem_set_o2_nonstationary_RHS_mat() {
  sem_set_o2_nonstationary_RHS_soma();

  for(auto& p_junc : p_neuron->p_junctions) { // Looping on junctions, setting descendants

    size_t& desc_mRNA_ind = p_junc -> p_to -> mRNA_ind;
    size_t& desc_prot_ind = p_junc -> p_to -> prot_ind;

    if(p_junc->type() != DEN_SYN)
      (*o2_nonstationary_RHS_mat)(sem_o2_ind(desc_mRNA_ind, desc_prot_ind), desc_mRNA_ind) += p_junc->p_to->translation_rate;
  }
}

void Analytic_engine::sem_set_o2_soma() {
  const auto& soma = *p_neuron->p_soma;
  auto& o2_mat = *p_o2_mat;
  auto& prot_ind = p_neuron->p_soma->prot_ind;
  
  for(size_t i=0; i<o1_dim; ++i) {
    o2_mat(sem_o2_ind(0,i), sem_o2_ind(0,i)) -= (soma.gene_activation_rate + soma.gene_deactivation_rate);
    o2_mat(sem_o2_ind(1,i), sem_o2_ind(0,i)) += soma.transcription_rate;
    o2_mat(sem_o2_ind(1,i), sem_o2_ind(1,i)) -= soma.mRNA_decay_rate;
    o2_mat(sem_o2_ind(prot_ind,i), sem_o2_ind(1,i)) += soma.translation_rate;
    o2_mat(sem_o2_ind(prot_ind,i), sem_o2_ind(prot_ind,i)) -= soma.protein_decay_rate;
  }
}

void Analytic_engine::set_o2_matrix() {
  set_o2_soma();

  auto& o2_mat = *p_o2_mat;
  
  for(auto& p_junc : p_neuron->p_junctions) { // Looping on junctions, setting descendants

    size_t& parent_start_ind = p_junc -> p_from -> o1_index;
    size_t& desc_start_ind = p_junc -> p_to -> o1_index;

    if(p_junc->type() == DEN_SYN)
      for(size_t i=0; i<o1_dim; ++i) {
        o2_mat(o2_ind(parent_start_ind+1, i), o2_ind(parent_start_ind+1, i)) -= p_junc->fwd_prot_hop_rate;
        o2_mat(o2_ind(desc_start_ind, i), o2_ind(parent_start_ind+1, i)) += p_junc->fwd_prot_hop_rate;
        o2_mat(o2_ind(desc_start_ind, i), o2_ind(desc_start_ind, i)) -= p_junc->bkwd_prot_hop_rate;
        o2_mat(o2_ind(parent_start_ind+1, i), o2_ind(desc_start_ind, i)) += p_junc->bkwd_prot_hop_rate;

        o2_mat(o2_ind(desc_start_ind, i), o2_ind(desc_start_ind, i)) -= p_junc->p_to->protein_decay_rate;
      }
    else if(p_junc->type() == DEN_DEN)
      for(size_t i=0; i<o1_dim; ++i) {
        o2_mat(o2_ind(parent_start_ind, i), o2_ind(parent_start_ind, i)) -= p_junc->fwd_mRNA_hop_rate;
        o2_mat(o2_ind(desc_start_ind, i), o2_ind(parent_start_ind, i)) += p_junc->fwd_mRNA_hop_rate;
        o2_mat(o2_ind(desc_start_ind, i), o2_ind(desc_start_ind, i)) -= p_junc->bkwd_mRNA_hop_rate;
        o2_mat(o2_ind(parent_start_ind, i), o2_ind(desc_start_ind, i)) += p_junc->bkwd_mRNA_hop_rate;       
        
        o2_mat(o2_ind(parent_start_ind+1, i), o2_ind(parent_start_ind+1, i)) -= p_junc->fwd_prot_hop_rate;
        o2_mat(o2_ind(desc_start_ind+1, i), o2_ind(parent_start_ind+1, i)) += p_junc->fwd_prot_hop_rate;
        o2_mat(o2_ind(desc_start_ind+1, i), o2_ind(desc_start_ind+1, i)) -= p_junc->bkwd_prot_hop_rate;
        o2_mat(o2_ind(parent_start_ind+1, i), o2_ind(desc_start_ind+1, i)) += p_junc->bkwd_prot_hop_rate;

        o2_mat(o2_ind(desc_start_ind+1, i), o2_ind(desc_start_ind, i)) += p_junc->p_to->translation_rate;
        o2_mat(o2_ind(desc_start_ind+1, i), o2_ind(desc_start_ind+1, i)) -= p_junc->p_to->protein_decay_rate;
        o2_mat(o2_ind(desc_start_ind, i), o2_ind(desc_start_ind, i)) -= p_junc->p_to->mRNA_decay_rate;
      }
    else if(p_junc->type() == SOM_DEN)
      for(size_t i=0; i<o1_dim; ++i) {
        o2_mat(o2_ind(parent_start_ind+1, i), o2_ind(parent_start_ind+1, i)) -= p_junc->fwd_mRNA_hop_rate;
        o2_mat(o2_ind(desc_start_ind, i), o2_ind(parent_start_ind+1, i)) += p_junc->fwd_mRNA_hop_rate;
        o2_mat(o2_ind(desc_start_ind, i), o2_ind(desc_start_ind, i)) -= p_junc->bkwd_mRNA_hop_rate;
        o2_mat(o2_ind(parent_start_ind+1, i), o2_ind(desc_start_ind, i)) += p_junc->bkwd_mRNA_hop_rate;       
        
        o2_mat(o2_ind(parent_start_ind+2, i), o2_ind(parent_start_ind+2, i)) -= p_junc->fwd_prot_hop_rate;
        o2_mat(o2_ind(desc_start_ind+1, i), o2_ind(parent_start_ind+2, i)) += p_junc->fwd_prot_hop_rate;
        o2_mat(o2_ind(desc_start_ind+1, i), o2_ind(desc_start_ind+1, i)) -= p_junc->bkwd_prot_hop_rate;
        o2_mat(o2_ind(parent_start_ind+2, i), o2_ind(desc_start_ind+1, i)) += p_junc->bkwd_prot_hop_rate;

        o2_mat(o2_ind(desc_start_ind+1, i), o2_ind(desc_start_ind, i)) += p_junc->p_to->translation_rate;
        o2_mat(o2_ind(desc_start_ind+1, i), o2_ind(desc_start_ind+1, i)) -= p_junc->p_to->protein_decay_rate;
        o2_mat(o2_ind(desc_start_ind, i), o2_ind(desc_start_ind, i)) -= p_junc->p_to->mRNA_decay_rate;
      }
    else {
      std::cerr << "-----------------------------------------\n"
                << "ERROR: UNKNOWN TYPE JUNCTION FOUND\n"
                << "-----------------------------------------\n";
      exit(1);
    }
  }
}

void Analytic_engine::sem_set_o2_matrix() {
  sem_set_o2_soma();

  auto& o2_mat = *p_o2_mat;
  
  for(auto& p_junc : p_neuron->p_junctions) { // Looping on junctions, setting descendants

    size_t& par_mRNA_ind = p_junc -> p_from -> mRNA_ind;
    size_t& par_prot_ind = p_junc -> p_from -> prot_ind;
    size_t& desc_mRNA_ind = p_junc -> p_to -> mRNA_ind;
    size_t& desc_prot_ind = p_junc -> p_to -> prot_ind;

    if(p_junc->type() == DEN_SYN)
      for(size_t i=0; i<o1_dim; ++i) {
        o2_mat(sem_o2_ind(par_prot_ind, i), sem_o2_ind(par_prot_ind, i)) -= p_junc->fwd_prot_hop_rate;
        o2_mat(sem_o2_ind(desc_prot_ind, i), sem_o2_ind(par_prot_ind, i)) += p_junc->fwd_prot_hop_rate;
        o2_mat(sem_o2_ind(desc_prot_ind, i), sem_o2_ind(desc_prot_ind, i)) -= p_junc->bkwd_prot_hop_rate;
        o2_mat(sem_o2_ind(par_prot_ind, i), sem_o2_ind(desc_prot_ind, i)) += p_junc->bkwd_prot_hop_rate;

        o2_mat(sem_o2_ind(desc_prot_ind, i), sem_o2_ind(desc_prot_ind, i)) -= p_junc->p_to->protein_decay_rate;
      }
    else if(p_junc->type() == DEN_DEN)
      for(size_t i=0; i<o1_dim; ++i) {
        o2_mat(sem_o2_ind(par_mRNA_ind, i), sem_o2_ind(par_mRNA_ind, i)) -= p_junc->fwd_mRNA_hop_rate;
        o2_mat(sem_o2_ind(desc_mRNA_ind, i), sem_o2_ind(par_mRNA_ind, i)) += p_junc->fwd_mRNA_hop_rate;
        o2_mat(sem_o2_ind(desc_mRNA_ind, i), sem_o2_ind(desc_mRNA_ind, i)) -= p_junc->bkwd_mRNA_hop_rate;
        o2_mat(sem_o2_ind(par_mRNA_ind, i), sem_o2_ind(desc_mRNA_ind, i)) += p_junc->bkwd_mRNA_hop_rate;       
        
        o2_mat(sem_o2_ind(par_prot_ind, i), sem_o2_ind(par_prot_ind, i)) -= p_junc->fwd_prot_hop_rate;
        o2_mat(sem_o2_ind(desc_prot_ind, i), sem_o2_ind(par_prot_ind, i)) += p_junc->fwd_prot_hop_rate;
        o2_mat(sem_o2_ind(desc_prot_ind, i), sem_o2_ind(desc_prot_ind, i)) -= p_junc->bkwd_prot_hop_rate;
        o2_mat(sem_o2_ind(par_prot_ind, i), sem_o2_ind(desc_prot_ind, i)) += p_junc->bkwd_prot_hop_rate;

        o2_mat(sem_o2_ind(desc_prot_ind, i), sem_o2_ind(desc_mRNA_ind, i)) += p_junc->p_to->translation_rate;
        o2_mat(sem_o2_ind(desc_prot_ind, i), sem_o2_ind(desc_prot_ind, i)) -= p_junc->p_to->protein_decay_rate;
        o2_mat(sem_o2_ind(desc_mRNA_ind, i), sem_o2_ind(desc_mRNA_ind, i)) -= p_junc->p_to->mRNA_decay_rate;
      }
    else if(p_junc->type() == SOM_DEN)
      for(size_t i=0; i<o1_dim; ++i) {
        o2_mat(sem_o2_ind(par_mRNA_ind, i), sem_o2_ind(par_mRNA_ind, i)) -= p_junc->fwd_mRNA_hop_rate;
        o2_mat(sem_o2_ind(desc_mRNA_ind, i), sem_o2_ind(par_mRNA_ind, i)) += p_junc->fwd_mRNA_hop_rate;
        o2_mat(sem_o2_ind(desc_mRNA_ind, i), sem_o2_ind(desc_mRNA_ind, i)) -= p_junc->bkwd_mRNA_hop_rate;
        o2_mat(sem_o2_ind(par_mRNA_ind, i), sem_o2_ind(desc_mRNA_ind, i)) += p_junc->bkwd_mRNA_hop_rate;       
        
        o2_mat(sem_o2_ind(par_prot_ind, i), sem_o2_ind(par_prot_ind, i)) -= p_junc->fwd_prot_hop_rate;
        o2_mat(sem_o2_ind(desc_prot_ind, i), sem_o2_ind(par_prot_ind, i)) += p_junc->fwd_prot_hop_rate;
        o2_mat(sem_o2_ind(desc_prot_ind, i), sem_o2_ind(desc_prot_ind, i)) -= p_junc->bkwd_prot_hop_rate;
        o2_mat(sem_o2_ind(par_prot_ind, i), sem_o2_ind(desc_prot_ind, i)) += p_junc->bkwd_prot_hop_rate;

        o2_mat(sem_o2_ind(desc_prot_ind, i), sem_o2_ind(desc_mRNA_ind, i)) += p_junc->p_to->translation_rate;
        o2_mat(sem_o2_ind(desc_prot_ind, i), sem_o2_ind(desc_prot_ind, i)) -= p_junc->p_to->protein_decay_rate;
        o2_mat(sem_o2_ind(desc_mRNA_ind, i), sem_o2_ind(desc_mRNA_ind, i)) -= p_junc->p_to->mRNA_decay_rate;
      }
    else {
      std::cerr << "-----------------------------------------\n"
                << "ERROR: UNKNOWN TYPE JUNCTION FOUND\n"
                << "-----------------------------------------\n";
      exit(1);
    }
  }
}  


void Analytic_engine::set_o2_RHS() {
  auto& o2_RHS = *p_o2_RHS;
  const auto& soma = *p_neuron->p_soma;

  o2_RHS(o2_ind(0,0)) += soma.gene_activation_rate*expectations(0);
  o2_RHS(o2_ind(0,1)) -= soma.transcription_rate*expectations(0);
  o2_RHS(o2_ind(1,2)) -= soma.translation_rate*expectations(1);
  for(size_t i=0; i<o1_dim; ++i)
    o2_RHS(o2_ind(0,i)) -= soma.gene_activation_rate*soma.number_of_gene_copies*expectations(i);

  for(auto& ds : p_neuron->p_dend_segments)
    o2_RHS(o2_ind(ds->o1_index,ds->o1_index+1)) -= ds->translation_rate*expectations(ds->o1_index);
}

void Analytic_engine::sem_set_o2_RHS() {
  auto& o2_RHS = *p_o2_RHS;
  const auto& soma = *p_neuron->p_soma;

  o2_RHS(sem_o2_ind(0,0)) += soma.gene_activation_rate*expectations(0);
  o2_RHS(sem_o2_ind(0,1)) -= soma.transcription_rate*expectations(0);
  o2_RHS(sem_o2_ind(1,soma.prot_ind)) -= soma.translation_rate*expectations(1);
  for(size_t i=0; i<o1_dim; ++i)
    o2_RHS(sem_o2_ind(0,i)) -= soma.gene_activation_rate*soma.number_of_gene_copies*expectations(i);

  for(auto& ds : p_neuron->p_dend_segments)
    o2_RHS(sem_o2_ind(ds->mRNA_ind,ds->prot_ind)) -= ds->translation_rate*expectations(ds->mRNA_ind);
}

const Compartment* Analytic_engine::set_o2_gene_mRNA_soma() {
  Soma& soma = *(p_neuron->p_soma);
  size_t sz = 1+p_neuron->p_dend_segments.size();
  
  o2_gene_mRNA_RHS = new arma::vec(sz);
  o2_gene_mRNA_mat = new arma::mat(sz,sz);

  (*o2_gene_mRNA_RHS)(0) -= soma.gene_activation_rate*soma.number_of_gene_copies*soma.n_mRNA_expectation + soma.transcription_rate*(soma.n_active_genes_expectation + o2_gene_gene);

  (*o2_gene_mRNA_mat)(0,0) -= soma.gene_activation_rate + soma.gene_deactivation_rate + soma.mRNA_decay_rate;

  return &soma;
}

void Analytic_engine::set_o2_gene_mRNA_matrix(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty())
    return;

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;
    
    size_t parent_ind = p_junc -> p_from -> mRNA_ind-1;
    size_t desc_ind = p_junc -> p_to -> mRNA_ind-1;
    
    if(p_junc->type() != DEN_SYN) {

      (*o2_gene_mRNA_mat)(desc_ind, desc_ind) -= p_neuron->p_soma->gene_activation_rate + p_neuron->p_soma->gene_deactivation_rate + p_junc->p_to->mRNA_decay_rate;

      (*o2_gene_mRNA_mat)(parent_ind, parent_ind) -= p_junc->fwd_mRNA_hop_rate;
      (*o2_gene_mRNA_mat)(desc_ind, parent_ind) += p_junc->fwd_mRNA_hop_rate;
      (*o2_gene_mRNA_mat)(desc_ind, desc_ind) -= p_junc->bkwd_mRNA_hop_rate;
      (*o2_gene_mRNA_mat)(parent_ind, desc_ind) += p_junc->bkwd_mRNA_hop_rate;

      (*o2_gene_mRNA_RHS)(desc_ind) -= p_neuron->p_soma->gene_activation_rate * p_neuron->p_soma->number_of_gene_copies * p_junc->p_to->n_mRNA_expectation;
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_o2_gene_mRNA_matrix(*((*it_p_junc)->p_to));
}

Analytic_engine& Analytic_engine::gene_mRNA_stationary_covariances() {

  set_o2_gene_mRNA_matrix(*set_o2_gene_mRNA_soma());
  
  o2_gene_mRNA = new arma::vec((*o2_gene_mRNA_mat).i()*(*o2_gene_mRNA_RHS));

  delete o2_gene_mRNA_mat; o2_gene_mRNA_mat=nullptr;
  delete o2_gene_mRNA_RHS; o2_gene_mRNA_RHS=nullptr;
  
  return *this;
}

std::vector<double> Analytic_engine::stationary_gene_mRNA_covariances() {

  gene_mRNA_stationary_covariances();

  size_t dim = 1+p_neuron->p_dend_segments.size();
  
  std::vector<double> correlations(dim);
  for(size_t i=0; i<dim; ++i)
    correlations[i] = (*o2_gene_mRNA)(i);

  return correlations;
}

const Compartment* Analytic_engine::set_o2_gene_prot_soma() {
  Soma& soma = *(p_neuron->p_soma);
  size_t sz = 1+p_neuron->p_dend_segments.size()+p_neuron->p_synapses.size();
  
  o2_gene_prot_RHS = new arma::vec(sz);
  o2_gene_prot_mat = new arma::mat(sz,sz);

  (*o2_gene_prot_RHS)(0) -= soma.gene_activation_rate*soma.number_of_gene_copies*soma.n_prot_expectation + soma.translation_rate*(*o2_gene_mRNA)(0);

  (*o2_gene_prot_mat)(0,0) -= soma.gene_activation_rate + soma.gene_deactivation_rate + soma.protein_decay_rate;

  return &soma;
}

void Analytic_engine::set_o2_gene_prot_matrix(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty())
    return;

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;
    
    size_t parent_ind = p_junc -> p_from -> id;
    size_t desc_ind = p_junc -> p_to -> id;
    

    (*o2_gene_prot_mat)(desc_ind, desc_ind) -= p_neuron->p_soma->gene_activation_rate + p_neuron->p_soma->gene_deactivation_rate + p_junc->p_to->protein_decay_rate;

    (*o2_gene_prot_mat)(parent_ind, parent_ind) -= p_junc->fwd_prot_hop_rate;
    (*o2_gene_prot_mat)(desc_ind, parent_ind) += p_junc->fwd_prot_hop_rate;
    (*o2_gene_prot_mat)(desc_ind, desc_ind) -= p_junc->bkwd_prot_hop_rate;
    (*o2_gene_prot_mat)(parent_ind, desc_ind) += p_junc->bkwd_prot_hop_rate;
    
    if(p_junc->type() != DEN_SYN)
      (*o2_gene_prot_RHS)(desc_ind) -= p_neuron->p_soma->gene_activation_rate * p_neuron->p_soma->number_of_gene_copies * p_junc->p_to->n_prot_expectation + p_junc->p_to->translation_rate*(*o2_gene_mRNA)(p_junc->p_to->mRNA_ind-1);
    else
      (*o2_gene_prot_RHS)(desc_ind) -= p_neuron->p_soma->gene_activation_rate * p_neuron->p_soma->number_of_gene_copies * p_junc->p_to->n_prot_expectation;
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_o2_gene_prot_matrix(*((*it_p_junc)->p_to));
}

Analytic_engine& Analytic_engine::gene_protein_stationary_covariances() {

  set_o2_gene_prot_matrix(*set_o2_gene_prot_soma());
  
  o2_gene_prot = new arma::vec((*o2_gene_prot_mat).i()*(*o2_gene_prot_RHS));

  delete o2_gene_prot_mat; o2_gene_prot_mat=nullptr;
  delete o2_gene_prot_RHS; o2_gene_prot_RHS=nullptr;

  return *this;
}

std::vector<double> Analytic_engine::stationary_gene_protein_covariances() {

  gene_protein_stationary_covariances();

  size_t dim = p_neuron->p_dend_segments.size() + p_neuron->p_synapses.size() + 1;  
  std::vector<double> correlations(dim);
  for(size_t i=0; i<dim; ++i)
    correlations[i] = (*o2_gene_prot)(i);

  return correlations;
}

const Compartment* Analytic_engine::set_o2_mRNA_mRNA_soma() {
  Soma& soma = *(p_neuron->p_soma);
  size_t sz = 1+p_neuron->p_dend_segments.size();
  
  o2_mRNA_mRNA_RHS = new arma::vec(sz*(sz+1)/2);
  o2_mRNA_mRNA_mat = new arma::mat(sz*(sz+1)/2, sz*(sz+1)/2);

  for(size_t i=0; i<sz; ++i) {
    (*o2_mRNA_mRNA_RHS)(i) -= soma.transcription_rate*(*o2_gene_mRNA)(i);
    (*o2_mRNA_mRNA_mat)(i,i) -= soma.mRNA_decay_rate;
  }

  return &soma;
}

void Analytic_engine::set_o2_mRNA_mRNA_matrix(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty())
    return;

  size_t sz = 1 + p_neuron->p_dend_segments.size();

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;

    if(p_junc->type() != DEN_SYN) {
      
      size_t parent_ind = p_junc -> p_from -> mRNA_ind-1;
      size_t desc_ind = p_junc -> p_to -> mRNA_ind-1;
    
      for(size_t i=0; i<sz; ++i) {
        (*o2_mRNA_mRNA_mat)(o2_ind(desc_ind, i, sz), o2_ind(desc_ind, i, sz)) -= p_junc->p_to->mRNA_decay_rate;

        (*o2_mRNA_mRNA_mat)(o2_ind(parent_ind, i, sz), o2_ind(parent_ind, i, sz)) -= p_junc->fwd_mRNA_hop_rate;
        (*o2_mRNA_mRNA_mat)(o2_ind(desc_ind, i, sz), o2_ind(parent_ind, i, sz)) += p_junc->fwd_mRNA_hop_rate;
        (*o2_mRNA_mRNA_mat)(o2_ind(desc_ind, i, sz), o2_ind(desc_ind, i, sz)) -= p_junc->bkwd_mRNA_hop_rate;
        (*o2_mRNA_mRNA_mat)(o2_ind(parent_ind, i, sz), o2_ind(desc_ind, i, sz)) += p_junc->bkwd_mRNA_hop_rate;
      }
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_o2_mRNA_mRNA_matrix(*((*it_p_junc)->p_to));
}

Analytic_engine& Analytic_engine::mRNA_mRNA_stationary_covariances() {
  if(!(o2_gene_mRNA)) {
    std::cerr << "ERROR in computing mRNA-mRNA covariances:\n"
              << "Prior to computing mRNA-mRNA covariances coompute\n"
              << "mRNA expectations and gene-mRNA stationary covariances()\n";
    exit(1);
  }
    
  set_o2_mRNA_mRNA_matrix(*set_o2_mRNA_mRNA_soma());
  
  o2_mRNA_mRNA = new arma::vec((*o2_mRNA_mRNA_mat).i()*(*o2_mRNA_mRNA_RHS));

  delete o2_mRNA_mRNA_mat; o2_mRNA_mRNA_mat=nullptr;
  delete o2_mRNA_mRNA_RHS; o2_mRNA_mRNA_RHS=nullptr;

  return *this;
}

std::vector<std::vector<double>> Analytic_engine::stationary_mRNA_mRNA_covariances() {

  mRNA_mRNA_stationary_covariances();

  size_t dim = p_neuron->p_dend_segments.size() + 1;  
  std::vector<std::vector<double>> correlations(dim, std::vector<double>(dim));
  for(size_t i=0; i<dim; ++i)
    for(size_t j=0; j<i; ++j)
      correlations[i][j] = correlations[j][i] = (*o2_mRNA_mRNA)(o2_ind(i,j,dim));
  for(size_t i=0; i<dim; ++i)
    correlations[i][i] = (*o2_mRNA_mRNA)(o2_ind(i,i,dim)) + mRNA_expectations(i);

  return correlations;
}

const Compartment* Analytic_engine::set_o2_mRNA_prot_soma() {
  Soma& soma = *(p_neuron->p_soma);
  size_t sz_mRNA = 1+p_neuron->p_dend_segments.size();
  size_t sz_prot = sz_mRNA + p_neuron->p_synapses.size();
  
  o2_mRNA_prot_RHS = new arma::vec(sz_mRNA*sz_prot);
  o2_mRNA_prot_mat = new arma::mat(sz_mRNA*sz_prot, sz_mRNA*sz_prot);

  (*o2_mRNA_prot_RHS)(0) -= soma.translation_rate*soma.n_mRNA_expectation; // 0==o2_ind_asym(0,0,sz_prot)

  for(size_t i=0; i<sz_mRNA; ++i) {
    (*o2_mRNA_prot_RHS)(i*sz_prot) -= soma.translation_rate*(*o2_mRNA_mRNA)(i); //i==o2_ind(0,i,sz_mRNA)
    (*o2_mRNA_prot_mat)(i*sz_prot, i*sz_prot) -= soma.protein_decay_rate; // i*sz_prot==o2_ind_asym(i,0,sz_mRNA)
  }
  for(size_t i=0; i<sz_prot; ++i) {
    (*o2_mRNA_prot_RHS)(i) -= soma.transcription_rate*(*o2_gene_prot)(i);
    (*o2_mRNA_prot_mat)(i, i) -= soma.mRNA_decay_rate; // i==o2_ind_asym(0,i,sz_prot)
  }
  return &soma;
}

void Analytic_engine::set_o2_mRNA_prot_matrix(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty())
    return;

  size_t sz_mRNA = 1 + p_neuron->p_dend_segments.size(),
    sz_prot = sz_mRNA + p_neuron->p_synapses.size();
  
  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;

    size_t parent_mRNA_ind = p_junc -> p_from -> mRNA_ind-1,
      desc_mRNA_ind = p_junc -> p_to -> mRNA_ind-1,
      parent_prot_ind = p_junc -> p_from -> id,
      desc_prot_ind = p_junc -> p_to -> id;

    if(p_junc->type() != DEN_SYN) {

      (*o2_mRNA_prot_RHS)(o2_ind_asym(desc_mRNA_ind, desc_prot_ind, sz_prot)) -= (p_junc->p_to->translation_rate)*(p_junc->p_to->n_mRNA_expectation);
    
      for(size_t i=0; i<sz_mRNA; ++i) {
        (*o2_mRNA_prot_RHS)(o2_ind_asym(i, desc_prot_ind, sz_prot)) -= (p_junc->p_to->translation_rate)*(*o2_mRNA_mRNA)(o2_ind(desc_mRNA_ind, i, sz_mRNA));
        
        (*o2_mRNA_prot_mat)(o2_ind_asym(i,desc_prot_ind,sz_prot), o2_ind_asym(i,desc_prot_ind,sz_prot)) -= p_junc->p_to->protein_decay_rate;

        (*o2_mRNA_prot_mat)(o2_ind_asym(i, parent_prot_ind, sz_prot), o2_ind_asym(i, parent_prot_ind, sz_prot)) -= p_junc->fwd_prot_hop_rate;
        (*o2_mRNA_prot_mat)(o2_ind_asym(i, desc_prot_ind, sz_prot), o2_ind_asym(i, parent_prot_ind, sz_prot)) += p_junc->fwd_prot_hop_rate;
        (*o2_mRNA_prot_mat)(o2_ind_asym(i, desc_prot_ind, sz_prot), o2_ind_asym(i, desc_prot_ind, sz_prot)) -= p_junc->bkwd_prot_hop_rate;
        (*o2_mRNA_prot_mat)(o2_ind_asym(i, parent_prot_ind, sz_prot), o2_ind_asym(i, desc_prot_ind, sz_prot)) += p_junc->bkwd_prot_hop_rate;
      }
      for(size_t i=0; i<sz_prot; ++i) {
        (*o2_mRNA_prot_mat)(o2_ind_asym(desc_mRNA_ind,i,sz_prot), o2_ind_asym(desc_mRNA_ind,i,sz_prot)) -= p_junc->p_to->mRNA_decay_rate;

        (*o2_mRNA_prot_mat)(o2_ind_asym(parent_mRNA_ind, i, sz_prot), o2_ind_asym(parent_mRNA_ind, i, sz_prot)) -= p_junc->fwd_mRNA_hop_rate;
        (*o2_mRNA_prot_mat)(o2_ind_asym(desc_mRNA_ind, i, sz_prot), o2_ind_asym(parent_mRNA_ind, i, sz_prot)) += p_junc->fwd_mRNA_hop_rate;
        (*o2_mRNA_prot_mat)(o2_ind_asym(desc_mRNA_ind, i, sz_prot), o2_ind_asym(desc_mRNA_ind, i, sz_prot)) -= p_junc->bkwd_mRNA_hop_rate;
        (*o2_mRNA_prot_mat)(o2_ind_asym(parent_mRNA_ind, i, sz_prot), o2_ind_asym(desc_mRNA_ind, i, sz_prot)) += p_junc->bkwd_mRNA_hop_rate;
      }
    }
    else {
      for(size_t i=0; i<sz_mRNA; ++i) {
        (*o2_mRNA_prot_mat)(o2_ind_asym(i,desc_prot_ind,sz_prot), o2_ind_asym(i,desc_prot_ind,sz_prot)) -= p_junc->p_to->protein_decay_rate;

        (*o2_mRNA_prot_mat)(o2_ind_asym(i, parent_prot_ind, sz_prot), o2_ind_asym(i, parent_prot_ind, sz_prot)) -= p_junc->fwd_prot_hop_rate;
        (*o2_mRNA_prot_mat)(o2_ind_asym(i, desc_prot_ind, sz_prot), o2_ind_asym(i, parent_prot_ind, sz_prot)) += p_junc->fwd_prot_hop_rate;
        (*o2_mRNA_prot_mat)(o2_ind_asym(i, desc_prot_ind, sz_prot), o2_ind_asym(i, desc_prot_ind, sz_prot)) -= p_junc->bkwd_prot_hop_rate;
        (*o2_mRNA_prot_mat)(o2_ind_asym(i, parent_prot_ind, sz_prot), o2_ind_asym(i, desc_prot_ind, sz_prot)) += p_junc->bkwd_prot_hop_rate;
      }
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_o2_mRNA_prot_matrix(*((*it_p_junc)->p_to));
}

Analytic_engine& Analytic_engine::mRNA_protein_stationary_covariances() {
  set_o2_mRNA_prot_matrix(*set_o2_mRNA_prot_soma());
  
  o2_mRNA_prot = new arma::vec((*o2_mRNA_prot_mat).i()*(*o2_mRNA_prot_RHS));

  delete o2_mRNA_prot_mat; o2_mRNA_prot_mat=nullptr;
  delete o2_mRNA_prot_RHS; o2_mRNA_prot_RHS=nullptr;

  return *this;
}

std::vector<std::vector<double>> Analytic_engine::stationary_mRNA_protein_covariances() {

  mRNA_protein_stationary_covariances();

  size_t mRNA_dim = p_neuron->p_dend_segments.size() + 1,
    prot_dim = mRNA_dim + p_neuron->p_synapses.size();
  std::vector<std::vector<double>> correlations(mRNA_dim, std::vector<double>(prot_dim));
  for(size_t i=0; i<mRNA_dim; ++i)
    for(size_t j=0; j<prot_dim; ++j)
      correlations[i][j] = (*o2_mRNA_prot)(o2_ind_asym(i,j,prot_dim));
  
  return correlations;
}

const Compartment* Analytic_engine::set_o2_prot_prot_soma() {
  Soma& soma = *(p_neuron->p_soma);
  size_t sz = 1 + p_neuron->p_dend_segments.size() + p_neuron->p_synapses.size();
  
  o2_prot_prot_RHS = new arma::vec(sz*(sz+1)/2);
  o2_prot_prot_mat = new arma::mat(sz*(sz+1)/2, sz*(sz+1)/2);

  for(size_t i=0; i<sz; ++i) {
    (*o2_prot_prot_RHS)(i) -= soma.translation_rate*(*o2_mRNA_prot)(i);
    (*o2_prot_prot_mat)(i,i) -= soma.protein_decay_rate;
  }

  return &soma;
}

void Analytic_engine::set_o2_prot_prot_matrix(const Compartment& parent) {
  if(parent.it_p_out_junctions.empty())
    return;

  size_t sz = 1 + p_neuron->p_dend_segments.size() + p_neuron->p_synapses.size();

  for(auto& it_p_junc : parent.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;

    size_t parent_ind = p_junc -> p_from -> id;
    size_t desc_ind = p_junc -> p_to -> id;
    
    if(p_junc->type() != DEN_SYN) {
      size_t desc_mRNA_ind = p_junc -> p_to -> mRNA_ind-1;
    
      for(size_t i=0; i<sz; ++i) {
        (*o2_prot_prot_RHS)(o2_ind(desc_ind, i, sz)) -= p_junc->p_to->translation_rate*(*o2_mRNA_prot)(o2_ind_asym(desc_mRNA_ind, i, sz));
        
        (*o2_prot_prot_mat)(o2_ind(desc_ind, i, sz), o2_ind(desc_ind, i, sz)) -= p_junc->p_to->protein_decay_rate;

        (*o2_prot_prot_mat)(o2_ind(parent_ind, i, sz), o2_ind(parent_ind, i, sz)) -= p_junc->fwd_prot_hop_rate;
        (*o2_prot_prot_mat)(o2_ind(desc_ind, i, sz), o2_ind(parent_ind, i, sz)) += p_junc->fwd_prot_hop_rate;
        (*o2_prot_prot_mat)(o2_ind(desc_ind, i, sz), o2_ind(desc_ind, i, sz)) -= p_junc->bkwd_prot_hop_rate;
        (*o2_prot_prot_mat)(o2_ind(parent_ind, i, sz), o2_ind(desc_ind, i, sz)) += p_junc->bkwd_prot_hop_rate;
      }
    }
    else {
      for(size_t i=0; i<sz; ++i) {        
        (*o2_prot_prot_mat)(o2_ind(desc_ind, i, sz), o2_ind(desc_ind, i, sz)) -= p_junc->p_to->protein_decay_rate;

        (*o2_prot_prot_mat)(o2_ind(parent_ind, i, sz), o2_ind(parent_ind, i, sz)) -= p_junc->fwd_prot_hop_rate;
        (*o2_prot_prot_mat)(o2_ind(desc_ind, i, sz), o2_ind(parent_ind, i, sz)) += p_junc->fwd_prot_hop_rate;
        (*o2_prot_prot_mat)(o2_ind(desc_ind, i, sz), o2_ind(desc_ind, i, sz)) -= p_junc->bkwd_prot_hop_rate;
        (*o2_prot_prot_mat)(o2_ind(parent_ind, i, sz), o2_ind(desc_ind, i, sz)) += p_junc->bkwd_prot_hop_rate;
      }      
    }
  }

  for(auto& it_p_junc : parent.it_p_out_junctions)
    set_o2_prot_prot_matrix(*((*it_p_junc)->p_to));
}

Analytic_engine& Analytic_engine::protein_protein_stationary_covariances() {
  if(!(o2_gene_mRNA && o2_gene_prot && o2_mRNA_mRNA && o2_mRNA_prot)) {
    std::cerr << "-----------------------------------------------------------------------\n"
              << "*** ERROR in Analytic_engine::stationary_protein_protein_covariances():\n"
              << "* Prior to calling stationary_protein_protein_covariances() call\n"
              << "* stationary_mRNA_expectations()\n"
              << "* stationary_protein_expectations()\n"
              << "* gene_mRNA_stationary_covariances()\n"
              << "* mRNA_mRNA_stationary_covariances()\n"
              << "* gene_protein_stationary_covariances()\n"
              << "* mRNA_protein_stationary_covariances()\n"
              << "This is needed because protein-protein correleations depend on mRNA-protein correlations and active genes, mRNA and protein expectations.\n"
              << "-----------------------------------------------------------------------\n";
    exit(1);
  }
  
  set_o2_prot_prot_matrix(*set_o2_prot_prot_soma());
  
  o2_prot_prot = new arma::vec((*o2_prot_prot_mat).i()*(*o2_prot_prot_RHS));

  delete o2_prot_prot_mat; o2_prot_prot_mat=nullptr;
  delete o2_prot_prot_RHS; o2_prot_prot_RHS=nullptr;

  return *this;
}

std::vector<std::vector<double>> Analytic_engine::stationary_protein_protein_covariances() {

  protein_protein_stationary_covariances();

  size_t dim = p_neuron->p_dend_segments.size() + p_neuron->p_synapses.size() + 1;
  std::vector<std::vector<double>> correlations(dim, std::vector<double>(dim));
  for(size_t i=0; i<dim; ++i)
    for(size_t j=0; j<i; ++j)
      correlations[i][j] = correlations[j][i] = (*o2_prot_prot)(o2_ind(i,j,dim));
  for(size_t i=0; i<dim; ++i)
    correlations[i][i] = (*o2_prot_prot)(o2_ind(i,i,dim)) + protein_expectations(i);
  
  return correlations;
}

Analytic_engine& Analytic_engine::protein_protein_stationary_covariances(std::ofstream &ofs) {
  std::cout << "Computing protein_protein_stationary_covariances...\n";
  set_o2_prot_prot_matrix(*set_o2_prot_prot_soma());
  
  o2_prot_prot = new arma::vec((*o2_prot_prot_mat).i()*(*o2_prot_prot_RHS));

  std::cerr << "protein_protein_covariances:\n" << *o2_prot_prot << std::endl;

  delete o2_prot_prot_mat; o2_prot_prot_mat=nullptr;
  delete o2_prot_prot_RHS; o2_prot_prot_RHS=nullptr;

  //////////// COMPUTING VARIANCES and PCCs //////////
  std::cout << "Protein means and standard deviations:\n";
  size_t sz = 1 + p_neuron->p_dend_segments.size() + p_neuron->p_synapses.size();
  std::vector<double> rmss(sz);
  for(size_t i=0; i<sz; ++i) {
    rmss[i] = sqrt((*o2_prot_prot)[o2_ind(i, i, sz)] - protein_expectations[i]*(protein_expectations[i]-1));
    if(rmss[i]>=0) {
      std::cout << o1_prot_names[i] + ": " << protein_expectations(i) << ", " << rmss[i] << std::endl;
    }
    else
      std::cout << o1_prot_names[i] + ": " << protein_expectations(i) << ", " << rmss[i] << " NEGATIVE!\n";
  }

  // Writing to file
  for(unsigned int i=0; i<sz; ++i)
    ofs << o1_prot_names[i] << ',';
  for(unsigned int i=0; i<sz; ++i) {
    ofs << std::endl << o1_prot_names[i] << ',';
    for(unsigned int j=0; j<sz; ++j)
      if(i != j)
        ofs << ((*o2_prot_prot)[o2_ind(i, j, sz)] - protein_expectations(i)*protein_expectations(j))/(rmss[i]*rmss[j]) << ',';
      else
        ofs << "1,";
  }

  return *this;
}

Analytic_engine& Analytic_engine::stationary_covariances(bool write_covariance_matrix) {
  initialise_o2();
  set_o2_matrix();
  set_o2_RHS();
  // std::cerr << "(*p_o2_RHS)(j) = " << (*p_o2_RHS) << std::endl;
  auto& covariances = *p_covariances;
  auto& o2_mat = *p_o2_mat;
  auto& o2_RHS = *p_o2_RHS;
  std::cerr << "Computing covariances...\n";
  covariances = o2_mat.i()*o2_RHS;

  // std::vector<double> rmss(o1_dim);
  // for(size_t i=0; i<o1_dim; ++i) {  
  //   rmss[i] = sqrt(covariances(o2_ind(i,i)) - expectations(i)*(expectations(i)-1));
  //   std::cerr << o1_var_names[i] + ": " << expectations(i) << ", " << rmss[i] << std::endl;
  // }

  if(p_cov_mat && write_covariance_matrix) // Setting covariance matrix for direct ODE solvers
    for(size_t i=0; i<o1_dim; ++i)
      for(size_t j=0; j<i; ++j)
        (*p_cov_mat)(j,i) = (*p_cov_mat)(i,j) = covariances(o2_ind(i,j));
  for(size_t i=0; i<o1_dim; ++i)
    (*p_cov_mat)(i,i) = (*p_cov_mat)(i,i) = covariances(o2_ind(i,i)) + expectations(i);
  
  return *this;
}

Analytic_engine& Analytic_engine::sem_stationary_covariances() {
  sem_initialise_o2();
  sem_set_o2_matrix();
  sem_set_o2_RHS();
  
  auto& covariances = *p_covariances;
  auto& o2_mat = *p_o2_mat;
  auto& o2_RHS = *p_o2_RHS;
  std::cerr << "Computing covariances...\n";
  covariances = o2_mat.i()*o2_RHS;

  std::vector<double> rmss(o1_dim);
  for(size_t i=0; i<o1_dim; ++i) {
    rmss[i] = sqrt(covariances(sem_o2_ind(i,i)) - expectations(i)*(expectations(i)-1));
    std::cerr << o1_var_names[i] + ": " << expectations(i) << ", " << rmss[i]<< ", " << rmss[i]/expectations(i) << std::endl;
  }

  // Writing covariance matrix
  for(size_t i=0; i<o1_dim; ++i)
    for(size_t j=i+1; j<o1_dim; ++j)
      (*p_cov_mat)(j,i) = (*p_cov_mat)(i,j) = covariances(sem_o2_ind(i,j));
  for(size_t i=0; i<o1_dim; ++i)
    (*p_cov_mat)(i,i) = covariances(sem_o2_ind(i,i)) + expectations(i);

  std::cerr << "--,";
  for(unsigned int i=0; i<o1_dim; ++i)
    std::cerr << o1_var_names[i] << ',';
  for(unsigned int i=0; i<o1_dim; ++i) {
    std::cerr << std::endl << o1_var_names[i] << ',';
    for(unsigned int j=0; j<o1_dim; ++j)
      std::cerr << covariances(sem_o2_ind(i,j)) << ',';
  }
  std::cerr << std::endl;
  
  return *this;
}

Analytic_engine& Analytic_engine::sem_stationary_pearson_correlations() {

  sem_stationary_covariances();
  auto& covariances = *p_covariances;

  std::vector<double> rmss(o1_dim);
  for(size_t i=0; i<o1_dim; ++i) {  
    rmss[i] = sqrt(covariances(sem_o2_ind(i,i)) - expectations(i)*(expectations(i)-1));
    std::cerr << o1_var_names[i] + ": " << expectations(i) << ", " << rmss[i]<< ", " << rmss[i]/expectations(i) << std::endl;
  }

  std::cerr << "\n----------------------------------------------------\n";

  std::cerr << "--,";
  for(unsigned int i=0; i<o1_dim; ++i)
    std::cerr << o1_var_names[i] << ',';
  for(unsigned int i=0; i<o1_dim; ++i) {
    std::cerr << std::endl << o1_var_names[i] << ',';
    for(unsigned int j=0; j<o1_dim; ++j)
      if(i != j)
        std::cerr << (covariances(sem_o2_ind(i,j)) - expectations(i)*expectations(j))/(rmss[i]*rmss[j]) << ',';
      else
        std::cerr << "1,";
  }
  std::cerr << std::endl;
  
  return *this;
}

Analytic_engine& Analytic_engine::sem_stationary_time_correlations(const std::list<double>& times) {
  std::cout << "sem_stationary_time_correlations method\n";
  std::cerr << "Setting o1_matrix...\n";
  sem_set_o1_matrix(*sem_set_o1_soma());
  auto& o1_mat = *p_o1_mat;
  o1_mat = -o1_mat;
    
  std::cerr << "Computing eigen decomposition...\n";
  arma::cx_vec eigval_c;
  arma::vec eigval(o1_dim);
  arma::cx_mat eigvec_c;
  arma::mat tm(o1_dim, o1_dim); // Transition matrix
  arma::eig_gen(eigval_c, eigvec_c, o1_mat);
  
  for(size_t i=0; i<o1_dim; ++i) {
    eigval(i) = eigval_c(i).real();
    for(size_t j=0; j<o1_dim; ++j)
      tm(i,j) = eigvec_c(i,j).real();
  }

  std::cerr << "Inverting transition matrix...\n";
  arma::mat inv_tm = tm.i(); // Transition matrix
  
  std::cerr << "Eigenvalues:\n"
            << eigval << std::endl;
  std::cerr << "Eigenvectors:\n"
            << tm << std::endl;

  for(auto& o1_var_name : o1_var_names)
    std::cout << o1_var_name << ',';
  std::cout << std::endl;
  
  // Computing stationary part
  for(size_t i=0; i<o1_dim; ++i) {
    double sum = 0;
    for(size_t k=0; k<o1_dim; ++k)
      sum += 1/eigval[k]*tm(i,k)*inv_tm(k,0);
    expectations(i) = p_neuron->p_soma->gene_activation_rate * p_neuron->p_soma->number_of_gene_copies * sum;
  }

  std::cout << "Stationary expectations from nonstationary covariances algorithm:\n";
  for(size_t i=0; i<o1_dim; ++i)
    std::cout << o1_var_names[i] << ": " << expectations(i) << std::endl;

  // Computing covariances
  sem_stationary_covariances();

  std::vector<double> rmss(o1_dim);
  for(size_t i=0; i<o1_dim; ++i) {
    rmss[i] = sqrt((*p_covariances)(sem_o2_ind(i,i)) - expectations(i)*(expectations(i)-1));
    std::cerr << o1_var_names[i] + ": " << expectations(i) << ", " << rmss[i]<< ", " << rmss[i]/expectations(i) << std::endl;
  }

  for(auto t : times) {
    // Computing matrix exponent
    arma::mat matrix_exp(o1_dim, o1_dim);
    for(size_t i=0; i<o1_dim; ++i)
      for(size_t j=0; j<o1_dim; ++j)
        for(size_t k=0; k<o1_dim; ++k)
          matrix_exp(i,j) += exp(-eigval(k)*t) * tm(i,k) * inv_tm(k,j);

    arma::mat correlator = matrix_exp * (*p_cov_mat);
    for(size_t i=0; i<o1_dim; ++i)
      for(size_t j=0; j<o1_dim; ++j) {
        arma::vec v = matrix_exp*expectations;
        correlator(i,j) += - v(i)*expectations(j);
        correlator(i,j) /= (rmss[i]*rmss[j]); // Pearson correlation coefficient
        // correlator(i,j) += (expectations(i) - v(i))*expectations(j);
      }
    
    // std::cout << "t= " << t << std::endl
    //           << "correlator:\n" << correlator << std::endl;

    // std::cerr << t << ',' << correlator(6,6) << ',' << correlator(5,6) << ',' << correlator(4,6) << std::endl;
    // std::cerr << t << ',' << correlator(14,14) << ',' << correlator(13,14) << ',' << correlator(11,14)  << ',' << correlator(10,14) << ',' << correlator(8,14) << ',' << correlator(7,14) //Syn
    //           << ',' << correlator(12,14) << ',' << correlator(9,14)  << ',' << correlator(6,14) << ',' << correlator(5,14) << std::endl; //Dend segments for 3 ds

    std::cerr << t << ',' << correlator(16,16) << ',' << correlator(15,16) << ',' << correlator(14,16)  << ',' << correlator(13,16) << ',' << correlator(12,16) << ',' << correlator(11,16) << ',' << correlator(10,16) << ',' << correlator(9,16)  << ',' << correlator(8,16) << ',' << correlator(7,16) << ',' << correlator(6,16) << std::endl; //Dend segments for 4 ds
  }
                
  return *this;
}

double Analytic_engine::active_genes_expectation() {
  return p_neuron->p_soma->n_active_genes_expectation;
}
double Analytic_engine::mRNA_expectation(const Compartment& comp) {
  if(comp.type() != SPINE)
    return mRNA_expectations(comp.mRNA_ind-1);
  else
    return 0;
}
double Analytic_engine::protein_expectation(const Compartment& comp) {
  return protein_expectations(comp.id);
}
double Analytic_engine::gene_mRNA_correlation(const Compartment& comp) {
  if(!o2_gene_mRNA) {
    std::cerr << "--------------------------------------------------------------\n"
              << "ERROR in gene_mRNA_correlation: Correlations are not computed.\n"
              << "--------------------------------------------------------------\n";
  }
  if(comp.type() != SPINE)
    return (*o2_gene_mRNA)(comp.mRNA_ind-1);
  else
    return 0;
}
double Analytic_engine::mRNA_mRNA_correlation(const Compartment& comp1, const Compartment& comp2) {
  if(!o2_mRNA_mRNA) {
    std::cerr << "--------------------------------------------------------------\n"
              << "ERROR in mRNA_mRNA_correlation: Correlations are not computed.\n"
              << "--------------------------------------------------------------\n";
  }
  size_t mRNA_dim = 1+p_neuron->p_dend_segments.size();
  if(comp1.type()!=SPINE && comp1.type()!=SPINE) {
    if(comp1.mRNA_ind != comp2.mRNA_ind)
      return (*o2_mRNA_mRNA)(o2_ind(comp1.mRNA_ind-1,comp2.mRNA_ind-1,mRNA_dim));
    else
      return (*o2_mRNA_mRNA)(o2_ind(comp1.mRNA_ind-1,comp2.mRNA_ind-1,mRNA_dim)) + mRNA_expectations(comp1.mRNA_ind-1);
  }
  else
    return 0;
}
double Analytic_engine::gene_protein_correlation(const Compartment& comp) {
  if(!o2_gene_prot) {
    std::cerr << "-----------------------------------------------------------------\n"
              << "ERROR in gene_protein_correlation: Correlations are not computed.\n"
              << "-----------------------------------------------------------------\n";
    exit(1);
  }
  return (*o2_gene_prot)(comp.id);
}
double Analytic_engine::mRNA_protein_correlation(const Compartment& comp1, const Compartment& comp2) {
  if(!o2_mRNA_prot) {
    std::cerr << "-----------------------------------------------------------------\n"
              << "ERROR in mRNA_protein_correlation: Correlations are not computed.\n"
              << "-----------------------------------------------------------------\n";
    exit(1);
  }
  size_t prot_dim = 1+p_neuron->p_dend_segments.size()+p_neuron->p_synapses.size();
  if(comp1.type() != SPINE)
    return (*o2_mRNA_prot)(o2_ind_asym(comp1.mRNA_ind-1,comp2.id, prot_dim));
  else
    return 0;
}
double Analytic_engine::protein_protein_correlation(const Compartment& comp1, const Compartment& comp2) {
  if(!o2_prot_prot) {
    std::cerr << "--------------------------------------------------------------------\n"
              << "ERROR in protein_protein_correlation: Correlations are not computed.\n"
              << "--------------------------------------------------------------------\n";
    exit(1);
  }
  size_t prot_dim = 1+p_neuron->p_dend_segments.size()+p_neuron->p_synapses.size();
  if(comp1.id != comp2.id)
      return (*o2_prot_prot)(o2_ind(comp1.id, comp2.id, prot_dim));
    else
      return (*o2_prot_prot)(o2_ind(comp1.id, comp2.id, prot_dim)) + protein_expectations(comp1.id);
}

Analytic_engine& Analytic_engine::clear_o1_mat_and_RHS() {
  if(p_o1_mat) {
    delete p_o1_mat;
    p_o1_mat = nullptr;
  }
  if(p_o1_RHS) {
    delete p_o1_RHS;
    p_o1_RHS = nullptr;
  }
  return *this;
}

Analytic_engine& Analytic_engine::clear_As_and_bs() {
  if(p_Ap) {
    delete p_Ap;
    p_Ap = nullptr;
  }
  if(p_Am) {
    delete p_Am;
    p_Am = nullptr;
  }
  if(p_b) {
    delete p_b;
    p_b = nullptr;
  }
  return *this;
}

Analytic_engine& Analytic_engine::clear_mRNA_As() {
  if(p_mRNA_Ap) {
    delete p_mRNA_Ap;
    p_mRNA_Ap = nullptr;
  }
  if(p_mRNA_Am) {
    delete p_mRNA_Am;
    p_mRNA_Am = nullptr;
  }
  return *this;
}

Analytic_engine& Analytic_engine::clear_prot_As() {
  if(p_prot_Ap) {
    delete p_prot_Ap;
    p_prot_Ap = nullptr;
  }
  if(p_prot_Am) {
    delete p_prot_Am;
    p_prot_Am = nullptr;
  }
  if(p_PM) {
    delete p_PM;
    p_PM = nullptr;
  }
  return *this;
}

Analytic_engine& Analytic_engine::clear_hopping_rate_matrix() {
  if(p_H) {
    delete p_H;
    p_H = nullptr;
  }
  return *this;
}

Analytic_engine& Analytic_engine::clear_mRNA_hopping_rate_matrix() {
  if(p_mRNA_H) {
    delete p_mRNA_H;
    p_mRNA_H = nullptr;
  }
  return *this;
}

Analytic_engine& Analytic_engine::clear_prot_hopping_rate_matrix() {
  if(p_prot_H) {
    delete p_prot_H;
    p_prot_H = nullptr;
  }
  return *this;
}

Analytic_engine& Analytic_engine::clear_o1() {
  clear_o1_mat_and_RHS();
  clear_As_and_bs();
  clear_mRNA_As();
  clear_prot_As();
  clear_hopping_rate_matrix();
  clear_mRNA_hopping_rate_matrix();
  clear_prot_hopping_rate_matrix();
  // Other stuff like o1_var_names that is currently in stack
  return *this;
}

Analytic_engine& Analytic_engine::clear_o2_mat_and_RHS() {
  if(p_o2_mat) {delete p_o2_mat; p_o2_mat = nullptr;}
  if(p_o2_RHS) {delete p_o2_RHS; p_o2_RHS = nullptr;}

  return *this;
}

Analytic_engine& Analytic_engine::clear_o2() {
  clear_o2_mat_and_RHS();
  if(p_o2_var_names) { delete p_o2_var_names; p_o2_var_names = nullptr;}
  if(p_covariances) {delete p_covariances; p_covariances=nullptr;}
  if(p_cov_mat) {delete p_cov_mat; p_cov_mat=nullptr;}
  if(p_mRNA_mRNA_cov_mat) {delete p_mRNA_mRNA_cov_mat; p_mRNA_mRNA_cov_mat=nullptr;}
  if(p_prot_prot_cov_mat) {delete p_prot_prot_cov_mat; p_prot_prot_cov_mat=nullptr;}
  if(p_PM) {delete p_PM; p_PM=nullptr;}
  if(o2_gene_mRNA) {delete o2_gene_mRNA; o2_gene_mRNA=nullptr;}
  if(o2_gene_mRNA_RHS) {delete o2_gene_mRNA_RHS; o2_gene_mRNA_RHS=nullptr;}
  if(o2_gene_mRNA_mat) {delete o2_gene_mRNA_mat; o2_gene_mRNA_mat=nullptr;}
  if(o2_gene_prot) {delete o2_gene_prot; o2_gene_prot=nullptr;}
  if(o2_gene_prot_RHS) {delete o2_gene_prot_RHS; o2_gene_prot_RHS=nullptr;}
  if(o2_gene_prot_mat) {delete o2_gene_prot_mat; o2_gene_prot_mat=nullptr;}
  if(o2_mRNA_mRNA) {delete o2_mRNA_mRNA; o2_mRNA_mRNA=nullptr;}
  if(o2_mRNA_mRNA_RHS) {delete o2_mRNA_mRNA_RHS; o2_mRNA_mRNA_RHS=nullptr;}
  if(o2_mRNA_mRNA_mat) {delete o2_mRNA_mRNA_mat; o2_mRNA_mRNA_mat=nullptr;}
  if(o2_mRNA_prot) {delete o2_mRNA_prot; o2_mRNA_prot=nullptr;}
  if(o2_mRNA_prot_RHS) {delete o2_mRNA_prot_RHS; o2_mRNA_prot_RHS=nullptr;}
  if(o2_mRNA_prot_mat) {delete o2_mRNA_prot_mat; o2_mRNA_prot_mat=nullptr;}
  if(o2_nonstationary_RHS_mat) {delete o2_nonstationary_RHS_mat; o2_nonstationary_RHS_mat=nullptr;}

  return *this;
}

Analytic_engine& Analytic_engine::clear_all() {
  clear_o1();
  clear_o2();
  
  return *this;
}

Analytic_engine::~Analytic_engine() {
  clear_all();
}
