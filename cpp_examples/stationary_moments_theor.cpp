#include <vector>
#include "../Neuron.hpp"
#include "../compartments/Soma.hpp"
#include "../compartments/Dendritic_segment.hpp"
#include "../compartments/Spine.hpp"
#include "../engines/analytic/Analytic_engine.hpp"

#define N_AVRG      10000
#define time_offset -4000
#define duration    7000
#define  dt        .01

void stationary_moments() {
  std::cerr << "Computing stationary moments...\n";

  ///// Branching neuron
  Soma soma("soma", 40);
  std::list<Spine*> p_synapses;
  Dendritic_segment ds(soma, "d_1");
  Spine syn_1_1(ds, "s_1_1", .6, 6);
  p_synapses.push_back(&syn_1_1);
  Spine syn_1_2(ds, "s_1_2", .6, 6);
  p_synapses.push_back(&syn_1_2);
  Dendritic_segment ds_1(ds, "d_1_1");
  Spine syn_11_1(ds_1, "s_1_1-1", .6, 6);
  p_synapses.push_back(&syn_11_1);
  Spine syn_11_2(ds_1, "s_1_1-2", .6, 6);
  p_synapses.push_back(&syn_11_2);
  Dendritic_segment ds_2(ds, "d_1_2");
  Spine syn_12_1(ds_2, "s_1_2-1", .6, 6);
  p_synapses.push_back(&syn_12_1);
  Spine syn_12_2(ds_2, "s_1_2-2", .6, 6);
  p_synapses.push_back(&syn_12_2);

  Neuron neuron(soma, "Y_neuron");

  Analytic_engine ae(neuron);
  size_t dim = ae.dim_o1();
  arma::mat covariances(dim,dim);
  arma::vec expectations(dim), stds(dim);
    
  std::ofstream ofs_expectations("stationary_expectations_2_genes"),
    ofs_stds("stationary_stds_2_genes"),
    ofs_covariances("stationary_covariances_2_genes"),
    ofs_pearson_correlations("stationary_pearson_correlations_2_genes"),
    ofs_var_names("var_names_2_genes");

  ae.stationary_expectations().stationary_covariances(true); // Initialising at stationarity
  
  auto var_names = *ae.o1_variable_names();
  ofs_var_names << var_names[0];
  for(size_t i=1; i<dim; ++i)
    ofs_var_names << ',' << var_names[i];
  ofs_var_names.close();

  expectations = ae.get_expectations();
  for(size_t i=0; i<dim; ++i)
    ofs_expectations << expectations(i) << ',';
  ofs_expectations << std::endl;

  covariances = ae.get_covariances();
  for(size_t i=0; i<dim; ++i) {
    ofs_covariances << covariances(i,0);
    for(size_t j=1; j<dim; ++j) 
      ofs_covariances << ',' << covariances(i,j);
    ofs_covariances << std::endl;
  }

  for(size_t i=0; i<dim; ++i)
    ofs_stds << (stds(i) = sqrt(covariances(i,i) - expectations(i)*expectations(i))) << ',';
  ofs_stds << std::endl;


  for(size_t i=0; i<dim; ++i) {
    ofs_pearson_correlations << (covariances(i,0) - expectations(i)*expectations(0))/(stds(i)*stds(0));
    for(size_t j=1; j<dim; ++j) 
      ofs_pearson_correlations << ',' << (covariances(i,j) - expectations(i)*expectations(j))/(stds(i)*stds(j));
    ofs_pearson_correlations << std::endl;
  }
    
  ofs_expectations.close();
  ofs_stds.close();
  ofs_covariances.close();
  ofs_pearson_correlations.close();
}

int main() {

  stationary_moments();
   
  return 0;
}
