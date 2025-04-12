#include <vector>
#include "../Neuron.hpp"
#include "../compartments/Soma.hpp"
#include "../compartments/Dendritic_segment.hpp"
#include "../compartments/Spine.hpp"
#include "../engines/stochastic/Gillespie_engine.hpp"
#include "../engines/analytic/Analytic_engine.hpp"

#define N_AVRG      10000
#define time_offset -4000
#define duration    7000
#define  dt        .01

void stationary_moments() {
  std::cerr << "Computing stationary moments...\n";

  ///// Branching neuron
  Soma soma("soma", 40);
  std::list<Spine*> p_spines;
  Dendritic_segment ds(soma, "d_1");
  Spine syn_1_1(ds, "s_1_1", .6, 6);
  p_spines.push_back(&syn_1_1);
  Spine syn_1_2(ds, "s_1_2", .6, 6);
  p_spines.push_back(&syn_1_2);
  Dendritic_segment ds_1(ds, "d_1_1");
  Spine syn_11_1(ds_1, "s_1_1-1", .6, 6);
  p_spines.push_back(&syn_11_1);
  Spine syn_11_2(ds_1, "s_1_1-2", .6, 6);
  p_spines.push_back(&syn_11_2);
  Dendritic_segment ds_2(ds, "d_1_2");
  Spine syn_12_1(ds_2, "s_1_2-1", .6, 6);
  p_spines.push_back(&syn_12_1);
  Spine syn_12_2(ds_2, "s_1_2-2", .6, 6);
  p_spines.push_back(&syn_12_2);

  Neuron neuron(soma, "Y_neuron");

  Analytic_engine ae(neuron);
  size_t dim = ae.dim_o1();
  arma::mat covariances(dim,dim);
  arma::vec expectations(dim), variances(dim);
    
  std::ofstream ofs_expectations("stationary_expectations_2_genes"),
    ofs_stds("stationary_stds_2_genes");

  ae.stationary_expectations().stationary_covariances(true); // Initialising at stationarity

  expectations = ae.get_expectations();
  for(size_t i=0; i<dim; ++i)
    ofs_expectations << expectations(i) << ',';
  ofs_expectations << std::endl;

  covariances = ae.get_covariances();

  for(size_t i=0; i<dim; ++i)
    ofs_stds << sqrt(covariances(i,i) - expectations(i)*expectations(i)) << ',';
  ofs_stds << std::endl;
    
  ofs_expectations.close();
  ofs_stds.close();
}

int main() {

  stationary_moments();
  
  PRNG rnd(1);
  
  std::string file_name = "../../data/gillespie/stationary_moments_for_cosyne_2_genes_more_more_more/SM_";

  std::list<double> times;
  for(double t=0; t<duration; t+=dt)
    times.push_back(t);  
  
  for(size_t run_ind=0; run_ind<N_AVRG; ++run_ind) {

    Soma soma("soma", 40);

    ///// Branching neuron
    std::list<Spine*> p_spines;
    Dendritic_segment ds(soma, "d_1");
    Spine syn_1_1(ds, "s_1_1", .6, 6);
    p_spines.push_back(&syn_1_1);
    Spine syn_1_2(ds, "s_1_2", .6, 6);
    p_spines.push_back(&syn_1_2);
    Dendritic_segment ds_1(ds, "d_1_1");
    Spine syn_11_1(ds_1, "s_1_1-1", .6, 6);
    p_spines.push_back(&syn_11_1);
    Spine syn_11_2(ds_1, "s_1_1-2", .6, 6);
    p_spines.push_back(&syn_11_2);
    Dendritic_segment ds_2(ds, "d_1_2");
    Spine syn_12_1(ds_2, "s_1_2-1", .6, 6);
    p_spines.push_back(&syn_12_1);
    Spine syn_12_2(ds_2, "s_1_2-2", .6, 6);
    p_spines.push_back(&syn_12_2);

    Neuron neuron(soma, "Y_neuron");
  

    std::ofstream ofs_gillespie(file_name + std::to_string(run_ind));

    std::cerr << "Writing Gillespie results to: " << file_name + std::to_string(run_ind) << '\n';
    
    Gillespie_engine(neuron, rnd).run_Gillespie(times, ofs_gillespie, time_offset);
    
    ofs_gillespie.close();
    
  }
  
  return 0;
}
