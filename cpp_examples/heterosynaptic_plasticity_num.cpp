#include <vector>
#include "../Neuron.hpp"
#include "../compartments/Soma.hpp"
#include "../compartments/Dendritic_segment.hpp"
#include "../compartments/Spine.hpp"
#include "../engines/analytic/Analytic_engine.hpp"

using namespace std;

#define syn_dec_rate 0 //1.2e-5*3600
#define transcription_rate_multiplier 2
#define binding_rate_multiplier 10
#define SYN syn_1_1

int main() {

  Soma soma("soma" /*,Parameters of the soma*/);

  ///// Branching neuron
  std::list<Spine*> p_spines;
  Dendritic_segment ds(soma, "d_1");
  Spine syn_1_1(ds, "s_1_1", .6, 6, syn_dec_rate);
  p_spines.push_back(&syn_1_1);
  Spine syn_1_2(ds, "s_1_2", .6, 6, syn_dec_rate);
  p_spines.push_back(&syn_1_2);
  Dendritic_segment ds_1(ds, "d_1_1");
  Spine syn_11_1(ds_1, "s_1_1-1", .6, 6, syn_dec_rate);
  p_spines.push_back(&syn_11_1);
  Spine syn_11_2(ds_1, "s_1_1-2", .6, 6, syn_dec_rate);
  p_spines.push_back(&syn_11_2);
  Dendritic_segment ds_2(ds, "d_1_2");
  Spine syn_12_1(ds_2, "s_1_2-1", .6, 6, syn_dec_rate);
  p_spines.push_back(&syn_12_1);
  Spine syn_12_2(ds_2, "s_1_2-2", .6, 6, syn_dec_rate);
  p_spines.push_back(&syn_12_2);



  Neuron neuron(soma, "Test_neuron");
  
  cout << neuron << endl;

  // cout << "----------------- ANALYTIC ENGINE -----------------\n";
  Analytic_engine ae(neuron);

#define dim 15
  
  arma::mat covariances(dim,dim);
  arma::vec expectations(dim), variances(dim);

  ofstream ofs_expectations("2_genes_expectations_"+ SYN.get_name() + "_transcription_rate_multiplier_" + std::to_string(transcription_rate_multiplier) + "_bind_r_mult_" + std::to_string(binding_rate_multiplier)),
    // ofs_covariances("covariances_"+ SYN.get_name() + "_transcription_rate_multiplier_" + std::to_string(transcription_rate_multiplier) + "_bind_r_mult_" + std::to_string(binding_rate_multiplier)),
    ofs_variances("2_genes_stds_"+ SYN.get_name() + "_transcription_rate_multiplier_" + std::to_string(transcription_rate_multiplier) + "_bind_r_mult_" + std::to_string(binding_rate_multiplier));
    //ofs_correlations("correlations_"+ SYN.get_name() + "_transcription_rate_multiplier_" + std::to_string(transcription_rate_multiplier) + "_bind_r_mult_" + std::to_string(binding_rate_multiplier));

#define t_start 0
#define t1 10
#define t2 500
  double dt = .001;

  ae.stationary_expectations().stationary_covariances(true); // Initialising at stationarity
  
  std::cerr << "------------------- Loop_1 -----------------------\n";
  for(double t=t_start; t<t1; t+=dt) {
    expectations = ae.get_expectations();
    ofs_expectations << t;
    for(size_t i=0; i<dim; ++i)
      ofs_expectations  << ',' << expectations(i);
    ofs_expectations << endl;
    
    if(t==0)
      ae.nonstationary_expectations_direct_ODE_solver_step(dt, true);
    else
      ae.nonstationary_expectations_direct_ODE_solver_step(dt);

    // ofs_covariances << "t=" << t << endl
    //                 << "covariances:\n" << (covariances = ae.get_covariances()) << endl;

    covariances = ae.get_covariances();
    
    // ofs_correlations << "t=" << t << std::endl;
    // for(size_t i=0; i<dim; ++i) {
    //   for(size_t j=0; j<dim; ++j)
    //     ofs_correlations << covariances(i,j) - expectations(i)*expectations(j) << ',';
    //   ofs_correlations << std::endl;
    // }        

    ofs_variances << t;
    for(size_t i=0; i<dim; ++i)
      ofs_variances  << ',' << sqrt(covariances(i,i) - expectations(i)*expectations(i));
    ofs_variances << endl;

    ae.nonstationary_covariances_direct_ODE_solver_step(dt);
  }

  std::cerr << "------------------- Loop_2 -----------------------\n";

  SYN.set_protein_binding_rate(.6*binding_rate_multiplier);
  soma.set_transcription_rate(3.*200/10000*0.001*3600*transcription_rate_multiplier);

  neuron.refresh();
  std::cerr << "neuron:\n" << neuron << std::endl;
  
  for(double t=t1; t<t2; t+=dt) {
    ofs_expectations << t;
    expectations = ae.get_expectations();
    for(size_t i=0; i<dim; ++i)
      ofs_expectations  << ',' << expectations(i);
    ofs_expectations << endl;
    
    if(t==t1)
      ae.nonstationary_expectations_direct_ODE_solver_step(dt, true);
    else
      ae.nonstationary_expectations_direct_ODE_solver_step(dt);

    // ofs_covariances << "t=" << t << endl
    //                 << "covariances:\n" << (covariances = ae.get_covariances()) << endl;

    covariances = ae.get_covariances();
      
    // ofs_correlations << "t=" << t << std::endl;
    // for(size_t i=0; i<dim; ++i) {
    //   for(size_t j=0; j<dim; ++j)
    //     ofs_correlations << covariances(i,j) - expectations(i)*expectations(j) << ',';
    //   ofs_correlations << std::endl;
    // }        

    ofs_variances << t;
    for(size_t i=0; i<dim; ++i)
      ofs_variances  << ',' << sqrt(covariances(i,i) - expectations(i)*expectations(i));
    ofs_variances << endl;

    ae.nonstationary_covariances_direct_ODE_solver_step(dt);
  }

  ofs_expectations.close();
  // ofs_covariances.close();
  ofs_variances.close();

  ae.stationary_expectations().stationary_covariances();
  
  return 0;
}
