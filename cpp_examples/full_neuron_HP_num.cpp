#include <vector>
#include "../Neuron.hpp"
#include "../compartments/Soma.hpp"
#include "../compartments/Dendritic_segment.hpp"
#include "../compartments/Spine.hpp"
#include "../engines/analytic/Analytic_engine.hpp"

using namespace std;

int main() {

  Soma soma("soma" /*,Parameters of the soma*/);

  ///// Branching neuron
  Dendritic_segment ds(soma, "d_1");
  Spine syn_1_1(ds, "s_1_1", .6, 6, 1.2e-5*3600 * 10);
  Spine syn_1_2(ds, "s_1_2", .6, 6, 1.2e-5*3600 * 10);
  Dendritic_segment ds_1(ds, "d_1-1");
  Spine syn_11_1(ds_1, "s_1_1-1", .6, 6, 1.2e-5*3600 * 10);
  Spine syn_11_2(ds_1, "s_1_1-2", .6, 6, 1.2e-5*3600 * 10);
  Dendritic_segment ds_2(ds, "d_1-1");
  Spine syn_12_1(ds_2, "s_1_2-1", .6, 6, 1.2e-5*3600 * 10);
  Spine syn_12_2(ds_2, "s_1_2-2", .6, 6, 1.2e-5*3600 * 10);

  Neuron neuron(soma, "Test_neuron");
  
  cout << neuron << endl;

  // cout << "----------------- ANALYTIC ENGINE -----------------\n";
  Analytic_engine ae(neuron);

#define dim 15
  
  arma::mat covariances(dim,dim);
  arma::vec expectations(dim), variances(dim);

  ofstream ofs_expectations("expectations"),
    ofs_covariances("covariances"),
    ofs_variances("variances"),
    ofs_correlations("correlations");

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

    ofs_covariances << "t=" << t << endl
                    << "covariances:\n" << (covariances = ae.get_covariances()) << endl;
    
    ofs_correlations << "t=" << t << std::endl;
    for(size_t i=0; i<dim; ++i) {
      for(size_t j=0; j<dim; ++j)
        ofs_correlations << covariances(i,j) - expectations(i)*expectations(j) << ',';
      ofs_correlations << std::endl;
    }        

    ofs_variances << t;
    for(size_t i=0; i<dim; ++i)
      ofs_variances  << ',' << sqrt(covariances(i,i) - expectations(i)*expectations(i));
    ofs_variances << endl;

    ae.nonstationary_covariances_direct_ODE_solver_step(dt);
  }

  std::cerr << "------------------- Loop_2 -----------------------\n";

  syn_11_2.set_protein_binding_rate(.6  * 10);
  soma.set_transcription_rate(3.*200/10000*0.001*3600  * 2);
  // ds.set_translation_rate(0.021*3600*10);

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

    ofs_covariances << "t=" << t << endl
                    << "covariances:\n" << (covariances = ae.get_covariances()) << endl;
    
    ofs_correlations << "t=" << t << std::endl;
    for(size_t i=0; i<dim; ++i) {
      for(size_t j=0; j<dim; ++j)
        ofs_correlations << covariances(i,j) - expectations(i)*expectations(j) << ',';
      ofs_correlations << std::endl;
    }        

    ofs_variances << t;
    for(size_t i=0; i<dim; ++i)
      ofs_variances  << ',' << sqrt(covariances(i,i) - expectations(i)*expectations(i));
    ofs_variances << endl;

    ae.nonstationary_covariances_direct_ODE_solver_step(dt);
  }

  ofs_expectations.close();
  ofs_covariances.close();
  ofs_variances.close();

  ae.stationary_expectations().stationary_covariances();
  
  return 0;
}
