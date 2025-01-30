#include <vector>
#include "../Neuron.hpp"
#include "../compartments/Soma.hpp"
#include "../compartments/Dendritic_segment.hpp"
#include "../compartments/Spine.hpp"
#include "../engines/analytic/Analytic_engine.hpp"
#include "../engines/stochastic/Gillespie_engine.hpp"

using namespace std;

#define N_FORKS 1 // Note that it is (2^N_FORKS - 1)*3 compartments (if 2 synapses on each dend seg)!
void fork_dendrite(Dendritic_segment* ds, size_t depth=0) {
  if (depth < N_FORKS) {
    auto ds1 = new Dendritic_segment(*ds, ds->get_name() + "-1");
    new Spine(*ds1, "s_" + ds1->get_name() + "_1");
    new Spine(*ds1, "s_" + ds1->get_name() + "_2", .6, 6);
    fork_dendrite(ds1, depth+1);

    auto ds2 = new Dendritic_segment(*ds, ds->get_name() + "-2");
    new Spine(*ds2, "s_" + ds2->get_name() + "_1");
    new Spine(*ds2, "s_" + ds2->get_name() + "_2", .6, 6);
    fork_dendrite(ds2, depth+1);
  }
}

int main() {

  Soma soma("soma" /*,Parameters of the soma*/);

  ///// Branching neuron
  Dendritic_segment* p_ds = new Dendritic_segment(soma, "d_1");
  Spine *p_syn_1_1 = new Spine(*p_ds, "s_1_1", .6, 6);
  Spine *p_syn_1_2 = new Spine(*p_ds, "s_1_2", .6, 6);
  // fork_dendrite(p_ds);


  Neuron neuron(soma, "Test_neuron");
  
  cout << neuron << endl;

  // cout << "----------------- ANALYTIC ENGINE -----------------\n";
  Analytic_engine ae(neuron);

#define dim 7
  
  arma::mat covariances(dim,dim);
  arma::vec expectations(dim), variances(dim);

  ofstream ofs_expectations("expectations"),
    ofs_covariances("covariances"),
    ofs_variances("variances"),
    ofs_correlations("correlations");

#define t1 5000
#define t2 10000
  double dt = .1;

  // ae.stationary_expectations().stationary_covariances();
  
  std::cerr << "------------------- Loop_1 -----------------------\n";
  for(double t=0; t<t1; t+=dt) {
    ofs_expectations << t << ',' << (expectations = ae.get_expectations()).t() << endl;
    
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

  ae.stationary_expectations().stationary_covariances();

  // Gillespie simulation
  PRNG rnd(1);

  std::string file_name = "gillespie.csv";
  std::list<double> times;
  for(double t=0; t<t1; t+=dt)
    times.push_back(t);

  std::cerr << "Writing Gillespie results to\n" << file_name << '\n';
  
  std::ofstream ofs_gillespie(file_name);
  Gillespie_engine(neuron, rnd).run_Gillespie(times, ofs_gillespie);
  

  std::cerr << "------------------- Loop_2 -----------------------\n";

  p_syn_1_2 -> set_protein_binding_rate(1.2);
  // p_ds -> set_translation_rate(0.021*3600*10);

  neuron.refresh();
  std::cerr << "neuron:\n" << neuron << std::endl;

  for(double t=t1; t<t2; t+=dt) {

    ofs_expectations << t << ',' << (expectations = ae.get_expectations()).t() << endl;
    
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

  // Gillespie simulation
  std::list<double> times_;
  for(double t=t1; t<t2; t+=dt)
    times_.push_back(t);

  std::cerr << "Writing Gillespie results to\n" << file_name << '\n';

  Gillespie_engine(neuron, rnd).run_Gillespie(times_, ofs_gillespie);
  
  ofs_gillespie.close();
  
  return 0;
}
