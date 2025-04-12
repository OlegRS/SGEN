#include <vector>
#include "../Neuron.hpp"
#include "../compartments/Soma.hpp"
#include "../compartments/Dendritic_segment.hpp"
#include "../compartments/Spine.hpp"
#include "../engines/stochastic/Gillespie_engine.hpp"

#define N_AVRG  10000
#define event_time 5000

#define syn_dec_rate 0 //1.2e-5*3600
#define transcription_rate_multiplier 2
#define binding_rate_multiplier 10
#define SYN syn_1_1
std::string syn_name = "syn_1_1";


int main() {

  PRNG rnd(1);

  double dt = .01;  
  
  std::string file_name = "../../data/gillespie/2_genes_heterosynaptic_plasticity_full/HP_" + syn_name + "_transcription_rate_multiplier_" + std::to_string(transcription_rate_multiplier) + "_bind_rate_multiplier_" + std::to_string(binding_rate_multiplier) + '_';

  std::list<double> times;
  for(double t=0; t<event_time; t+=dt)
    times.push_back(t);
  
  for(size_t i=0; i<N_AVRG; ++i) {

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

    std::ofstream ofs_gillespie(file_name + std::to_string(i));

    std::cerr << "Writing Gillespie results to: " << file_name + std::to_string(i) << '\n';
  
    std::cerr << "------------------- Loop_1 -----------------------\n";
    std::cout << neuron << std::endl;
 
    Gillespie_engine(neuron, rnd).run_Gillespie(times, ofs_gillespie);
  
    std::cerr << "------------------- Loop_2 -----------------------\n";

    SYN.set_protein_binding_rate(.6*binding_rate_multiplier);
    soma.set_transcription_rate(3.*200/10000*0.001*3600*transcription_rate_multiplier);
    neuron.refresh();

    std::cout << neuron << std::endl;
    
    Gillespie_engine(neuron, rnd).run_Gillespie(times, ofs_gillespie, event_time);
  
    ofs_gillespie.close();
  }
  
  return 0;
}
