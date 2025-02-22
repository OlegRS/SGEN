#include <math.h>
#include "../../../include/engines/stochastic/Gillespie_engine.hpp"

Compartment* Gillespie_engine::initialise_soma() {

  auto& soma = *p_neuron->p_soma;

  p_events[0] = &soma.gene_activation.set_rate(soma.gene_activation_rate*(soma.number_of_gene_copies-soma.n_active_genes));
  p_neuron->total_rate = soma.gene_activation.rate;
  p_events[1] = &soma.gene_deactivation.set_rate(soma.gene_deactivation_rate*soma.n_active_genes);
  p_neuron->total_rate += soma.gene_deactivation.rate;
  p_events[2] = &soma.mRNA_creation.set_rate(soma.transcription_rate*soma.n_active_genes);
  p_neuron->total_rate += soma.mRNA_creation.rate;
  p_events[3] = &soma.mRNA_decay.set_rate(soma.mRNA_decay_rate*soma.n_mRNAs);
  p_neuron->total_rate += soma.mRNA_decay.rate;
  p_events[4] = &soma.protein_creation.set_rate(soma.translation_rate*soma.n_mRNAs);
  p_neuron->total_rate += soma.protein_creation.rate;
  p_events[5] = &soma.protein_decay.set_rate(soma.protein_decay_rate*soma.n_proteins);
  p_neuron->total_rate += soma.protein_decay.rate;
  
  ev_ind = 6;

  return &soma;
}

void Gillespie_engine::initialise_from(Compartment& comp) {

  for(auto& it_p_junc : comp.it_p_out_junctions) {
    auto& p_junc = *it_p_junc;
    auto& par = *(p_junc->p_from); // Parent
    auto& desc = *(p_junc->p_to);  // Descendant

    if(p_junc->type() == DEN_SYN) {
      p_events[ev_ind++] = &p_junc->prot_hop_forward.set_rate(p_junc->fwd_prot_hop_rate*par.n_proteins);
      p_neuron->total_rate += p_junc->prot_hop_forward.rate;
      p_events[ev_ind++] = &p_junc->prot_hop_backward.set_rate(p_junc->bkwd_prot_hop_rate*desc.n_proteins);
      p_neuron->total_rate += p_junc->prot_hop_backward.rate;

      // No mRNAs creation/decay and protein creation in synapses (descendant is a synapse)
      desc.mRNA_creation.rate = 0;
      desc.mRNA_decay.rate = 0;
      desc.protein_creation.rate = 0;
      p_events[ev_ind++] = &desc.protein_decay.set_rate(desc.protein_decay_rate*desc.n_proteins);
      p_neuron->total_rate += desc.protein_decay.rate; 
    }
    else if(p_junc->type() == DEN_DEN || p_junc->type() == SOM_DEN) {
      p_events[ev_ind++] = &p_junc->mRNA_hop_forward.set_rate(p_junc->fwd_mRNA_hop_rate*par.n_mRNAs);
      p_neuron->total_rate += p_junc->mRNA_hop_forward.rate;
      p_events[ev_ind++] = &p_junc->mRNA_hop_backward.set_rate(p_junc->bkwd_mRNA_hop_rate*desc.n_mRNAs);
      p_neuron->total_rate += p_junc->mRNA_hop_backward.rate;
      p_events[ev_ind++] = &p_junc->prot_hop_forward.set_rate(p_junc->fwd_prot_hop_rate*par.n_proteins);
      p_neuron->total_rate += p_junc->prot_hop_forward.rate;
      p_events[ev_ind++] = &p_junc->prot_hop_backward.set_rate(p_junc->bkwd_prot_hop_rate*desc.n_proteins);
      p_neuron->total_rate += p_junc->prot_hop_backward.rate;

      desc.mRNA_creation.rate = 0;
      p_events[ev_ind++] = &desc.mRNA_decay.set_rate(desc.mRNA_decay_rate*desc.n_mRNAs);
      p_neuron->total_rate += desc.mRNA_decay.rate;
      p_events[ev_ind++] = &desc.protein_creation.set_rate(desc.translation_rate*desc.n_mRNAs);
      p_neuron->total_rate += desc.protein_creation.rate;
      p_events[ev_ind++] = &desc.protein_decay.set_rate(desc.protein_decay_rate*desc.n_proteins);
      p_neuron->total_rate += desc.protein_decay.rate;
    }
    else {
      std::cerr << "-----------------------------------------\n"
                << "ERROR: UNKNOWN TYPE JUNCTION FOUND\n"
                << "-----------------------------------------\n";
      exit(1);
    }
    
    initialise_from(*((*it_p_junc)->p_to));
  }
}

inline double Gillespie_engine::draw_delta_t() {  
  // Sampling the time of an event
  double delta_t = -log(1-rnd())/p_neuron->total_rate;

  return delta_t; 
}

// std::default_random_engine generator;

void Gillespie_engine::update_Gillespie() {  
  // Sampling the event
  
  // std::vector<double> probabilities(p_events.size());
  // for(size_t i=0; i<p_events.size(); ++i)
  //   probabilities[i] = p_events[i]->rate;

  // std::discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());          
  // int i = distribution(generator);

  // double s=0;
  // for(auto& p_event : p_events)
  //   s+=p_event->rate;
  // if(abs(p_neuron->total_rate - s) > 1e-5) {
  //   std::cerr << "total_rate - s = " << p_neuron->total_rate - s << std::endl;
  // }

  size_t i=0;
  double r = rnd()*p_neuron->total_rate;
  for(double sum=p_events[0]->rate; sum<r; sum += p_events[i]->rate)
    ++i;
  if(i<p_events.size())
    (*p_events[i])(); // Triggering the event
  else { // Recalculate total rate when numerical errors accumulate
    std::cerr <<"\n-----------------------------------------\n"
              <<" WARNING: GILLESPIE ENGINE ANOMALY DETECTED"
              <<"\n-----------------------------------------\n";
    double s=0;
    for(auto& p_event : p_events)
      s+=p_event->rate;
    std::cerr << "total_rate - s = " << p_neuron->total_rate - s << std::endl;
  }
}

Gillespie_engine& Gillespie_engine::run_Gillespie(const double& time) {
  
  // FOR PRELIMINARY TESTING ONLY!!!
  std::cout << "t," << "Soma_AG," << "Soma_mRNA,"<< "Soma_Prot,";
  for(auto& p_ds : p_neuron->p_dend_segments)
    std::cout << p_ds->name + "_mRNA," << p_ds->name+"_Prot,";
  for(auto& p_s : p_neuron->p_synapses)
    std::cout << p_s->name + "_Prot,";
  std::cout << std::endl;

  double t_prev=-1.1, delta_t_write = 1;
  for(double t = 0; t<time; t+=draw_delta_t()) {
    if(t-t_prev > delta_t_write) {
      std::cout << t << ',' << p_neuron->p_soma->n_active_genes << ',' << p_neuron->p_soma->n_mRNAs << ',' << p_neuron->p_soma->n_proteins << ',';
      for(auto& p_ds : p_neuron->p_dend_segments)
        std::cout << p_ds->n_mRNAs << ',' << p_ds->n_proteins << ',';
      for(auto& p_s : p_neuron->p_synapses)
        std::cout << p_s->n_proteins << ',';
      std::cout << std::endl;
      t_prev = t;
    }
    update_Gillespie();
  }
  return *this;
}

Gillespie_engine& Gillespie_engine::run_Gillespie(const std::list<double>& times, std::ostream& os, const double& time_offset) {
  std::cout << "Running Gillespie...\n";
  size_t n_jumps = 0;
  
  double t = 0;
  double time_written;
  for(auto it_times=times.begin(); it_times!=times.end(); ) {
    if(*it_times > t) {
      t += draw_delta_t();
      while(*it_times <= t && it_times!=times.end()) {
        os << t << ',' << time_offset + (time_written=*it_times) << ',' << p_neuron->p_soma->n_active_genes << ',' << p_neuron->p_soma->n_mRNAs << ',' << p_neuron->p_soma->n_proteins;
        for(auto& p_ds : p_neuron->p_dend_segments)
          os << ',' << p_ds->n_mRNAs << ',' << p_ds->n_proteins;
        for(auto& p_s : p_neuron->p_synapses)
          os << ',' << p_s->n_proteins;
        os << std::endl;
        ++it_times;
      }
      update_Gillespie();
           
      n_jumps++;
    }
    else {
      if(time_written != *it_times) {
        os << t << ',' << time_offset + *it_times << ',' << p_neuron->p_soma->n_active_genes << ',' << p_neuron->p_soma->n_mRNAs << ',' << p_neuron->p_soma->n_proteins;
        for(auto& p_ds : p_neuron->p_dend_segments)
          os << ',' << p_ds->n_mRNAs << ',' << p_ds->n_proteins;
        for(auto& p_s : p_neuron->p_synapses)
          os << ',' << p_s->n_proteins;
        os << std::endl;
      }
      ++it_times;
    }
  }
  std::cout << "n_jumps = " << n_jumps << std::endl;
  return *this;
}

std::vector<std::vector<double>> Gillespie_engine::run_Gillespie(const std::vector<double>& times, const std::string& file_name, const double& time_offset) {
  std::vector<std::vector<double>> results(times.size(), std::vector<double>(dim+1)); //+1 for time variable
  
  for(auto it_times=times.begin(); it_times!=times.end()-1;) {
    if(*it_times < 0) {
      std::cerr << "\n--------------------------------------------------------\n"
                << "- ERROR in run_Gillespie: record_times should be positive!"
                << "\n--------------------------------------------------------\n";
      return results;
    }
    if(*it_times == *(it_times+1)) {
      std::cerr << "\n------------------------------------------------------------\n"
                << "- WARNING from run_Gillespie: there are repeated record_times!"
                << "\n------------------------------------------------------------\n";
      return results;
    }
    if(*it_times > *++it_times) {
      std::cerr << "\n--------------------------------------------------------------------------\n"
                << "- ERROR in run_Gillespie: record_times should be ordered in ascending order!"
                << "\n--------------------------------------------------------------------------\n";
      return results;
    }
  }

  size_t n_jumps = 0;
  
  if(file_name.empty()) {
    double t = 0;
    for(size_t t_ind=0; t_ind < times.size(); ) {
      if(times[t_ind] >= t) {
        t += draw_delta_t();
        while(t_ind != times.size() && times[t_ind] <= t) {
          write_results(time_offset + times[t_ind], results[t_ind]);
          ++t_ind;
        }

        if(t_ind != times.size()) {// This would rarely not be the case with high event rates
          update_Gillespie();           
          ++n_jumps;
        }
      }
      else {// This should never happen. Keeping for debugging as it does not affect performance.
        std::cerr << "\n----------------------------------------------------------------\n"
                  << "- WARNING: Something weird just happened with Gillespie algorithm!"
                  << "\n----------------------------------------------------------------\n";
        ++t_ind;
      }
    }
  }
  else {
    std::ofstream ofs(file_name);
    print_variable_names(ofs << "t,") << std::endl;
    
    double t = 0;
    for(size_t t_ind=0; t_ind < times.size(); ) {
      if(times[t_ind] >= t) {
        t += draw_delta_t();
        while(t_ind != times.size() && times[t_ind] <= t) {
          print_variables(ofs << time_offset + times[t_ind] << ',') << std::endl;
          write_results(time_offset + times[t_ind], results[t_ind]);
          ++t_ind;
        }

        if(t_ind != times.size()) {// This would rarely not be the case with high event rates
          update_Gillespie();           
          ++n_jumps;
        }
      }
      else {// This should never happen. Keeping for debugging as it does not affect performance.
        std::cerr << "\n----------------------------------------------------------------\n"
                  << "- WARNING: Something weird just happened with Gillespie algorithm!"
                  << "\n----------------------------------------------------------------\n";
        ++t_ind;
      }
    }
    ofs.close();
  }

  std::cout << "n_jumps = " << n_jumps << std::endl;
  
  return results;
}

inline void Gillespie_engine::write_results(const double& time, std::vector<double> &results) {
  
  results[0] = time;
  results[1] = p_neuron->p_soma->n_active_genes;
  results[2] = p_neuron->p_soma->n_mRNAs;
  results[3] = p_neuron->p_soma->n_proteins;

  size_t i=4;

  for(auto& p_ds : p_neuron->p_dend_segments) {
    results[i++] = p_ds->n_mRNAs;
    results[i++] = p_ds->n_proteins;
  }
  
  for(auto& p_s : p_neuron->p_synapses)
    results[i++] = p_s->n_proteins;
}

inline std::ostream& Gillespie_engine::print_variables(std::ostream& os) {
  os << p_neuron->p_soma->n_active_genes << ',' << p_neuron->p_soma->n_mRNAs << ',' << p_neuron->p_soma->n_proteins;
  for(auto& p_ds : p_neuron->p_dend_segments)
    os << ',' << p_ds->n_mRNAs << ',' << p_ds->n_proteins;
  for(auto& p_s : p_neuron->p_synapses)
    os << ',' << p_s->n_proteins;

  return os;
}

inline std::ostream& Gillespie_engine::print_variable_names(std::ostream& os) {
  os << "Soma_AG," << "Soma_mRNA,"<< "Soma_Prot";
  for(auto& p_ds : p_neuron->p_dend_segments)
    os << ',' << p_ds->name + "_mRNA," << p_ds->name+"_Prot";
  for(auto& p_s : p_neuron->p_synapses)
    os << ',' << p_s->name + "_Prot";
  
  return os;
}

std::vector<std::string> Gillespie_engine::variable_names() {

  std::vector<std::string> var_names(dim);

  var_names[0] = p_neuron->p_soma->name + "_gene";
  var_names[1] = p_neuron->p_soma->name + "_mRNA";
  var_names[2] = p_neuron->p_soma->name + "_prot";

  size_t i=3;

  for(auto& p_ds : p_neuron->p_dend_segments) {
    var_names[i++] = p_ds->name + "_mRNA";
    var_names[i++] = p_ds->name + "_prot";
  }
  
  for(auto& p_s : p_neuron->p_synapses)
    var_names[i++] = p_s->name + "_prot";

  
  return var_names;
}
