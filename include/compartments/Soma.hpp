#ifndef __SOMA_HPP__
#define __SOMA_HPP__

#include "../Neuron.hpp"

class Soma : public Compartment {
  friend class Neuron::Som_den_junction;
  friend class Analytic_engine;
  friend class Gillespie_engine;
  friend class Dendritic_segment;
  
  // Parameters
  size_t number_of_gene_copies = 1; //For CaMKIIa it is 2 (Fonkeu)
  double gene_activation_rate = 1/12.;
  double gene_deactivation_rate = 1/12.;

  // For Monte Carlo engines
  struct Gene_activation : public Event {
    Gene_activation(Soma* p_loc) : Event(p_loc) {}
    Event::Type type() {return GENE_ACTIVATION;}
    void operator()();
  } gene_activation;

  struct Gene_deactivation : public Event {
    Gene_deactivation(Soma* p_loc) : Event(p_loc) {}
    Event::Type type() {return GENE_DEACTIVATION;}
    void operator()();
  } gene_deactivation;

  double n_active_genes_expectation = gene_activation_rate/(gene_activation_rate + gene_deactivation_rate)*number_of_gene_copies;
  double n_active_genes_variance = n_active_genes_expectation*(1+gene_activation_rate/(gene_activation_rate + gene_deactivation_rate)*(number_of_gene_copies-1));
  size_t n_active_genes = 0;
  size_t n_descending_DS = 0;

  double transcription_rate = (3.*200/*dend_length*//10000)*.001*3600; // /hour; mRNA transcription rate (0.001/s CaMKII Fonkeu) // THE FACTOR IN () ACCOUNTS FOR THE REDUCED LENGTH OF THE SIMPLE MODEL DENDRITE COMPARED TO THE REAL NEURONS
  
public:
  Soma(const std::string& name="no_name", const double& length=20, const double& x=0, const double& y=0, const double& z=0, const double& r=10, const size_t& n_gene_copies=1, const double& gene_activation_rate=1/12., const double& gene_deactivation_rate=1/12., const double& mRNA_decay_rate=0.0432, const double& transcription_rate=.216, const double& translation_rate=75.6, const double& protein_decay_rate=0.004356, const double& mRNA_diffusion_constant=3.4e-3, const double& protein_diffusion_constant=.24, const double& mRNA_forward_trafficking_velocity=.5e-2, const double& mRNA_backward_trafficking_velocity=.1e-2, const double& protein_forward_trafficking_velocity=0, const double& protein_backward_trafficking_velocity=0) :
    Compartment(length, name, x, y, z, r, /*theta=*/0, /*phi=*/0, mRNA_decay_rate, translation_rate, protein_decay_rate, mRNA_diffusion_constant, protein_diffusion_constant, mRNA_forward_trafficking_velocity, mRNA_backward_trafficking_velocity, protein_forward_trafficking_velocity, protein_backward_trafficking_velocity),
    gene_activation(this), gene_deactivation(this), number_of_gene_copies(n_gene_copies), gene_activation_rate(gene_activation_rate), gene_deactivation_rate(gene_deactivation_rate), transcription_rate(transcription_rate) {}

  Soma(unsigned int &number_of_gene_copies, double &gene_activation_rate, double &gene_deactivation_rate, double &transcription_rate, double &mRNA_decay_rate, double &translation_rate, double &protein_decay_rate, unsigned int active_genes_number = 0, unsigned int protein_number = 0, unsigned int mRNA_number = 0);

  Soma& set_gene_activation_rate(const double& rate) {gene_activation_rate=rate; return *this;}
  Soma& set_gene_deactivation_rate(const double& rate) {gene_deactivation_rate = rate; return *this;}
  Soma& set_number_of_gene_copies(const unsigned int& N) {number_of_gene_copies = N; return *this;}
  Compartment& set_transcription_rate(const double& rate) {transcription_rate=rate; return *this;}

  double get_gene_activation_rate() const {return gene_activation_rate;}
  double get_gene_deactivation_rate() const {return gene_deactivation_rate;}

  double get_transcription_rate() const {return transcription_rate;}

  size_t get_n_gene_copies() const {return number_of_gene_copies;}
  size_t get_n_active_genes() const {return n_active_genes;}
  
  Compartment::Type type() const {return SOMA;}
  
  Compartment* operator+=(Compartment&);

  friend std::ostream& operator<<(std::ostream&, const Soma&);
 
  Compartment* add(Compartment&);
};

#endif
