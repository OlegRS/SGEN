#include "../../include/compartments/Soma.hpp"

std::ostream& operator<<(std::ostream& os, const Soma& soma) {
  if(!soma.p_neuron)
    os << "\n-------------------------------------------------\n"
       << "WARNING: Soma is not associated to a neuron, \n"
       << "         so junctions may be uninitialised.\n"
       << "         Use SGEN_Py.Neuron(.) on this object to\n"
       << "         create neuron with this soma.\n"
       << "-------------------------------------------------\n";
    
  os << "Name: " << soma.get_name() << std::endl
     << "- Expected number of active genes: " << soma.n_active_genes_expectation << std::endl
     << "- Expected number of mRNAs: " << soma.n_mRNA_expectation << std::endl
     << "- Expected number of proteins: " << soma.n_prot_expectation << std::endl
     << "- out_junctions:\n";
  for(auto& it_p_out_junc : soma.it_p_out_junctions)
    os << **it_p_out_junc << std::endl;
    
  return os;
}
