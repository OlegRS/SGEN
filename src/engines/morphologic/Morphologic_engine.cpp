#include "../../../include/engines/morphologic/Morphologic_engine.hpp"

std::vector<std::vector<std::vector<double>>> Morphologic_engine::segments() {

  std::vector<std::vector<std::vector<double>>> segs(p_neuron->p_dend_segments.size() + p_neuron->p_synapses.size() + 1, std::vector<std::vector<double>>(2, std::vector<double>(4, 0)));

  Soma &soma = *p_neuron->p_soma;
  segs[0][0] = {soma.x, soma.y, soma.z-soma.length, soma.r};
  segs[0][1] = {soma.x, soma.y, soma.z, soma.r};

  for(auto& p_junc : p_neuron->p_junctions) {
    segs[p_junc->p_to->id][0] = {p_junc->p_from->x,p_junc->p_from->y,p_junc->p_from->z,p_junc->p_from->r};
    segs[p_junc->p_to->id][1] = {p_junc->p_to->x,p_junc->p_to->y,p_junc->p_to->z,p_junc->p_to->r};
  }

  return segs;
}

std::vector<double> Morphologic_engine::volumes() {

  std::vector<double> vols(p_neuron->p_dend_segments.size() + p_neuron->p_synapses.size() + 1);

  vols[0] = p_neuron->p_soma->cross_section() * p_neuron->p_soma->length;
  
  for(auto& p_junc : p_neuron->p_junctions)
    vols[p_junc->p_to->id] = p_junc->p_to->cross_section() * p_junc->p_to->length;

  return vols;
}
