#include "../../../include/engines/morphologic/Morphologic_engine.hpp"

std::vector<std::vector<std::vector<double>>> Morphologic_engine::segments() {

  std::vector<std::vector<std::vector<double>>> segs(p_neuron->p_dend_segments.size() + p_neuron->p_synapses.size() + 1, std::vector<std::vector<double>>(2, std::vector<double>(4, 0)));

  for(auto& p_junc : p_neuron->p_junctions) {
    segs[p_junc->p_from->id][0] = {p_junc->p_from->x,p_junc->p_from->y,p_junc->p_from->z,p_junc->p_from->r};
    segs[p_junc->p_from->id][1] = {p_junc->p_to->x,p_junc->p_to->y,p_junc->p_to->z,p_junc->p_to->r};

    // Handling the ends of the tree
    if(p_junc->p_to->p_descendants.empty()) {
      auto& comp = *p_junc->p_to;
      segs[comp.id][0] = {comp.x, comp.y, comp.z, comp.r};
      segs[comp.id][1] = {comp.x + comp.length*sin(comp.theta)*cos(comp.phi),comp.y + comp.length*sin(comp.theta)*sin(comp.phi),comp.z + comp.length*cos(comp.theta),comp.r};
    }
  }

  return segs;
}

std::vector<double> Morphologic_engine::volumes() {

  std::vector<double> vols(p_neuron->p_dend_segments.size() + p_neuron->p_synapses.size() + 1);

  for(auto& p_junc : p_neuron->p_junctions) {
    vols[p_junc->p_from->id] = PI * p_junc->p_from->r*p_junc->p_from->r * p_junc->p_from->length;

    // Handling the ends of the tree
    if(p_junc->p_to->p_descendants.empty())
      vols[p_junc->p_to->id] = PI * p_junc->p_to->r*p_junc->p_to->r * p_junc->p_to->length;
  }

  return vols;
}
