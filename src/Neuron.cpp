#include "../include/Neuron.hpp"
#include "../include/compartments/Soma.hpp"
#include "../include/compartments/Dendritic_segment.hpp"
#include "../include/compartments/Spine.hpp"
#include "../include/junctions/Som_den_junction.hpp"
#include "../include/junctions/Den_den_junction.hpp"
#include "../include/junctions/Den_syn_junction.hpp"

void Neuron::associate(Compartment& compartment) {

  compartment.p_neuron = this;

  Compartment::Type p_type = compartment.type();
  compartment.id = comp_id++;
  
  if(p_type == SPINE) {
    p_synapses.push_back(&compartment);
    compartment.iterator = p_synapses.end();
    compartment.o1_index = o1_index;
    o1_index += 1;
  }
  else if(p_type == APICAL_DENDRITE || p_type == BASAL_DENDRITE) {
    p_dend_segments.push_back(&compartment);
    compartment.iterator = p_dend_segments.end();
    compartment.o1_index = o1_index;
    o1_index += 2;
    compartment.mRNA_ind = mRNA_ind++;
  }
  else if(p_type == SOMA) {
    compartment.o1_index = o1_index;
    o1_index += 3;
    compartment.mRNA_ind = mRNA_ind++;
  }

  // Writing junctions
  for(auto& p_d_comp : compartment.p_descendants){//Loop over descendants
    auto d_type = p_d_comp->type();
    bool p_dend = p_type == APICAL_DENDRITE || p_type == BASAL_DENDRITE;
    bool d_dend = d_type == APICAL_DENDRITE || d_type == BASAL_DENDRITE;
    if(p_type == SOMA && d_dend) {
      p_junctions.push_back(new Som_den_junction(&compartment, p_d_comp));
      compartment.it_p_out_junctions.push_back(--p_junctions.end());
      p_d_comp->it_p_in_junctions.push_back(--p_junctions.end());
      ++n_SDJ;
    }
    else if(p_dend && d_type==SPINE) {
      p_junctions.push_back(new Den_syn_junction(&compartment, p_d_comp));
      compartment.it_p_out_junctions.push_back(--p_junctions.end());
      p_d_comp->it_p_in_junctions.push_back(--p_junctions.end());
      ++n_DSJ;
    }
    else if(p_dend && d_dend) {
      p_junctions.push_back(new Den_den_junction(&compartment, p_d_comp));
      compartment.it_p_out_junctions.push_back(--p_junctions.end());
      p_d_comp->it_p_in_junctions.push_back(--p_junctions.end());
      ++n_DDJ;
    }
    else {
      std::cerr << "---------------------------------------------------\n"
                << "ERROR: " << d_type << " descending from " << p_type << std::endl
                << "---------------------------------------------------\n";
      exit(1);
    }

    associate(*p_d_comp);
  }
}

// void Neuron::dissociate(Compartment& compartment) {
//   compartment.p_neuron = nullptr;

//   Compartment::Type type = compartment.type();

//   if(type.id == SPINE) 
//     p_synapses.erase(compartment.iterator);
//   else if(type.id == APICAL_DENDRITE || type.id == BASAL_DENDRITE)
//     p_dend_segments.erase(compartment.iterator);
  
//   for (auto& comp : compartment.p_descendants)
//     dissociate(*comp);
// }

Neuron& Neuron::refresh() {

  clear_junctions();

  p_dend_segments.clear();
  p_synapses.clear();
  
  o1_index=0; mRNA_ind=1; comp_id=0;

  associate(static_cast<Compartment&>(*p_soma));
  
  return *this;
}

Soma& Neuron::soma() {
  if(!p_soma) {
    std::cerr << "-----------------------------------------------\n"
              << "- ERROR: Soma of uninitialised neuron requested\n"
              << "-----------------------------------------------\n";
    exit(1);
  }

  return *p_soma;
}

Neuron::Neuron(Soma &soma, const std::string &name) : name(name) {
  p_soma = &soma;
  associate(static_cast<Compartment&>(soma));
}




Neuron::Neuron(const std::string& file_name, const std::string& name) : name(name) {
    constexpr size_t OFFSET = 2;  // Lines to skip in the SWC file
    std::ifstream ifs(file_name);

    if (!ifs.is_open()) {
        std::cerr << "------------------------------------------\n"
                  << "ERROR: File not found: " << file_name << "\n"
                  << "------------------------------------------\n";
        std::exit(1);
    }

    ifs.seekg(-2, std::ios_base::end);
    bool keepLooping = true;
    while(keepLooping) {
      char ch;
      ifs.get(ch);

      if(ch == '\n')
        keepLooping = false;
      else
        ifs.seekg(-2, std::ios_base::cur);
    }

    //// Obtaining the total number of compartments as the last compartment's id
    size_t n_comp;
    ifs >> n_comp;
    n_comp -= 2;
    ifs.seekg(0);

    std::string line;
    Compartment** p_compartments = new Compartment*[n_comp];
    size_t id, type;
    double x, y, z, r;
    int parent_id;

    size_t line_count = 0;

    // Skip offset lines and read data
    while (std::getline(ifs, line)) {
      std::istringstream line_stream(line);
      if (line[0] == '#') // Skip comments in SWC files
        continue;
        
      if (line_count++ < OFFSET)
        continue;

      line_stream >> id >> type >> x >> y >> z >> r >> parent_id;

      // Assume the coordinates in .swc files are the ends of the compartments
      if (type == SOMA)
        p_compartments[id-OFFSET-1] = p_soma = new Soma("soma_" + std::to_string(id-OFFSET), r, 0, 0, 0, r); // Place soma to the origin and assume it is oriented as (0,0,1) with radius and length r
      else if (type == BASAL_DENDRITE || type == APICAL_DENDRITE) {
        if (parent_id == 1) // If parent is SOMA
          p_compartments[id-OFFSET-1] = new Dendritic_segment(*p_compartments[0], x, y, z, r, "ds_" + std::to_string(id), sqrt(x*x + y*y + z*z));
        else {
          Compartment* p_parent = p_compartments[parent_id-OFFSET-1];
          p_compartments[id-OFFSET-1] = new Dendritic_segment(*p_parent, x, y, z, r, "ds_" + std::to_string(id), sqrt((x-p_parent->x)*(x-p_parent->x)+(y-p_parent->y)*(y-p_parent->y)+(z-p_parent->z)*(z-p_parent->z)));
        }
      }
    }
    
    associate(*p_soma);
    delete[] p_compartments;
}

std::ostream& operator<<(std::ostream &os , const Neuron &neur) {

  os << "******* COMPARTMENTS *******:\n"
     << "* Soma:\n** " << *(Compartment*)neur.p_soma << std::endl
     << "\n* Dendritic segments:\n";

  for(auto& p_ds : neur.p_dend_segments)
    os << "** " << *p_ds << "\n";
  os << '\n';
  os << "* Spines:\n";
  for(auto& p_syn : neur.p_synapses)
    os << "** " << *p_syn << "\n";
  os << '\n';
  os << "******* JUNCTIONS *******:\n";
  for(auto& junct : neur.p_junctions)
    os << "** " << *junct << "\n";
  
  return os;
}

Neuron& Neuron::clear_junctions() {
  p_soma->clear_junctions();
  for(auto p_ds : p_dend_segments)
    p_ds -> clear_junctions();
  for(auto p_s : p_synapses)
    p_s -> clear_junctions();
  for(auto& junct : p_junctions)
    delete junct;
  p_junctions.clear();

  n_SDJ=0; n_DSJ=0; n_DDJ=0; // Numbers of different junctions

  return *this;
}

Neuron::~Neuron() {
  for(auto& junct : p_junctions)
    delete junct;  
}
