#include <pybind11/pybind11.h>
#include "../../include/compartments/Dendritic_segment.hpp"

namespace py = pybind11;

void bind_Dendritic_segment(py::module& m) {
  py::class_<Dendritic_segment, Compartment>(m, "Dendritic_segment")
    .def(py::init<Compartment&, const std::string&, const double&,  const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const bool&>(),
         py::arg("parent"), py::arg("name")="no_name_Dendritic_segment", py::arg("length")=10, py::arg("radius")=5, py::arg("d_theta")=0, py::arg("d_phi")=0, py::arg("mRNA_decay_rate")=0.0432, py::arg("translation_rate")=75.6, py::arg("protein_decay_rate")=0.004356, py::arg("mRNA_diffusion_constant")=3.4e-3, py::arg("protein_diffusion_constant")=.24, py::arg("mRNA_forward_trafficking_velocity")=.5e-2, py::arg("mRNA_backward_trafficking_velocity")=.1e-2, py::arg("protein_forward_trafficking_velocity")=0, py::arg("protein_backward_trafficking_velocity")=0, py::arg("middle_placement")=false, 
         "Create dendritic segment with specified parameters")
    .def("translation_rate", &Dendritic_segment::translation_rate,
         "Returns translation rate in the compartment [1/hour]")
    .def("set_translation_rate", &Dendritic_segment::set_translation_rate,
         py::arg("translation_rate")=75.6,
         "Sets translation rate in the compartment [1/hour]");
}
