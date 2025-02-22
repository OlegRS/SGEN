#include <pybind11/pybind11.h>
#include <sstream>
#include "../../include/compartments/Soma.hpp"

namespace py = pybind11;

// Helper function to invoke std::ostream <<
std::string to_string(const Soma& soma) {
    std::ostringstream oss;
    oss << soma;  // Uses the overloaded << operator
    return oss.str();
}

void bind_Soma(py::module& m) {
  py::class_<Soma, Compartment>(m, "Soma")
    .def(py::init<const std::string&, const double&>(),
         py::arg("name")="no_name_Soma", py::arg("length")=10,
         "Create soma with certain name and length")
    .def(py::init<const std::string&, const double&, const double&, const double&, const double&, const double&, const size_t&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&>(),
         py::arg("name")="no_name_Soma", py::arg("length")=10, py::arg("x")=0, py::arg("y")=0, py::arg("z")=0, py::arg("radius")=5, py::arg("n_gene_copies")=1, py::arg("gene_activation_rate")=1/12., py::arg("gene_deactivation_rate")=1/12., py::arg("mRNA_decay_rate")=0.0432, py::arg("translation_rate")=75.6, py::arg("protein_decay_rate")=0.004356, py::arg("mRNA_diffusion_constant")=1.08e-5, py::arg("protein_diffusion_constant")=7.6e-4, py::arg("mRNA_forward_trafficking_velocity")=1.6e-5, py::arg("mRNA_backward_trafficking_velocity")=3.18, py::arg("protein_forward_trafficking_velocity")=0, py::arg("protein_backward_trafficking_velocity")=0,
         "Create soma with certain name and length, location (x,y,z) and radius")
    .def("gene_activation_rate", &Soma::get_gene_activation_rate)
    .def("gene_deactivation_rate", &Soma::get_gene_deactivation_rate)
    .def("__str__", &to_string);
}
