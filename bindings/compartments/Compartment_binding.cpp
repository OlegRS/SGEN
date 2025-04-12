#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../include/compartments/Compartment.hpp"

namespace py = pybind11;

void bind_Compartment(py::module& m) {
  py::class_<Compartment>(m, "Compartment")
    .def("name", &Compartment::get_name,
         "Return compartment name")
    .def("mRNA_count", &Compartment::mRNA_count,
         "Return the current mRNA count obtained in simulation")
    .def("protein_count", &Compartment::protein_count,
         "Return the current protein count obtained in simulation")
    .def("radius", &Compartment::radius,
         "Return compartment radius")
    .def("position", &Compartment::position,
         "Return {x,y,z}-coordinates of the compartment end point")
    .def("orientation", &Compartment::orientation,
         "Return {theta,phi}-orientation of the compartment cylinder")
    .def("length", &Compartment::get_length,
         "Return compartment length")
    .def("set_radius", &Compartment::set_radius,
         "Set radius of the compartment");
}
