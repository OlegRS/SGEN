#include <pybind11/pybind11.h>
#include "../../include/compartments/Spine.hpp"

namespace py = pybind11;

void bind_Spine(py::module& m) {
  py::class_<Spine, Compartment>(m, "Spine")
    .def(py::init<Compartment&, const std::string&, const double&, const double&, const double&, const double&>(),
         py::arg("parent"), py::arg("name")="no_name", py::arg("length")=2, py::arg("radius")=1, py::arg("d_theta")=PI/2, py::arg("d_phi")=0,
         "Create dendritic spine with certain name and length, radius and angle (d_theta, d_phi)");
}
