#include <pybind11/pybind11.h>
#include "../../include/compartments/Compartment.hpp"

namespace py = pybind11;

void bind_Compartment(py::module& m) {
  py::class_<Compartment>(m, "Compartment");
}
