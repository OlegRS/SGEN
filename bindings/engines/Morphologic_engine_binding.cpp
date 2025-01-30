#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../include/engines/morphologic/Morphologic_engine.hpp"

namespace py = pybind11;

void bind_Morphologic_engine(py::module& m) {
    py::class_<Morphologic_engine>(m, "Morphologic_engine")
      .def(py::init<Neuron&>(),
           py::arg("neuron"),
           "Initialise Morphologic Engine for a given neuron")
      .def("segments", &Morphologic_engine::segments,
           "Return list of all segments (compartments) of the neuron in the form [..., [[x_i,y_i,z_i,r_i],[x_i+1,y_i+1,z_i+1,r_i+1]],...]")
      .def("volumes", &Morphologic_engine::volumes,
           "Return list of volumes of all compartments ordered consistently with segments()");

}
