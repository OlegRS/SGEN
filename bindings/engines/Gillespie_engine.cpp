#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // Enables automatic std::vector conversion
#include "../../include/engines/stochastic/Gillespie_engine.hpp"
#include "../../include/randomisation/PRNG.hpp"

namespace py = pybind11;

void bind_Gillespie_engine(py::module& m) {
    py::class_<Gillespie_engine>(m, "_Gillespie_engine")
      .def(py::init<Neuron&>(),
           py::arg("neuron"),
           "Initialise Gillespie Engine for a given neuron.")
      .def("run_Gillespie", py::overload_cast<const std::vector<double>&, const std::string&, const double&>(&Gillespie_engine::run_Gillespie),
           py::arg("record_times"), py::arg("file_name")="Gillespie_out.csv", py::arg("burn_in")=0,
           "Run Gillespie algorithm writing at record_times to file_name with offset time_offset")
      .def("variable_names", &Gillespie_engine::variable_names,
           "Return ordered array of variable names");
}
