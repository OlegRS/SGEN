#include <pybind11/pybind11.h>
#include <sstream>
#include "../include/Neuron.hpp"
#include "../include/compartments/Soma.hpp"

namespace py = pybind11;

// Helper function to invoke std::ostream <<
std::string to_string(const Neuron& neuron) {
    std::ostringstream oss;
    oss << neuron;  // Uses the overloaded << operator
    return oss.str();
}

void bind_Neuron(py::module& m) {    
  py::class_<Neuron>(m, "_Neuron")
    .def(py::init<Soma&, const std::string&>(),
         py::arg("soma"), py::arg("name")="no_name")
    .def(py::init<const std::string&, const std::string&>(),
         py::arg("swc_file_name"), py::arg("name")="no_name",
         "Load a neuron from .swc file")
    .def("soma", &Neuron::soma, py::return_value_policy::reference,
         "Return soma of the neuron")
    .def("refresh", &Neuron::refresh,
         "Update the parameters of the neuron with the current parameters of the compartments")
    .def("__str__", &to_string); // Link Python's str() to C++'s << 
}
