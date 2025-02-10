#include <pybind11/pybind11.h>

namespace py = pybind11;

// Declare binding functions for different components
void bind_Compartment(py::module&);
void bind_Neuron(py::module&);
void bind_Soma(py::module&);
void bind_Dendritic_segment(py::module&);
void bind_Spine(py::module&);
void bind_Analytic_engine(py::module&);
void bind_Morphologic_engine(py::module&);
void bind_Gillespie_engine(py::module&);

PYBIND11_MODULE(SGEN_Py, m) {
    m.doc() = "Python bindings for the Stochastic Gene Expression in Neurons (SGEN) library";

    // Bind components
    bind_Compartment(m);
    bind_Soma(m);
    bind_Dendritic_segment(m);
    bind_Spine(m);
    bind_Neuron(m);
    bind_Analytic_engine(m);
    bind_Morphologic_engine(m);
    bind_Gillespie_engine(m);
}
