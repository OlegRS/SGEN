#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // Enables automatic std::vector conversion
#include "../../include/engines/stochastic/Gillespie_engine.hpp"
#include "../../include/randomisation/PRNG.hpp"

namespace py = pybind11;

void bind_Gillespie_engine(py::module& m) {
    py::class_<Gillespie_engine>(m, "Gillespie_engine")
      .def(py::init<Neuron&>(),
           py::arg("neuron"),
           "Initialise Gillespie Engine for a given neuron.")
      .def("run_Gillespie", py::overload_cast<const std::vector<double>&, const std::string&, const double&>(&Gillespie_engine::run_Gillespie),
           py::arg("record_times"), py::arg("file_name")="Gillespie_out.csv", py::arg("time_offset")=0,
           "Run Gillespie algorithm writing at record_times to file_name with offset time_offset");

      // .def("run_Gillespie", &Gillespie_engine::run_Gillespie,
      //      py::arg("record_times"), py::arg("file_name")="Gillespie_out", py::arg("time_offset")=0,
      //      "Run Gillespie algorithm writing at record_times to file_name with offset time_offset");


      
      // .def("run_Gillespie", 
     // [&](Gillespie_engine& self, const std::vector<double>& record_times, const std::string& file_name, double time_offset) {
     //     std::ofstream file_out(file_name);
     //     if (!file_out.is_open()) {
     //         throw std::runtime_error("Failed to open file: " + file_name);
     //     }
     //     return self.run_Gillespie(record_times, file_out, time_offset);
     // },
     // py::arg("record_times"), 
     // py::arg("file_name") = "Gillespie_out", 
     // py::arg("time_offset") = 0,
     // "Run Gillespie algorithm writing at record_times to file_name with offset time_offset.");

      
      
      // .def("run_Gillespie",
      //      [](Gillespie_engine& self, const std::vector<double>& record_times, const std::string& file_name, double time_offset) {
      //            std::ofstream file_out(file_name);
      //            if (!file_out.is_open()) {
      //                throw std::runtime_error("Failed to open file: " + file_name);
      //            }
      //            return self.run_Gillespie(record_times, file_out, time_offset);
      //        },
      //        py::arg("record_times"), 
      //        py::arg("file_name") = "Gillespie_out", 
      //        py::arg("time_offset") = 0,
      //        "Run Gillespie algorithm writing at record_times to file_name with offset time_offset.");
      // .def("run_Gillespie", py::overload_cast<const std::vector<double>&, std::ofstream&, const double&>(&Gillespie_engine::run_Gillespie),
      //      py::arg("record_times"), py::arg("file_name")="Gillespie_out", py::arg("time_offset")=0,
      //      "Run Gillespie algorithm writing at record_times to file_name with offset time_offset");
}
