#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include "gmds/baptiste/RLBlockSet.h"
#include <pybind11/embed.h>

namespace py = pybind11;

PYBIND11_MODULE(environment, m) {
	py::class_<gmds::RLBlockSet>(m, "RLBlockSet")
	   .def(py::init<gmds::MeshModel>())
	   .def(py::init<>())
	   .def("LinearSpacedArray", &gmds::RLBlockSet::LinearSpacedArray);
}