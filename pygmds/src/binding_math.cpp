/*----------------------------------------------------------------------------*/
#include "gmds/math/Point.h"
#include "gmds/math/Vector.h"
/*----------------------------------------------------------------------------*/
#include <pybind11/pybind11.h>
/*----------------------------------------------------------------------------*/
namespace py = pybind11;
/*----------------------------------------------------------------------------*/
void bind_math(py::module &m){
	py::class_<gmds::math::Vector3d>(m, "Vector3d")
	   .def(py::init<gmds::TCoord &, gmds::TCoord &, gmds::TCoord &>())
	   .def("x", &gmds::math::Vector3d::X)
	   .def("y", &gmds::math::Vector3d::Y)
	   .def("z", &gmds::math::Vector3d::Z)
	   .def("norm", &gmds::math::Vector3d::norm)
	   .def("get_normalize", &gmds::math::Vector3d::getNormalize)
	   .def("dot", &gmds::math::Vector3d::dot)
	   .def("cross", &gmds::math::Vector3d::cross);

	py::class_<gmds::math::Point>(m, "Point")
	   .def(py::init<gmds::TCoord &, gmds::TCoord &, gmds::TCoord &>())
	   .def("distance", &gmds::math::Point::distance)
	   .def("x", static_cast<gmds::TCoord &(gmds::math::Point::*) ()>(&gmds::math::Point::X))
	   .def("y", static_cast<gmds::TCoord &(gmds::math::Point::*) ()>(&gmds::math::Point::Y))
	   .def("z", static_cast<gmds::TCoord &(gmds::math::Point::*) ()>(&gmds::math::Point::Z));
}