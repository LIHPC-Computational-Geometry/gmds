#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include "gmds/baptiste/RLBlockSet.h"
#include <pybind11/embed.h>

namespace py = pybind11;

PYBIND11_MODULE(environment, m)
{
	py::class_<gmds::RLBlockSet>(m, "RLBlockSet")
	   .def(py::init<gmds::MeshModel>())
	   .def(py::init<>())
	   .def_readwrite("m_mesh", &gmds::RLBlockSet::m_mesh)
	   .def_readwrite("xSize", &gmds::RLBlockSet::xSize)
	   .def_readwrite("ySize", &gmds::RLBlockSet::ySize)
	   .def("LinearSpacedArray", &gmds::RLBlockSet::LinearSpacedArray)
		.def("setFrame", &gmds::RLBlockSet::setFrame)
		.def("countBlocks", &gmds::RLBlockSet::countBlocks)
		.def("deleteBlock", &gmds::RLBlockSet::deleteBlock)
		.def("editCorner", &gmds::RLBlockSet::editCorner)
	   .def("saveMesh", &gmds::RLBlockSet::saveMesh)
	   .def("setFromFile", &gmds::RLBlockSet::setFromFile)
	   .def("getAllFaces", &gmds::RLBlockSet::getAllFaces)
	   ;

	py::class_<gmds::MeshModel>(m, "MeshModel")
	   .def(py::init<const int &>());

	py::class_<gmds::Mesh>(m, "Mesh")
	   .def(py::init<const gmds::MeshModel&>())
	   .def("faces", &gmds::Mesh::faces)
		.def("getNbNodes", &gmds::Mesh::getNbNodes)
		.def("getNbEdges", &gmds::Mesh::getNbEdges)
		.def("getNbFaces", &gmds::Mesh::getNbFaces)
		.def("getNbRegions", &gmds::Mesh::getNbRegions)
		.def("getNode", &gmds::Mesh::get<gmds::Node>)
		.def("getEdge", &gmds::Mesh::get<gmds::Edge>)
		.def("getFace", &gmds::Mesh::get<gmds::Face>)
		.def("getRegion", &gmds::Mesh::get<gmds::Region>);

	py::class_<gmds::FaceContainer>(m, "FaceContainer")
	   .def("__iter__", [](gmds::FaceContainer &s) { return py::make_iterator(s.begin(), s.end()); },
	      py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */);

	py::class_<gmds::Face>(m, "Face");
 }