#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "gmds/baptiste/RLBlockSet.h"
#include "gmds/baptiste/tools.h"
#include "gmds/baptiste/capsule.h"

namespace py = pybind11;

PYBIND11_MODULE(environment, m)
{
	py::class_<gmds::RLBlockSet>(m, "RLBlockSet")
	   .def(py::init<gmds::MeshModel>())
	   .def(py::init<>())
	   //.def("__copy__",  [](const gmds::RLBlockSet &self) { return gmds::RLBlockSet(self); })
	   //.def("__deepcopy__", [](const gmds::RLBlockSet &self, py::dict) { return gmds::RLBlockSet(self); })
	   .def_readwrite("m_mesh", &gmds::RLBlockSet::m_mesh)
	   .def_readwrite("xSize", &gmds::RLBlockSet::xSize)
	   .def_readwrite("ySize", &gmds::RLBlockSet::ySize)
		.def("setFrame", &gmds::RLBlockSet::setFrame)
		.def("countBlocks", &gmds::RLBlockSet::countBlocks)
		.def("deleteBlock", &gmds::RLBlockSet::deleteBlock)
		.def("editCorner", &gmds::RLBlockSet::editCorner)
	   .def("saveMesh", &gmds::RLBlockSet::saveMesh)
	   .def("setFromFile", &gmds::RLBlockSet::setFromFile)
	   .def("getAllFaces", &gmds::RLBlockSet::getAllFaces)
	   .def("getReward", &gmds::RLBlockSet::getReward)
	   .def("isValid", &gmds::RLBlockSet::isValid)
	   .def("getStateID", &gmds::RLBlockSet::getStateID)
	   .def("overlap", &gmds::RLBlockSet::overlap)
       .def("getLocalIou", &gmds::RLBlockSet::getLocalIou)
       .def("getMinMaxCoordinates", &gmds::RLBlockSet::getMinMaxCoordinates);

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
		.def("getRegion", &gmds::Mesh::get<gmds::Region>)
	   .def("newDouble", &gmds::Mesh::newVariable<double, gmds::GMDS_FACE>)
	   .def("getDouble", &gmds::Mesh::getVariable<double, gmds::GMDS_FACE>);

	py::class_<gmds::FaceContainer>(m, "FaceContainer")
	   .def("__iter__", [](gmds::FaceContainer &s) { return py::make_iterator(s.begin(), s.end()); },
	      py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */);

	py::class_<gmds::Face>(m, "Face")
	   .def("area", &gmds::Face::area)
	   .def("nbNodes", &gmds::Face::nbNodes);

	m.def("cloneBlockSet", &gmds::cloneBlockSet);

	py::class_<gmds::Capsule>(m, "Capsule")
	   .def(py::init<>())
	   .def_readwrite("m_mesh", &gmds::Capsule::m_mesh)
	   .def("readMesh", &gmds::Capsule::readMesh);
 }