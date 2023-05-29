/*----------------------------------------------------------------------------*/
#include "gmds/blocking/CurvedBlocking.h"
#include "gmds/blocking/CurvedBlockingClassifier.h"
#include "gmds/cadfac/FACManager.h"
#include "gmds/ig/Mesh.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKReader.h"
#include "gmds/utils/CommonTypes.h"
/*----------------------------------------------------------------------------*/
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/*----------------------------------------------------------------------------*/
namespace py = pybind11;
// ig is the submodule name
PYBIND11_MODULE(gmds, m)
{
	// current version
	m.attr("__version__") = "dev";
	// module description
	m.doc() = "gmds is the python binding to a set of features provided by the gmds project";

	py::enum_<gmds::ECellType>(m, "CellType")
	   .value("GMDS_NODE", gmds::ECellType::GMDS_NODE)
	   .value("GMDS_EDGE", gmds::ECellType::GMDS_EDGE)
	   .value("GMDS_FACE", gmds::ECellType::GMDS_FACE)
	   .value("GMDS_REGION", gmds::ECellType::GMDS_REGION)
	   .export_values();

	py::enum_<gmds::EMeshDefinition>(m, "ModelDefinition", py::arithmetic())
	   .value("DIM3", gmds::EMeshDefinition::DIM3)
	   .value("N", gmds::EMeshDefinition::N)
	   .value("E", gmds::EMeshDefinition::E)
	   .value("F", gmds::EMeshDefinition::F)
	   .value("R", gmds::EMeshDefinition::R)
	   .value("N2E", gmds::EMeshDefinition::N2E)
	   .value("N2F", gmds::EMeshDefinition::N2F)
	   .value("N2R", gmds::EMeshDefinition::N2R)
	   .value("E2N", gmds::EMeshDefinition::E2N)
	   .value("E2F", gmds::EMeshDefinition::E2F)
	   .value("E2R", gmds::EMeshDefinition::E2R)
	   .value("F2N", gmds::EMeshDefinition::F2N)
	   .value("F2E", gmds::EMeshDefinition::F2E)
	   .value("F2R", gmds::EMeshDefinition::F2R)
	   .value("R2N", gmds::EMeshDefinition::R2N)
	   .value("R2E", gmds::EMeshDefinition::R2E)
	   .value("R2F", gmds::EMeshDefinition::R2F)
	   .export_values();

	py::class_<gmds::MeshModel>(m, "MeshModel").def(py::init<const int &>());

	py::class_<gmds::math::Point>(m, "Point").def(py::init<gmds::TCoord &, gmds::TCoord &, gmds::TCoord &>());

	py::class_<gmds::Variable<int>>(m, "VariableInt")
	   .def(py::init<const std::string &>())
	   .def("set", &gmds::Variable<int>::set)
	   .def("value", static_cast<int &(gmds::Variable<int>::*) (const gmds::TCellID &)>(&gmds::Variable<int>::value));

	py::class_<gmds::Mesh>(m, "Mesh")
	   .def(py::init<const gmds::MeshModel &>())
	   .def("get_nb_nodes", &gmds::Mesh::getNbNodes)
	   .def("get_nb_edges", &gmds::Mesh::getNbEdges)
	   .def("get_nb_faces", &gmds::Mesh::getNbFaces)
	   .def("get_nb_regions", &gmds::Mesh::getNbRegions)
	   .def("get_node", &gmds::Mesh::get<gmds::Node>)
	   .def("get_edge", &gmds::Mesh::get<gmds::Edge>)
	   .def("get_face", &gmds::Mesh::get<gmds::Face>)
	   .def("get_region", &gmds::Mesh::get<gmds::Region>)
	   .def("new_int_variable_at_node", [](gmds::Mesh &AThis, std::string &AName) { return AThis.getOrCreateVariable<int, gmds::GMDS_NODE>("AName"); });

	py::class_<gmds::MeshDoctor>(m, "MeshDoctor")
	   .def(py::init<gmds::Mesh *>())
	   .def("build_faces_and_R2F", &gmds::MeshDoctor::buildFacesAndR2F)
	   .def("build_edges_and_X2E", &gmds::MeshDoctor::buildEdgesAndX2E)
	   .def("update_upward_connectivity", &gmds::MeshDoctor::updateUpwardConnectivity);

	py::class_<gmds::IMeshIOService>(m, "AbstractIOService");

	py::class_<gmds::IGMeshIOService, gmds::IMeshIOService>(m, "IOService").def(py::init<gmds::Mesh *>());

	py::class_<gmds::VTKReader>(m, "VTKReader")
	   .def(py::init<gmds::IMeshIOService *>())
	   .def("read", &gmds::VTKReader::read)
	   .def("set_cell_options", &gmds::VTKReader::setCellOptions)
	   .def("set_data_options", &gmds::VTKReader::setDataOptions);

	py::class_<gmds::cad::GeomManager>(m, "GeomManager");

	py::class_<gmds::cad::FACManager, gmds::cad::GeomManager>(m, "FACManager")
	   .def(py::init<>())
	   .def("get_nb_volumes", &gmds::cad::FACManager::getNbVolumes)
	   .def("get_nb_surfaces", &gmds::cad::FACManager::getNbSurfaces)
	   .def("get_nb_curves", &gmds::cad::FACManager::getNbCurves)
	   .def("get_nb_points", &gmds::cad::FACManager::getNbPoints)
	   .def("get_volumes", &gmds::cad::FACManager::getVolumes)
	   .def("get_surfaces", &gmds::cad::FACManager::getSurfaces)
	   .def("get_curves", &gmds::cad::FACManager::getCurves)
	   .def("get_points", &gmds::cad::FACManager::getPoints)
	   .def("get_volume", &gmds::cad::FACManager::getVolume)
	   .def("get_surface", &gmds::cad::FACManager::getSurface)
	   .def("get_curve", &gmds::cad::FACManager::getCurve)
	   .def("get_point", &gmds::cad::FACManager::getPoint)
	   .def("init_from_3d_mesh", &gmds::cad::FACManager::initFrom3DMesh);


	py::class_<gmds::blocking::CurvedBlocking>(m, "Blocking")
	   .def(py::init<gmds::cad::GeomManager *, bool>())
	   .def("get_nb_blocks", &gmds::blocking::CurvedBlocking::get_nb_cells<3>)
	   .def("get_nb_faces", &gmds::blocking::CurvedBlocking::get_nb_cells<2>)
	   .def("get_nb_edges", &gmds::blocking::CurvedBlocking::get_nb_cells<1>)
	   .def("get_nb_nodes", &gmds::blocking::CurvedBlocking::get_nb_cells<0>)
	   .def("get_all_blocks", &gmds::blocking::CurvedBlocking::get_all_blocks)
	   .def("get_all_faces", &gmds::blocking::CurvedBlocking::get_all_faces)
	   .def("get_all_edges", &gmds::blocking::CurvedBlocking::get_all_edges)
	   .def("get_all_nodes", &gmds::blocking::CurvedBlocking::get_all_nodes)
	   .def("get_edges_of_block", &gmds::blocking::CurvedBlocking::get_edges_of_block)
	   .def("get_nodes_of_edge", &gmds::blocking::CurvedBlocking::get_nodes_of_edge)
	   .def("get_nodes_of_face", &gmds::blocking::CurvedBlocking::get_nodes_of_face)
	   .def("get_nodes_of_block", &gmds::blocking::CurvedBlocking::get_nodes_of_block)
	   .def("get_all_sheet_edges", &gmds::blocking::CurvedBlocking::get_all_sheet_edges)
	   .def("cut_sheet",
	        static_cast<void (gmds::blocking::CurvedBlocking::*)(const gmds::blocking::CurvedBlocking::Edge)>(&gmds::blocking::CurvedBlocking::cut_sheet))
	   .def("cut_sheet_with_point", static_cast<void (gmds::blocking::CurvedBlocking::*)(const gmds::blocking::CurvedBlocking::Edge, const gmds::math::Point &)>(
	                                   &gmds::blocking::CurvedBlocking::cut_sheet))
	   .def("cut_sheet_with_param", static_cast<void (gmds::blocking::CurvedBlocking::*)(const gmds::blocking::CurvedBlocking::Edge, const double)>(
	                                   &gmds::blocking::CurvedBlocking::cut_sheet))
	   .def("create_block", static_cast<gmds::blocking::CurvedBlocking::Block (gmds::blocking::CurvedBlocking::*)(
	                           gmds::math::Point &, gmds::math::Point &, gmds::math::Point &, gmds::math::Point &, gmds::math::Point &, gmds::math::Point &,
	                           gmds::math::Point &, gmds::math::Point &)>(&gmds::blocking::CurvedBlocking::create_block))
	   .def("remove_block", &gmds::blocking::CurvedBlocking::remove_block)
	   .def("info", &gmds::blocking::CurvedBlocking::info)
	   .def("convert_to_mesh", &gmds::blocking::CurvedBlocking::convert_to_mesh);

		py::class_<gmds::blocking::CurvedBlocking::Block>(m, "Block");

	py::class_<gmds::blocking::ClassificationErrors>(m, "ClassificationErrors")
	   .def_readonly("non_captured_points", &gmds::blocking::ClassificationErrors::non_captured_points)
	   .def_readonly("non_captured_curves", &gmds::blocking::ClassificationErrors::non_captured_curves)
	   .def_readonly("non_captured_surfaces", &gmds::blocking::ClassificationErrors::non_captured_surfaces)
	   .def_readonly("non_classified_nodes", &gmds::blocking::ClassificationErrors::non_classified_nodes)
	   .def_readonly("non_classified_edges", &gmds::blocking::ClassificationErrors::non_classified_edges)
	   .def_readonly("non_classified_faces", &gmds::blocking::ClassificationErrors::non_classified_faces);

	py::class_<gmds::blocking::CurvedBlockingClassifier>(m, "BlockingClassifier")
	   .def(py::init<gmds::blocking::CurvedBlocking *>())
	   .def("clear_classification", &gmds::blocking::CurvedBlockingClassifier::clear_classification)
	   .def("classify", &gmds::blocking::CurvedBlockingClassifier::classify);


}