/*----------------------------------------------------------------------------*/
#include "gmds/blocking/CurvedBlocking.h"
#include "gmds/blocking/CurvedBlockingClassifier.h"
#include "gmds/cadfac/FACManager.h"
#include "gmds/cadfac/FACPoint.h"
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

	py::class_<gmds::math::Point>(m, "Point")
	   .def(py::init<gmds::TCoord &, gmds::TCoord &, gmds::TCoord &>())
	   .def("x", static_cast<gmds::TCoord &(gmds::math::Point::*)()>(&gmds::math::Point::X))
	   .def("y", static_cast<gmds::TCoord &(gmds::math::Point::*)()>(&gmds::math::Point::Y))
	   .def("z", static_cast<gmds::TCoord &(gmds::math::Point::*)()>(&gmds::math::Point::Z));

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



	py::enum_<gmds::cad::GeomCurve::CurvatureInfo>(m, "CurvatureInfo")
	   .value("CURVATURE_FLAT", gmds::cad::GeomCurve::Flat)
	   .value("CURVATURE_CONVEX", gmds::cad::GeomCurve::Convex)
	   .value("CURVATURE_CONCAVE", gmds::cad::GeomCurve::Concave)
	   .value("CURVATURE_SMOOTH_CONVEX", gmds::cad::GeomCurve::SmoothConvex)
	   .value("CURVATURE_SMOOTH_CONCAV", gmds::cad::GeomCurve::SmoothConcave)
	   .export_values();

	py::class_<gmds::cad::GeomManager>(m, "GeomManager");

		py::class_<gmds::cad::FACManager, gmds::cad::GeomManager>(m, "FACManager")
			.def(py::init<>())
			.def("get_nb_volumes", &gmds::cad::FACManager::getNbVolumes)
			.def("get_nb_surfaces", &gmds::cad::FACManager::getNbSurfaces)
			.def("get_nb_curves", &gmds::cad::FACManager::getNbCurves)
			.def("get_nb_points", &gmds::cad::FACManager::getNbPoints)
			.def("get_volumes", static_cast< std::vector<gmds::cad::GeomVolume*> (gmds::cad::FACManager::*) ()const>(&gmds::cad::FACManager::getVolumes), py::return_value_policy::reference)
			.def("get_surfaces", static_cast< std::vector<gmds::cad::GeomSurface*> (gmds::cad::FACManager::*) ()const>(&gmds::cad::FACManager::getSurfaces),py::return_value_policy::reference)
			.def("get_curves", static_cast< std::vector<gmds::cad::GeomCurve*> (gmds::cad::FACManager::*)() const>(&gmds::cad::FACManager::getCurves), py::return_value_policy::reference)
			.def("get_points", static_cast< std::vector<gmds::cad::GeomPoint*> (gmds::cad::FACManager::*)() const>(&gmds::cad::FACManager::getPoints), py::return_value_policy::reference)
			.def("get_volume", &gmds::cad::FACManager::getVolume)
			.def("get_surface", &gmds::cad::FACManager::getSurface)
			.def("get_curve", &gmds::cad::FACManager::getCurve)
			.def("get_point", &gmds::cad::FACManager::getPoint)
			.def("init_from_3d_mesh", &gmds::cad::FACManager::initFrom3DMesh);

	   py::class_<gmds::cad::GeomPoint>(m, "GeomPoint");
	   py::class_<gmds::cad::FACPoint, gmds::cad::GeomPoint>(m, "FACPoint")
	      .def("x", &gmds::cad::FACPoint::X)
	      .def("y", &gmds::cad::FACPoint::Y)
	      .def("z", &gmds::cad::FACPoint::Z)
	      .def("id", &gmds::cad::FACPoint::id)
	      .def("bbox", &gmds::cad::FACPoint::BBox)
	      .def("curves", &gmds::cad::FACPoint::curves, py::return_value_policy::reference)
	      .def("surfaces", &gmds::cad::FACPoint::surfaces, py::return_value_policy::reference)
	      .def("volumes", &gmds::cad::FACPoint::volumes, py::return_value_policy::reference);


	   py::class_<gmds::cad::GeomCurve>(m, "GeomCurve");
	   py::class_<gmds::cad::FACCurve, gmds::cad::GeomCurve>(m, "FACCurve")
	      .def("id", &gmds::cad::FACCurve::id)
	      .def("length", &gmds::cad::FACCurve::length)
	      .def("closest_point", &gmds::cad::FACCurve::closestPoint)
	      .def("bbox", &gmds::cad::FACCurve::BBox)
	      .def("is_loop", &gmds::cad::FACCurve::isALoop)
	      .def("curvature_info", &gmds::cad::FACCurve::getCurvatureInfo)
	      .def("dihedral_angle", &gmds::cad::FACCurve::computeDihedralAngle)
	      .def("points", &gmds::cad::FACCurve::points, py::return_value_policy::reference)
	      .def("surfaces", &gmds::cad::FACCurve::surfaces, py::return_value_policy::reference)
	      .def("volumes", &gmds::cad::FACCurve::volumes, py::return_value_policy::reference);

	   py::class_<gmds::cad::GeomSurface>(m, "GeomSurface");
	   py::class_<gmds::cad::FACSurface, gmds::cad::GeomSurface>(m, "FACSurface")
	      .def("id", &gmds::cad::FACSurface::id)
	      .def("area", &gmds::cad::FACSurface::computeArea)
	      .def("center", &gmds::cad::FACSurface::center)
	      .def("closest_point", &gmds::cad::FACSurface::closestPoint)
	      .def("bbox", &gmds::cad::FACSurface::BBox)
	      .def("points", &gmds::cad::FACSurface::points, py::return_value_policy::reference)
	      .def("curves", &gmds::cad::FACSurface::curves, py::return_value_policy::reference)
	      .def("volumes", &gmds::cad::FACSurface::volumes, py::return_value_policy::reference);

	   py::class_<gmds::cad::GeomVolume>(m, "GeomVolume");
	   py::class_<gmds::cad::FACVolume, gmds::cad::GeomVolume>(m, "FACVolume")
	      .def("id", &gmds::cad::FACVolume::id)
	      .def("bbox", &gmds::cad::FACVolume::BBox)
	      .def("points", &gmds::cad::FACVolume::points, py::return_value_policy::reference)
	      .def("curves", &gmds::cad::FACVolume::curves, py::return_value_policy::reference)
	      .def("surfaces", &gmds::cad::FACVolume::surfaces, py::return_value_policy::reference);


	         py::class_<gmds::blocking::CurvedBlocking>(m, "Blocking")
	   .def(py::init<gmds::cad::GeomManager *, bool>())
	   .def("get_node_info", &gmds::blocking::CurvedBlocking::get_node_info)
	   .def("get_edge_info", &gmds::blocking::CurvedBlocking::get_edge_info)
	   .def("get_face_info", &gmds::blocking::CurvedBlocking::get_face_info)
	   .def("get_block_info", &gmds::blocking::CurvedBlocking::get_block_info)
	   .def("get_all_id_nodes", &gmds::blocking::CurvedBlocking::get_all_id_nodes)
	   .def("get_all_id_edges", &gmds::blocking::CurvedBlocking::get_all_id_edges)
	   .def("get_all_id_faces", &gmds::blocking::CurvedBlocking::get_all_id_faces)
	   .def("get_all_id_blocks", &gmds::blocking::CurvedBlocking::get_all_id_blocks)
	   .def("get_nb_blocks", &gmds::blocking::CurvedBlocking::get_nb_cells<3>)
	   .def("get_nb_faces", &gmds::blocking::CurvedBlocking::get_nb_cells<2>)
	   .def("get_nb_edges", &gmds::blocking::CurvedBlocking::get_nb_cells<1>)
	   .def("get_nb_nodes", &gmds::blocking::CurvedBlocking::get_nb_cells<0>)
	   .def("get_all_blocks", &gmds::blocking::CurvedBlocking::get_all_blocks)
	   .def("get_all_faces", &gmds::blocking::CurvedBlocking::get_all_faces)
	   .def("get_all_edges", &gmds::blocking::CurvedBlocking::get_all_edges)
	   .def("get_all_nodes", &gmds::blocking::CurvedBlocking::get_all_nodes)
	   .def("get_edges_of_node", &gmds::blocking::CurvedBlocking::get_edges_of_node)
	   .def("get_faces_of_node", &gmds::blocking::CurvedBlocking::get_faces_of_node)
	   .def("get_blocks_of_node", &gmds::blocking::CurvedBlocking::get_blocks_of_node)
	   .def("get_nodes_of_edge", &gmds::blocking::CurvedBlocking::get_nodes_of_edge)
	   .def("get_faces_of_edge", &gmds::blocking::CurvedBlocking::get_faces_of_edge)
	   .def("get_blocks_of_edge", &gmds::blocking::CurvedBlocking::get_blocks_of_edge)
	   .def("get_nodes_of_face", &gmds::blocking::CurvedBlocking::get_nodes_of_face)
	   .def("get_edges_of_face", &gmds::blocking::CurvedBlocking::get_edges_of_face)
	   .def("get_blocks_of_face", &gmds::blocking::CurvedBlocking::get_blocks_of_face)
	   .def("get_edges_of_block", &gmds::blocking::CurvedBlocking::get_edges_of_block)
	   .def("get_faces_of_block", &gmds::blocking::CurvedBlocking::get_faces_of_block)
	   .def("get_nodes_of_block", &gmds::blocking::CurvedBlocking::get_nodes_of_block)
	   .def("cut_sheet",
	        static_cast<void (gmds::blocking::CurvedBlocking::*)(const gmds::blocking::CurvedBlocking::Edge)>(&gmds::blocking::CurvedBlocking::cut_sheet))
	   .def("cut_sheet_with_point", static_cast<void (gmds::blocking::CurvedBlocking::*)(const gmds::blocking::CurvedBlocking::Edge, const gmds::math::Point &)>(
	                                   &gmds::blocking::CurvedBlocking::cut_sheet))
	   .def("cut_sheet_with_param", static_cast<void (gmds::blocking::CurvedBlocking::*)(const gmds::blocking::CurvedBlocking::Edge, const double)>(
	                                   &gmds::blocking::CurvedBlocking::cut_sheet))
	   .def("get_all_sheet_edge_sets",&gmds::blocking::CurvedBlocking::get_all_sheet_edge_sets)
	   .def("get_projection_info",&gmds::blocking::CurvedBlocking::get_projection_info)
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