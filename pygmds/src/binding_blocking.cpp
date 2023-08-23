/*----------------------------------------------------------------------------*/
#include "gmds/math/Point.h"
#include "gmds/blocking/CurvedBlocking.h"
#include "gmds/blocking/CurvedBlockingClassifier.h"
/*----------------------------------------------------------------------------*/
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/*----------------------------------------------------------------------------*/
namespace py = pybind11;
/*----------------------------------------------------------------------------*/
void bind_blocking(py::module &m){

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
            .def("get_normal_of_face", &gmds::blocking::CurvedBlocking::get_normal_of_face)
            .def("cut_sheet",
                 static_cast<void (gmds::blocking::CurvedBlocking::*)(const gmds::blocking::CurvedBlocking::Edge)>(&gmds::blocking::CurvedBlocking::cut_sheet))
            .def("cut_sheet_with_point", static_cast<void (gmds::blocking::CurvedBlocking::*)(const gmds::blocking::CurvedBlocking::Edge, const gmds::math::Point &)>(
                    &gmds::blocking::CurvedBlocking::cut_sheet))
            .def("cut_sheet_with_param", static_cast<void (gmds::blocking::CurvedBlocking::*)(const gmds::TCellID, const double)>(
                    &gmds::blocking::CurvedBlocking::cut_sheet))
            .def("check_capt_element",&gmds::blocking::CurvedBlocking::check_capt_element)
            .def("capt_element",&gmds::blocking::CurvedBlocking::capt_element)
            .def("get_all_sheet_edge_sets", &gmds::blocking::CurvedBlocking::get_all_sheet_edge_sets)
            .def("get_projection_info", &gmds::blocking::CurvedBlocking::get_projection_info)
            .def("create_block", static_cast<gmds::blocking::CurvedBlocking::Block (gmds::blocking::CurvedBlocking::*)(
                    gmds::math::Point &, gmds::math::Point &, gmds::math::Point &, gmds::math::Point &, gmds::math::Point &, gmds::math::Point &,
                    gmds::math::Point &, gmds::math::Point &)>(&gmds::blocking::CurvedBlocking::create_block))
            .def("remove_block", static_cast<void (gmds::blocking::CurvedBlocking::*)(gmds::blocking::CurvedBlocking::Block)>(
                    &gmds::blocking::CurvedBlocking::remove_block))
            .def("remove_block_with_id", static_cast<void (gmds::blocking::CurvedBlocking::*)(const gmds::TCellID)>(
                    &gmds::blocking::CurvedBlocking::remove_block))
            .def("info", &gmds::blocking::CurvedBlocking::info)
            .def("convert_to_mesh", &gmds::blocking::CurvedBlocking::convert_to_mesh)
            .def("save_vtk_blocking", &gmds::blocking::CurvedBlocking::save_vtk_blocking)
            .def("get_block_id", &gmds::blocking::CurvedBlocking::get_block_id);

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
            .def("classify", &gmds::blocking::CurvedBlockingClassifier::classify)
            .def("detect_classification_errors", &gmds::blocking::CurvedBlockingClassifier::detect_classification_errors)
            .def("actions_cut_current_state", &gmds::blocking::CurvedBlockingClassifier::list_Possible_Cuts);
}