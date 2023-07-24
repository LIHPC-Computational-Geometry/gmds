/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/blocking/CurvedBlocking.h>
#include <gmds/blocking/CurvedBlockingClassifier.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <unit_test_config.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
/*----------------------------------------------------------------------------*/
/**@brief setup function that initialize a geometric model using the faceted
 * representation and an input vtk file name. The vtk file must contain a
 * tetrahedral mesh
 *
 * @param AGeomModel geometric model we initialize
 * @param AFileName vtk filename
 */
void set_up_file(gmds::cad::FACManager* AGeomModel, const std::string AFileName)
{
	gmds::Mesh vol_mesh(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::F | gmds::E | gmds::N | gmds::R2N | gmds::R2F | gmds::R2E | gmds::F2N | gmds::F2R | gmds::F2E
													| gmds::E2F | gmds::E2N | gmds::N2E));
	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir +"/"+ AFileName;
	gmds::IGMeshIOService ioService(&vol_mesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::R);
	vtkReader.read(vtk_file);
	gmds::MeshDoctor doc(&vol_mesh);
	doc.buildFacesAndR2F();
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
	AGeomModel->initFrom3DMesh(&vol_mesh);

}

TEST(ExecutionActionsTestSuite,cube){

	gmds::cad::FACManager geom_model;
	set_up_file(&geom_model,"Cube.vtk");
	gmds::blocking::CurvedBlocking bl(&geom_model,true);
	gmds::blocking::CurvedBlockingClassifier classifier(&bl);


	classifier.clear_classification();


	auto errors = classifier.classify();


	//Check nb points of the geometry and nb nodes of the blocking
	ASSERT_EQ(8,geom_model.getNbPoints());
	ASSERT_EQ(12,geom_model.getNbCurves());
	ASSERT_EQ(6,geom_model.getNbSurfaces());
	ASSERT_EQ(8,bl.get_all_nodes().size());
	ASSERT_EQ(12,bl.get_all_edges().size());
	ASSERT_EQ(6,bl.get_all_faces().size());

	//Do 1 cut
	auto e = bl.get_all_edges()[0];

	auto e2 = bl.gmap()->attribute<1>(bl.gmap()->alpha<1>(e->dart()));
	bl.cut_sheet(e);
	auto coloredFaces = classifier.blocking_color_faces();
	//Check nb points of the geometry and nb nodes of the blocking after the split
	ASSERT_EQ(8,geom_model.getNbPoints());
	ASSERT_EQ(12,geom_model.getNbCurves());
	ASSERT_EQ(6,geom_model.getNbSurfaces());
	ASSERT_EQ(12,bl.get_all_nodes().size());
	ASSERT_EQ(20,bl.get_all_edges().size());
	ASSERT_EQ(11,bl.get_all_faces().size());



	errors = classifier.classify();
	//Check nb nodes/edges/faces no classified after the first split
	ASSERT_EQ(0,errors.non_classified_nodes.size());
	ASSERT_EQ(0,errors.non_classified_edges.size());
	ASSERT_EQ(1,errors.non_classified_faces.size());

	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(0,errors.non_captured_points.size());
	ASSERT_EQ(0,errors.non_captured_curves.size());
	//ASSERT_EQ(0,errors.non_captured_surfaces.size());

	bl.cut_sheet(e2);
	//Check nb points of the geometry and nb nodes of the blocking after the split
	ASSERT_EQ(8,geom_model.getNbPoints());
	ASSERT_EQ(12,geom_model.getNbCurves());
	ASSERT_EQ(6,geom_model.getNbSurfaces());
	ASSERT_EQ(18,bl.get_all_nodes().size());
	ASSERT_EQ(33,bl.get_all_edges().size());
	ASSERT_EQ(20,bl.get_all_faces().size());

	errors = classifier.classify();


	//Check nb nodes/edges/faces no classified after the second split
	ASSERT_EQ(0,errors.non_classified_nodes.size());
	ASSERT_EQ(1,errors.non_classified_edges.size());
	ASSERT_EQ(4,errors.non_classified_faces.size());



	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::E|gmds::F|gmds::R|gmds::E2N|gmds::F2N|gmds::R2N));
	bl.convert_to_mesh(m);


	gmds::IGMeshIOService ios(&m);
	gmds::VTKWriter vtk_writer(&ios);
	vtk_writer.setCellOptions(gmds::N|gmds::R);
	vtk_writer.setDataOptions(gmds::N|gmds::R);
	vtk_writer.write("debug_blocking.vtk");
	gmds::VTKWriter vtk_writer_edges(&ios);
	vtk_writer_edges.setCellOptions(gmds::N|gmds::E);
	vtk_writer_edges.setDataOptions(gmds::N|gmds::E);
	vtk_writer_edges.write("debug_blocking_edges.vtk");
}

TEST(ExecutionActionsTestSuite,cb1){

	gmds::cad::FACManager geom_model;
	set_up_file(&geom_model,"cb1.vtk");
	gmds::blocking::CurvedBlocking bl(&geom_model,true);
	gmds::blocking::CurvedBlockingClassifier classifier(&bl);


	classifier.clear_classification();


	auto errors = classifier.classify();


	//Check nb points of the geometry and nb nodes of the blocking
	ASSERT_EQ(12,geom_model.getNbPoints());
	ASSERT_EQ(18,geom_model.getNbCurves());
	ASSERT_EQ(8,geom_model.getNbSurfaces());
	ASSERT_EQ(8,bl.get_all_nodes().size());
	ASSERT_EQ(12,bl.get_all_edges().size());
	ASSERT_EQ(6,bl.get_all_faces().size());


	//Check elements class and captured
	//Check nb nodes/edges/faces no classified
	ASSERT_EQ(2,errors.non_classified_nodes.size());
	ASSERT_EQ(5,errors.non_classified_edges.size());
	ASSERT_EQ(6,errors.non_classified_faces.size());

	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(6,errors.non_captured_points.size());
	ASSERT_EQ(11,errors.non_captured_curves.size());
	ASSERT_EQ(8,errors.non_captured_surfaces.size());


	//Do 1 cut
	auto e = bl.get_all_edges()[1];

	auto e2 = bl.gmap()->attribute<1>(bl.gmap()->alpha<1>(e->dart()));
	bl.cut_sheet(e);
	//Check nb points of the geometry and nb nodes of the blocking after the split
	ASSERT_EQ(12,geom_model.getNbPoints());
	ASSERT_EQ(18,geom_model.getNbCurves());
	ASSERT_EQ(8,geom_model.getNbSurfaces());
	ASSERT_EQ(12,bl.get_all_nodes().size());
	ASSERT_EQ(20,bl.get_all_edges().size());
	ASSERT_EQ(11,bl.get_all_faces().size());



	errors = classifier.classify();
	//Check nb nodes/edges/faces no classified after the first split
	ASSERT_EQ(2,errors.non_classified_nodes.size());
	ASSERT_EQ(8,errors.non_classified_edges.size());
	ASSERT_EQ(11,errors.non_classified_faces.size());

	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(4,errors.non_captured_points.size());
	ASSERT_EQ(8,errors.non_captured_curves.size());
	ASSERT_EQ(8,errors.non_captured_surfaces.size());


	//Second Split
	bl.cut_sheet(e2);


	//Check nb points of the geometry and nb nodes of the blocking after the split
	ASSERT_EQ(12,geom_model.getNbPoints());
	ASSERT_EQ(18,geom_model.getNbCurves());
	ASSERT_EQ(8,geom_model.getNbSurfaces());
	ASSERT_EQ(18,bl.get_all_nodes().size());
	ASSERT_EQ(33,bl.get_all_edges().size());
	ASSERT_EQ(20,bl.get_all_faces().size());






	errors = classifier.classify();


	//Check nb nodes/edges/faces no classified after the second split
	ASSERT_EQ(2,errors.non_classified_nodes.size()); //sans delete 2, avec delete 0
	ASSERT_EQ(5,errors.non_classified_edges.size()); //sans delete 5, avec delete 0
	ASSERT_EQ(8,errors.non_classified_faces.size()); //sans delete 8, avec delete 2

	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(0,errors.non_captured_points.size());
	ASSERT_EQ(0,errors.non_captured_curves.size());
	ASSERT_EQ(2,errors.non_captured_surfaces.size());



	// To delete the useless block
	auto eS = bl.get_all_edges()[5];
	auto ed = bl.get_edges_of_node(bl.get_all_nodes()[7]);
	auto fe = bl.get_faces_of_edge(ed[1]);

	bl.remove_block(bl.get_blocks_of_face(fe[0])[0]);



	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::E|gmds::F|gmds::R|gmds::E2N|gmds::F2N|gmds::R2N));
	bl.convert_to_mesh(m);


	gmds::IGMeshIOService ios(&m);
	gmds::VTKWriter vtk_writer(&ios);
	vtk_writer.setCellOptions(gmds::N|gmds::R);
	vtk_writer.setDataOptions(gmds::N|gmds::R);
	vtk_writer.write("debug_blocking_cb1.vtk");
	gmds::VTKWriter vtk_writer_edges(&ios);
	vtk_writer_edges.setCellOptions(gmds::N|gmds::E);
	vtk_writer_edges.setDataOptions(gmds::N|gmds::E);
	vtk_writer_edges.write("debug_blocking_edges_cb1.vtk");
}



TEST(ExecutionActionsTestSuite,cb2){

	gmds::cad::FACManager geom_model;
	set_up_file(&geom_model,"cb2.vtk");
	gmds::blocking::CurvedBlocking bl(&geom_model,true);
	gmds::blocking::CurvedBlockingClassifier classifier(&bl);


	classifier.clear_classification();

	auto errors = classifier.classify();


	//Check nb points of the geometry and nb nodes of the blocking
	ASSERT_EQ(16,geom_model.getNbPoints());
	ASSERT_EQ(24,geom_model.getNbCurves());
	ASSERT_EQ(10,geom_model.getNbSurfaces());
	ASSERT_EQ(8,bl.get_all_nodes().size());
	ASSERT_EQ(12,bl.get_all_edges().size());
	ASSERT_EQ(6,bl.get_all_faces().size());



	//Check elements class and captured
	//Check nb nodes/edges/faces no classified
	ASSERT_EQ(0,errors.non_classified_nodes.size());
	ASSERT_EQ(0,errors.non_classified_edges.size());
	ASSERT_EQ(6,errors.non_classified_faces.size());

	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(8,errors.non_captured_points.size());
	ASSERT_EQ(12,errors.non_captured_curves.size());
	ASSERT_EQ(10,errors.non_captured_surfaces.size());

	auto listEdgesPara = bl.get_all_sheet_edge_sets();
	auto paramCut = classifier.get_cut_info(errors.non_captured_points[0],listEdgesPara);
	//Do 1 cut
	bl.cut_sheet(paramCut.first,paramCut.second);


	//============================================
	listEdgesPara = bl.get_all_sheet_edge_sets();
	paramCut = classifier.get_cut_info(errors.non_captured_points[1],listEdgesPara);

	//Do 2 cut
	bl.cut_sheet(paramCut.first,paramCut.second);

	//============================================
	listEdgesPara = bl.get_all_sheet_edge_sets();
	paramCut = classifier.get_cut_info(errors.non_captured_points[2],listEdgesPara);

	//Do 3 cut
	bl.cut_sheet(paramCut.first,paramCut.second);



	//============================================
	listEdgesPara = bl.get_all_sheet_edge_sets();
	paramCut = classifier.get_cut_info(errors.non_captured_points[6],listEdgesPara);

	//Do 4 cut
	bl.cut_sheet(paramCut.first,paramCut.second);



	//Try class after spliting

	errors = classifier.classify();


	//Check nb points of the geometry and nb nodes of the blocking
	ASSERT_EQ(16,geom_model.getNbPoints());
	ASSERT_EQ(24,geom_model.getNbCurves());
	ASSERT_EQ(10,geom_model.getNbSurfaces());
	ASSERT_EQ(32,bl.get_all_nodes().size());
	ASSERT_EQ(64,bl.get_all_edges().size());
	ASSERT_EQ(42,bl.get_all_faces().size());


	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(0,errors.non_captured_points.size());
	ASSERT_EQ(0,errors.non_captured_curves.size());
	ASSERT_EQ(4,errors.non_captured_surfaces.size());




	auto listB = bl.get_all_blocks();

	/*
	for(int i =0; i<listB.size();i++){
		std::cout<<"Block : "<<listB[i]->info().topo_id<<" a la position : "<<i<<std::endl;
	}*/


	//Delete the useless block
	bl.remove_block(listB[0]);


	//Try class after spliting
	errors = classifier.classify();


	//Check nb points of the geometry and nb nodes of the blocking
	ASSERT_EQ(16,geom_model.getNbPoints());
	ASSERT_EQ(24,geom_model.getNbCurves());
	ASSERT_EQ(10,geom_model.getNbSurfaces());
	ASSERT_EQ(32,bl.get_all_nodes().size());
	ASSERT_EQ(64,bl.get_all_edges().size());
	ASSERT_EQ(40,bl.get_all_faces().size());


	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(0,errors.non_captured_points.size());
	ASSERT_EQ(0,errors.non_captured_curves.size());
	ASSERT_EQ(0,errors.non_captured_surfaces.size());



	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::E|gmds::F|gmds::R|gmds::E2N|gmds::F2N|gmds::R2N));
	bl.convert_to_mesh(m);


	gmds::IGMeshIOService ios(&m);
	gmds::VTKWriter vtk_writer(&ios);
	vtk_writer.setCellOptions(gmds::N|gmds::R);
	vtk_writer.setDataOptions(gmds::N|gmds::R);
	vtk_writer.write("debug_blocking_cb2.vtk");
	gmds::VTKWriter vtk_writer_edges(&ios);
	vtk_writer_edges.setCellOptions(gmds::N|gmds::E);
	vtk_writer_edges.setDataOptions(gmds::N|gmds::E);
	vtk_writer_edges.write("debug_blocking_edges_cb2.vtk");
	gmds::VTKWriter vtk_writer_faces(&ios);
	vtk_writer_faces.setCellOptions(gmds::N|gmds::F);
	vtk_writer_faces.setDataOptions(gmds::N|gmds::F);
	vtk_writer_faces.write("debug_blocking_faces_cb2.vtk");
}




TEST(ExecutionActionsTestSuite,cb3){

	gmds::cad::FACManager geom_model;
	set_up_file(&geom_model,"cb3.vtk");
	gmds::blocking::CurvedBlocking bl(&geom_model,true);
	gmds::blocking::CurvedBlockingClassifier classifier(&bl);


	classifier.clear_classification();

	auto errors = classifier.classify();


	//Check nb points of the geometry and nb nodes of the blocking
	ASSERT_EQ(16,geom_model.getNbPoints());
	ASSERT_EQ(24,geom_model.getNbCurves());
	ASSERT_EQ(10,geom_model.getNbSurfaces());
	ASSERT_EQ(8,bl.get_all_nodes().size());
	ASSERT_EQ(12,bl.get_all_edges().size());
	ASSERT_EQ(6,bl.get_all_faces().size());


	//Check elements class and captured
	//Check nb nodes/edges/faces no classified
	ASSERT_EQ(0,errors.non_classified_nodes.size());
	ASSERT_EQ(2,errors.non_classified_edges.size());
	ASSERT_EQ(6,errors.non_classified_faces.size());

	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(8,errors.non_captured_points.size());
	ASSERT_EQ(14,errors.non_captured_curves.size());
	ASSERT_EQ(10,errors.non_captured_surfaces.size());

	//===================================================================

	auto listEdgesPara = bl.get_all_sheet_edge_sets();
	auto paramCut = classifier.get_cut_info(errors.non_captured_points[3],listEdgesPara);
	//Do 1 cut
	bl.cut_sheet(paramCut.first,paramCut.second);


	errors = classifier.classify();
	//Check elements class and captured
	//Check nb nodes/edges/faces no classified
	ASSERT_EQ(0,errors.non_classified_nodes.size());
	ASSERT_EQ(5,errors.non_classified_edges.size());
	ASSERT_EQ(11,errors.non_classified_faces.size());

	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(6,errors.non_captured_points.size());
	ASSERT_EQ(11,errors.non_captured_curves.size());
	ASSERT_EQ(10,errors.non_captured_surfaces.size());

	//========================================================================
	listEdgesPara = bl.get_all_sheet_edge_sets();
	paramCut = classifier.get_cut_info(errors.non_captured_points[3],listEdgesPara);
	//Do 2 cut
	bl.cut_sheet(paramCut.first,paramCut.second);


	errors = classifier.classify();
	//Check elements class and captured
	//Check nb nodes/edges/faces no classified
	ASSERT_EQ(0,errors.non_classified_nodes.size());
	ASSERT_EQ(8,errors.non_classified_edges.size());
	ASSERT_EQ(16,errors.non_classified_faces.size());

	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(4,errors.non_captured_points.size());
	ASSERT_EQ(8,errors.non_captured_curves.size());
	ASSERT_EQ(10,errors.non_captured_surfaces.size());

	//========================================================================
	listEdgesPara = bl.get_all_sheet_edge_sets();
	paramCut = classifier.get_cut_info(errors.non_captured_points[0],listEdgesPara);
	//Do 2 cut
	bl.cut_sheet(paramCut.first,paramCut.second);


	errors = classifier.classify();
	//Check elements class and captured
	//Check nb nodes/edges/faces no classified
	ASSERT_EQ(0,errors.non_classified_nodes.size());
	ASSERT_EQ(2,errors.non_classified_edges.size());
	ASSERT_EQ(10,errors.non_classified_faces.size());

	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(0,errors.non_captured_points.size());
	ASSERT_EQ(0,errors.non_captured_curves.size());
	ASSERT_EQ(3,errors.non_captured_surfaces.size());


	//======================================================
	//==================== DELETE ==============================

	auto listB = bl.get_all_blocks();
	//Delete the useless block
	bl.remove_block(listB[1]);

	errors = classifier.classify();
	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(0,errors.non_captured_points.size());
	ASSERT_EQ(0,errors.non_captured_curves.size());
	ASSERT_EQ(0,errors.non_captured_surfaces.size());



	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::E|gmds::F|gmds::R|gmds::E2N|gmds::F2N|gmds::R2N));
	bl.convert_to_mesh(m);


	gmds::IGMeshIOService ios(&m);
	gmds::VTKWriter vtk_writer(&ios);
	vtk_writer.setCellOptions(gmds::N|gmds::R);
	vtk_writer.setDataOptions(gmds::N|gmds::R);
	vtk_writer.write("debug_blocking_cb3.vtk");
	gmds::VTKWriter vtk_writer_edges(&ios);
	vtk_writer_edges.setCellOptions(gmds::N|gmds::E);
	vtk_writer_edges.setDataOptions(gmds::N|gmds::E);
	vtk_writer_edges.write("debug_blocking_edges_cb3.vtk");
}


TEST(ExecutionActionsTestSuite,cb4){

	gmds::cad::FACManager geom_model;
	set_up_file(&geom_model,"cb4.vtk");
	gmds::blocking::CurvedBlocking bl(&geom_model,true);
	gmds::blocking::CurvedBlockingClassifier classifier(&bl);


	classifier.clear_classification();

	auto errors = classifier.classify();


	//Check nb points of the geometry and nb nodes of the blocking
	ASSERT_EQ(14,geom_model.getNbPoints());
	ASSERT_EQ(21,geom_model.getNbCurves());
	ASSERT_EQ(9,geom_model.getNbSurfaces());

	ASSERT_EQ(8,bl.get_all_nodes().size());
	ASSERT_EQ(12,bl.get_all_edges().size());
	ASSERT_EQ(6,bl.get_all_faces().size());


	//Check elements class and captured
	//Check nb nodes/edges/faces no classified
	ASSERT_EQ(1,errors.non_classified_nodes.size());
	ASSERT_EQ(3,errors.non_classified_edges.size());
	ASSERT_EQ(6,errors.non_classified_faces.size());

	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(7,errors.non_captured_points.size());
	ASSERT_EQ(12,errors.non_captured_curves.size());
	ASSERT_EQ(9,errors.non_captured_surfaces.size());



	//===================================================================

	auto listEdgesPara = bl.get_all_sheet_edge_sets();
	auto paramCut = classifier.get_cut_info(errors.non_captured_points[0],listEdgesPara);
	//Do 1 cut
	bl.cut_sheet(paramCut.first,paramCut.second);


	errors = classifier.classify();
	//Check elements class and captured
	//Check nb nodes/edges/faces no classified
	ASSERT_EQ(1,errors.non_classified_nodes.size());
	ASSERT_EQ(7,errors.non_classified_edges.size());
	ASSERT_EQ(11,errors.non_classified_faces.size());

	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(6,errors.non_captured_points.size());
	ASSERT_EQ(11,errors.non_captured_curves.size());
	ASSERT_EQ(9,errors.non_captured_surfaces.size());


	//===================================================================
	listEdgesPara = bl.get_all_sheet_edge_sets();
	paramCut = classifier.get_cut_info(errors.non_captured_points[0],listEdgesPara);
	//Do 2 cut
	bl.cut_sheet(paramCut.first,paramCut.second);


	errors = classifier.classify();
	//Check elements class and captured
	//Check nb nodes/edges/faces no classified
	ASSERT_EQ(1,errors.non_classified_nodes.size());
	ASSERT_EQ(14,errors.non_classified_edges.size());
	ASSERT_EQ(20,errors.non_classified_faces.size());

	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(4,errors.non_captured_points.size());
	ASSERT_EQ(8,errors.non_captured_curves.size());
	ASSERT_EQ(9,errors.non_captured_surfaces.size());

	//===================================================================
	listEdgesPara = bl.get_all_sheet_edge_sets();
	paramCut = classifier.get_cut_info(errors.non_captured_points[0],listEdgesPara);
	//Do 3 cut
	bl.cut_sheet(paramCut.first,paramCut.second);


	errors = classifier.classify();
	//Check elements class and captured
	//Check nb nodes/edges/faces no classified
	ASSERT_EQ(1,errors.non_classified_nodes.size());
	ASSERT_EQ(6,errors.non_classified_edges.size());
	ASSERT_EQ(15,errors.non_classified_faces.size());

	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(0,errors.non_captured_points.size());
	ASSERT_EQ(0,errors.non_captured_curves.size());
	ASSERT_EQ(3,errors.non_captured_surfaces.size());




	//======================================================
	//==================== DELETE ==============================

	auto listB = bl.get_all_blocks();
	//Delete the useless block
	bl.remove_block(listB[4]);

	errors = classifier.classify();
	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(0,errors.non_captured_points.size());
	ASSERT_EQ(0,errors.non_captured_curves.size());
	ASSERT_EQ(0,errors.non_captured_surfaces.size());



	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::E|gmds::F|gmds::R|gmds::E2N|gmds::F2N|gmds::R2N));
	bl.convert_to_mesh(m);


	gmds::IGMeshIOService ios(&m);
	gmds::VTKWriter vtk_writer(&ios);
	vtk_writer.setCellOptions(gmds::N|gmds::R);
	vtk_writer.setDataOptions(gmds::N|gmds::R);
	vtk_writer.write("debug_blocking.vtk");
	gmds::VTKWriter vtk_writer_edges(&ios);
	vtk_writer_edges.setCellOptions(gmds::N|gmds::E);
	vtk_writer_edges.setDataOptions(gmds::N|gmds::E);
	vtk_writer_edges.write("debug_blocking_edges.vtk");
}

TEST(ExecutionActionsTestSuite,cb5){

	gmds::cad::FACManager geom_model;
	set_up_file(&geom_model,"cb5.vtk");
	gmds::blocking::CurvedBlocking bl(&geom_model,true);
	gmds::blocking::CurvedBlockingClassifier classifier(&bl);


	classifier.clear_classification();

	auto errors = classifier.classify();


	//Check nb points of the geometry and nb nodes of the blocking
	ASSERT_EQ(28,geom_model.getNbPoints());
	ASSERT_EQ(42,geom_model.getNbCurves());
	ASSERT_EQ(16,geom_model.getNbSurfaces());

	ASSERT_EQ(8,bl.get_all_nodes().size());
	ASSERT_EQ(12,bl.get_all_edges().size());
	ASSERT_EQ(6,bl.get_all_faces().size());


	//Check elements class and captured
	//Check nb nodes/edges/faces no classified
	ASSERT_EQ(2,errors.non_classified_nodes.size());
	ASSERT_EQ(8,errors.non_classified_edges.size());
	ASSERT_EQ(6,errors.non_classified_faces.size());

	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(22,errors.non_captured_points.size());
	ASSERT_EQ(38,errors.non_captured_curves.size());
	ASSERT_EQ(16,errors.non_captured_surfaces.size());


	//======================= START CUTTING =======================
	//===================================================================

	auto listEdgesPara = bl.get_all_sheet_edge_sets();
	auto paramCut = classifier.get_cut_info(errors.non_captured_points[0],listEdgesPara);
	//Do 1 cut
	bl.cut_sheet(paramCut.first,paramCut.second);


	errors = classifier.classify();
	//Check elements class and captured
	//Check nb nodes/edges/faces no classified
	ASSERT_EQ(2,errors.non_classified_nodes.size());
	ASSERT_EQ(13,errors.non_classified_edges.size());
	ASSERT_EQ(11,errors.non_classified_faces.size());

	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(20,errors.non_captured_points.size());
	ASSERT_EQ(35,errors.non_captured_curves.size());
	ASSERT_EQ(16,errors.non_captured_surfaces.size());


	//=======================================================================
	listEdgesPara = bl.get_all_sheet_edge_sets();
	paramCut = classifier.get_cut_info(errors.non_captured_points[3],listEdgesPara);
	//Do 1 cut
	bl.cut_sheet(paramCut.first,paramCut.second);


	errors = classifier.classify();
	//=======================================================================
	listEdgesPara = bl.get_all_sheet_edge_sets();
	paramCut = classifier.get_cut_info(errors.non_captured_points[3],listEdgesPara);
	//Do 1 cut
	bl.cut_sheet(paramCut.first,paramCut.second);


	errors = classifier.classify();

	//=======================================================================
	listEdgesPara = bl.get_all_sheet_edge_sets();
	paramCut = classifier.get_cut_info(errors.non_captured_points[5],listEdgesPara);
	//Do 1 cut
	bl.cut_sheet(paramCut.first,paramCut.second);


	errors = classifier.classify();

	//=======================================================================
	listEdgesPara = bl.get_all_sheet_edge_sets();
	paramCut = classifier.get_cut_info(errors.non_captured_points[1],listEdgesPara);
	//Do 1 cut
	bl.cut_sheet(paramCut.first,paramCut.second);


	errors = classifier.classify();

	//=======================================================================
	listEdgesPara = bl.get_all_sheet_edge_sets();
	paramCut = classifier.get_cut_info(errors.non_captured_points[0],listEdgesPara);
	//Do 1 cut
	bl.cut_sheet(paramCut.first,paramCut.second);


	errors = classifier.classify();

	//=======================================================================
	listEdgesPara = bl.get_all_sheet_edge_sets();
	paramCut = classifier.get_cut_info(errors.non_captured_points[1],listEdgesPara);
	//Do 1 cut
	bl.cut_sheet(paramCut.first,paramCut.second);


	errors = classifier.classify();

	srand (time(NULL));

	//===================================================================
	while(errors.non_captured_points.size()!=0) {

		int sizeVect = errors.non_captured_points.size()-1;
		int random = rand() % sizeVect;
		listEdgesPara = bl.get_all_sheet_edge_sets();
		paramCut = classifier.get_cut_info(errors.non_captured_points[random], listEdgesPara);
		// Do 1 cut
		bl.cut_sheet(paramCut.first, paramCut.second);

		errors = classifier.classify();

		std::cout<<"CUT"<<std::endl;
	}
	/*
	//Check elements class and captured
	//Check nb nodes/edges/faces no classified
	ASSERT_EQ(27,errors.non_classified_nodes.size());
	ASSERT_EQ(311,errors.non_classified_edges.size());
	ASSERT_EQ(312,errors.non_classified_faces.size());

	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(0,errors.non_captured_points.size());
	ASSERT_EQ(21,errors.non_captured_curves.size());
	ASSERT_EQ(23,errors.non_captured_surfaces.size());


	auto listB = bl.get_all_blocks();

	//Delete the useless block
	bl.remove_block(listB[82]);
	bl.remove_block(listB[75]);
	bl.remove_block(listB[66]);
	bl.remove_block(listB[64]);
	bl.remove_block(listB[63]);
	bl.remove_block(listB[49]);
	bl.remove_block(listB[48]);
	bl.remove_block(listB[41]);
	bl.remove_block(listB[40]);

	bl.remove_block(listB[73]);
	bl.remove_block(listB[67]);
	bl.remove_block(listB[29]);
	bl.remove_block(listB[28]);
	bl.remove_block(listB[26]);
	bl.remove_block(listB[25]);
	bl.remove_block(listB[22]);
	bl.remove_block(listB[19]);
	bl.remove_block(listB[15]);
	bl.remove_block(listB[12]);
	bl.remove_block(listB[8]);


	bl.remove_block(listB[46]);
	bl.remove_block(listB[14]);
	bl.remove_block(listB[13]);
	bl.remove_block(listB[6]);

	bl.remove_block(listB[44]);
	 */
	/*
	for(int i =0; i<listB.size();i++){
		std::cout<<"Block : "<<listB[i]->info().topo_id<<" a la position : "<<i<<std::endl;
	}*/
/*
	errors = classifier.classify();

	//Check elements class and captured
	//Check nb nodes/edges/faces no classified
	ASSERT_EQ(27,errors.non_classified_nodes.size());
	ASSERT_EQ(311,errors.non_classified_edges.size());
	ASSERT_EQ(312,errors.non_classified_faces.size());

	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(0,errors.non_captured_points.size());
	ASSERT_EQ(21,errors.non_captured_curves.size());
	ASSERT_EQ(23,errors.non_captured_surfaces.size());
*/

	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::E|gmds::F|gmds::R|gmds::E2N|gmds::F2N|gmds::R2N));
	bl.convert_to_mesh(m);


	gmds::IGMeshIOService ios(&m);
	gmds::VTKWriter vtk_writer(&ios);
	vtk_writer.setCellOptions(gmds::N|gmds::R);
	vtk_writer.setDataOptions(gmds::N|gmds::R);
	vtk_writer.write("debug_blocking.vtk");
	gmds::VTKWriter vtk_writer_edges(&ios);
	vtk_writer_edges.setCellOptions(gmds::N|gmds::E);
	vtk_writer_edges.setDataOptions(gmds::N|gmds::E);
	vtk_writer_edges.write("debug_blocking_edges.vtk");
}

