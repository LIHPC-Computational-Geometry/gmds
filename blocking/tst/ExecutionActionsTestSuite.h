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
	auto coloredFaces = classifier.blocking_color_faces();
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

	/*
	// For delete the useless block
	auto eS = bl.get_all_edges()[5];
	auto ed = bl.get_edges_of_node(bl.get_all_nodes()[7]);
	auto fe = bl.get_faces_of_edge(ed[1]);

	bl.remove_block(bl.get_blocks_of_face(fe[0])[0]);
	 */

	std::cout<<"AFAZFZAEFZAFZAFaz"<<std::endl;

	errors = classifier.classify();


	//Check nb nodes/edges/faces no classified after the second split
	ASSERT_EQ(2,errors.non_classified_nodes.size()); //sans delete 2, avec delete 0
	ASSERT_EQ(5,errors.non_classified_edges.size()); //sans delete 5, avec delete 0
	ASSERT_EQ(8,errors.non_classified_faces.size()); //sans delete 8, avec delete 2

	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(0,errors.non_captured_points.size());
	ASSERT_EQ(0,errors.non_captured_curves.size());
	ASSERT_EQ(2,errors.non_captured_surfaces.size());



// For delete the useless block
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
	vtk_writer.write("debug_blocking.vtk");
	gmds::VTKWriter vtk_writer_edges(&ios);
	vtk_writer_edges.setCellOptions(gmds::N|gmds::E);
	vtk_writer_edges.setDataOptions(gmds::N|gmds::E);
	vtk_writer_edges.write("debug_blocking_edges.vtk");
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
