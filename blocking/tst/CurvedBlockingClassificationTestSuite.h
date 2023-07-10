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
void set_up(gmds::cad::FACManager* AGeomModel, const std::string AFileName)
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
/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingClassifierTestSuite, simple_box)
{
	gmds::cad::FACManager geom_model;
	set_up(&geom_model,"tet_in_box.vtk");
	gmds::blocking::CurvedBlocking bl(&geom_model,true);

	gmds::blocking::CurvedBlockingClassifier classifier(&bl);
	classifier.clear_classification();
	classifier.classify();

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

/*----------------------------------------------------------------------------*/
TEST(CurvedBlockingClassifierTestSuite, b80)
{
	gmds::cad::FACManager geom_model;
	set_up(&geom_model,"B80.vtk");
	gmds::blocking::CurvedBlocking bl(&geom_model,true);

	gmds::blocking::CurvedBlockingClassifier classifier(&bl);
	classifier.clear_classification();
	classifier.classify(4,0.5);

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

TEST(CurvedBlockingClassifierTestSuite,ReturnErrors){

	gmds::cad::FACManager geom_model;
	set_up(&geom_model,"B0.vtk");
	gmds::blocking::CurvedBlocking bl(&geom_model,true);

	geom_model.write_surfaces("debug_classify.vtk");
	gmds::blocking::CurvedBlockingClassifier classifier(&bl);




	classifier.clear_classification();

	//Check nb points of the geometry and nb nodes of the blocking
	ASSERT_EQ(12,geom_model.getNbPoints());
	ASSERT_EQ(8,bl.get_all_nodes().size());


	auto errors = classifier.classify();


	//Check nb nodes/edges/faces no classified
	ASSERT_EQ(0,errors.non_classified_nodes.size());
	ASSERT_EQ(2,errors.non_classified_edges.size());
	//ASSERT_EQ(3,errors.non_classified_faces.size());

	//errors = classifier.detect_classification_errors();
	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(4,errors.non_captured_points.size());
	ASSERT_EQ(8,errors.non_captured_curves.size());
	//ASSERT_EQ(5,errors.non_captured_surfaces.size());
	auto all_faces = bl.get_all_faces();
	classifier.exterior_faces_coloration(all_faces);

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


TEST(CurvedBlockingClassifierTestSuite,classifyFace){

	gmds::cad::FACManager geom_model;
	set_up(&geom_model,"Cube.vtk");
	gmds::blocking::CurvedBlocking bl(&geom_model,true);
	gmds::blocking::CurvedBlockingClassifier classifier(&bl);




	classifier.clear_classification();

	//Check nb points of the geometry and nb nodes of the blocking
	ASSERT_EQ(8,geom_model.getNbPoints());
	ASSERT_EQ(12,geom_model.getNbCurves());
	ASSERT_EQ(6,geom_model.getNbSurfaces());
	ASSERT_EQ(8,bl.get_all_nodes().size());
	ASSERT_EQ(12,bl.get_all_edges().size());
	ASSERT_EQ(6,bl.get_all_faces().size());


	auto errors = classifier.classify();


	//Check nb nodes/edges/faces no classified
	ASSERT_EQ(0,errors.non_classified_nodes.size());
	ASSERT_EQ(0,errors.non_classified_edges.size());
	ASSERT_EQ(0,errors.non_classified_faces.size());

	//Check nb points/curves/surfaces no captured
	ASSERT_EQ(0,errors.non_captured_points.size());
	ASSERT_EQ(0,errors.non_captured_curves.size());
	//ASSERT_EQ(0,errors.non_captured_surfaces.size());
	auto all_faces = bl.get_all_faces();
	classifier.exterior_faces_coloration(all_faces);

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
