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
void set_up_geom_model(gmds::cad::FACManager* AGeomModel, const std::string AFileName)
{
	gmds::Mesh vol_mesh(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::F | gmds::E | gmds::N | gmds::R2N | gmds::R2F | gmds::R2E | gmds::F2N | gmds::F2R | gmds::F2E
	                                    | gmds::E2F | gmds::E2N | gmds::N2E));
	std::string vtk_file = AFileName;
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


TEST(ClassificationTestSuite, test_chord_collapse)
{
	gmds::cad::FACManager geom_model;
	set_up_geom_model(&geom_model,"tet_in_box.vtk");
	gmds::blocking::CurvedBlocking bl(&geom_model, false);

	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::R|gmds::R2N));

	gmds::GridBuilder gb(&m,3);
	gb.execute(4,1.0, 4, 1.0, 4, 1.0);
	bl.init_from_mesh(m);


	std::vector<gmds::blocking::CurvedBlocking::Face> bl_faces = bl.get_all_faces();

	ASSERT_EQ(bl.get_nb_cells<3>(),27);

	bool found_face = false;
	auto face_id = -1;
	gmds::math::Point seed(1.5,1.5,0.0);
	for(auto i=0; i<bl_faces.size() && !found_face;i++){
		gmds::blocking::CurvedBlocking::Face fi = bl_faces[i];
		gmds::math::Point ci = bl.get_center_of_face(fi);
		if(ci.distance2(seed)<0.1){
			found_face = true;
			face_id=i;
		}
	}

	std::vector<gmds::blocking::CurvedBlocking::Node> f_nodes = bl.get_nodes_of_face(bl_faces[face_id]);
	bl.collapse_chord(bl_faces[face_id],f_nodes[0],f_nodes[2]);

	ASSERT_EQ(bl.get_nb_cells<3>(),24);

	gmds::Mesh m_out(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E | gmds::F | gmds::R | gmds::E2N | gmds::F2N | gmds::R2N));
	bl.convert_to_mesh(m_out);

	gmds::IGMeshIOService ioService(&m_out);
	gmds::VTKWriter writer(&ioService);
	writer.setCellOptions(gmds::N | gmds::F);
	writer.setDataOptions(gmds::N |gmds::F );
	writer.write("collapse.vtk");
}