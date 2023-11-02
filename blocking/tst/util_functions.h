/*----------------------------------------------------------------------------*/
#ifndef GMDS_UTIL_FUNCTIONS_H
#define GMDS_UTIL_FUNCTIONS_H
/*----------------------------------------------------------------------------*/
#include <gmds/blocking/CurvedBlocking.h>
#include <gmds/blocking/CurvedBlockingClassifier.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
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
std::tuple<int,int,int,int> get_node_statistics(gmds::blocking::CurvedBlocking& ABlocking){
	auto nb_on_vertex=0;
	auto nb_on_curve=0;
	auto nb_on_surface=0;
	auto nb_in_volume=0;
	std::vector<gmds::blocking::CurvedBlocking::Node> all_nodes = ABlocking.get_all_nodes();
	for(auto n:all_nodes){
		if(n->info().geom_dim==0)
			nb_on_vertex++;
		else if(n->info().geom_dim==1)
			nb_on_curve++;
		else if(n->info().geom_dim==2)
			nb_on_surface++;
		else if(n->info().geom_dim==3)
			nb_in_volume++;
	}
	return std::make_tuple(nb_on_vertex,nb_on_curve,nb_on_surface,nb_in_volume);
}
/*----------------------------------------------------------------------------*/
std::tuple<int,int,int> get_edge_statistics(gmds::blocking::CurvedBlocking& ABlocking){
	auto nb_on_curve=0;
	auto nb_on_surface=0;
	auto nb_in_volume=0;
	std::vector<gmds::blocking::CurvedBlocking::Edge> all_edges = ABlocking.get_all_edges();
	for(auto e:all_edges){
		if(e->info().geom_dim==1)
			nb_on_curve++;
		else if(e->info().geom_dim==2)
			nb_on_surface++;
		else if(e->info().geom_dim==3)
			nb_in_volume++;
	}
	return std::make_tuple(nb_on_curve,nb_on_surface,nb_in_volume);
}
/*----------------------------------------------------------------------------*/
std::tuple<int,int> get_face_statistics(gmds::blocking::CurvedBlocking& ABlocking){
	auto nb_on_surface=0;
	auto nb_in_volume=0;
	std::vector<gmds::blocking::CurvedBlocking::Face> all_faces = ABlocking.get_all_faces();
	for(auto f:all_faces){
		if(f->info().geom_dim==2)
			nb_on_surface++;
		else if(f->info().geom_dim==3)
			nb_in_volume++;
	}
	return std::make_tuple(nb_on_surface,nb_in_volume);
}
/*----------------------------------------------------------------------------*/
void export_vtk(gmds::blocking::CurvedBlocking& ABlocking, int AModel, const std::string& AFileName){
	gmds::Mesh m_out(gmds::MeshModel(gmds::DIM3 | gmds::N | gmds::E | gmds::F | gmds::R | gmds::E2N | gmds::F2N | gmds::R2N));
	ABlocking.convert_to_mesh(m_out);
	gmds::IGMeshIOService ioService(&m_out);
	gmds::VTKWriter writer(&ioService);
	writer.setCellOptions(AModel);
	writer.setDataOptions(AModel);
	writer.write(AFileName);
}
/*----------------------------------------------------------------------------*/
#endif     // GMDS_UTIL_FUNCTIONS_H
/*----------------------------------------------------------------------------*/
