/*------------------------------------------------------------------------*/
#include "gmds/ig/MeshDoctor.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKReader.h"
#include "gmds/io/VTKWriter.h"
#include <gmds/morphMesh/MorphMesh.h>
#include <gmds/morphMesh/EllipticMorph.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/
int main(int argc, char* argv[]){

	std::string param_file(argv[1]);

	gmds::Mesh m_mesh(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::F | gmds::N | gmds::N2R | gmds::R2N | gmds::R2F | gmds::F2N | gmds::F2R));
	gmds::IGMeshIOService ioService(&m_mesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::R);
	vtkReader.setDataOptions(gmds::N|gmds::R);
	vtkReader.read(param_file);



	gmds::MeshDoctor doc(&m_mesh);
	doc.buildN2R(m_mesh.getModel());
	doc.buildFacesAndR2F();
	doc.updateUpwardConnectivity();

	//gmds::morphmesh::MorphMesh morph(&m_mesh,points,2.5);
	//morph.execute();

	morphmesh::EllipticMorph emorph(&m_mesh);

	std::vector<std::vector<double>> ellipses;
	ellipses.push_back({4,1,1,1});
	ellipses.push_back({5,1,1.2,1});
	ellipses.push_back({6,1,1,1});
	//ellipses.push_back({0.2,1,1,1});
	//ellipses.push_back({3,1.2,1.5,1});
	//ellipses.push_back({4,1.2,1.5,1});
	//ellipses.push_back({5,1.4,1.75,1.5});
	//ellipses.push_back({7,1.4,1.75,1.5});
	//ellipses.push_back({8,1,1,1});
	//ellipses.push_back({5.6,1.2,1.5,1});
	//ellipses.push_back({8.5,1.2,1.5,1});

	emorph.execute(ellipses);

	std::cout<<"test"<<std::endl;

	gmds::VTKWriter w(&ioService);
	w.setCellOptions(gmds::N|gmds::R);
	w.setDataOptions(gmds::N|gmds::R);
	w.write("morphingtest.vtk");


}