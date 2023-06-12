/*------------------------------------------------------------------------*/
#include "gmds/ig/MeshDoctor.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKReader.h"
#include "gmds/io/VTKWriter.h"
#include <gmds/morphMesh/MorphMesh.h>
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

	gmds::math::Point point(6,0.5,0);
	std::vector<gmds::math::Point> points;
	points.push_back(point);

	gmds::morphmesh::MorphMesh morph(&m_mesh,points,2.5);
	morph.execute();

	gmds::VTKWriter w(&ioService);
	w.setCellOptions(gmds::N|gmds::R);
	w.setDataOptions(gmds::N|gmds::R);
	w.write("morphingtest.vtk");


}