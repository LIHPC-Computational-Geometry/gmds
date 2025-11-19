//
// Created by rochec on 24/03/2022.
//

#include "gmds/aero/CGNSWriterND.h"
#include "gmds/ig/MeshDoctor.h"
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <iostream>

/*----------------------------------------------------------------------------*/
using namespace gmds;

/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
	std::string param_file(argv[1]);
	std::string output_file(argv[2]);

	std::cout << "=== TEST ====" << std::endl;

	Mesh *mesh = new Mesh(MeshModel(DIM3 | N | F | R | N2F | N2R | F2N | F2R | R2N | R2F));
	gmds::IGMeshIOService ioService(mesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::R);
	vtkReader.setDataOptions(gmds::N|gmds::R);
	vtkReader.read(param_file);

	std::cout<<"End reading"<<std::endl;

	MeshDoctor doc(mesh);
	doc.buildFacesAndR2F();
	doc.buildF2R(mesh->getModel());

	gmds::aero::CGNSWriterND writer(mesh, 3);
	writer.write("", output_file, "/home/calderans/");
}