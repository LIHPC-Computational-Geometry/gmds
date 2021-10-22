/*------------------------------------------------------------------------*/
// Created by Claire Roche on 21/10/2021.
/*------------------------------------------------------------------------*/
#include <gmds/claire/Smooth2D.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
	std::cout << "=== CLAIRE ALGO ====" << std::endl;
	Mesh m(MeshModel(DIM3 | R | F | N | F2N | N2F));

	//==================================================================
	// MESH READING
	//==================================================================
	std::cout << "Reading " << std::endl;
	std::string fIn, fOut;

	if (argc != 3)
		throw gmds::GMDSException("[Wrong parameters] usage should have two vtk files (in-out)");

	fIn = std::string(argv[1]);
	fOut = std::string(argv[2]);

	IGMeshIOService ioService(&m);
	VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::R);
	vtkReader.read(fIn);

	//==================================================================
	// MESH PREPARATION
	//==================================================================
	MeshDoctor doctor(&m);
	doctor.buildFacesAndR2F();
	doctor.buildEdgesAndX2E();
	doctor.updateUpwardConnectivity();

}