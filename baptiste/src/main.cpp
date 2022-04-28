#include <iostream>
#include "gmds/baptiste/RLBlockSet.h"
#include "gmds/ig/Mesh.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKWriter.h"

//using namespace std;
using namespace gmds;

int main()
{
	std::cout << "Hello World" << "\n";

	RLBlockSet blockSet(MeshModel(DIM3|F|N|F2N|N2F));

	double xMin = 0;
	double yMin = 0;
	double xMax = 9;
	double yMax = 9;
	blockSet.setFrame(xMin, yMin, xMax, yMax);

	int nbFaces = blockSet.getNumberOfBlocks();
	std::cout << "Number of blocks : " << nbFaces << "\n";

	for(auto faceID: blockSet.m_mesh.faces())
	{
		std::cout << faceID << "\n";
		if (faceID > 4)
		{
			blockSet.deleteBlock(faceID);
			std::cout << "Face " << faceID << " deleted" << "\n";
		}
	}

	nbFaces = blockSet.getNumberOfBlocks();
	std::cout << "Number of blocks : " << nbFaces << "\n";

	for(auto faceID: blockSet.m_mesh.faces())
	{
		blockSet.editCorner(faceID, false, "y", -2);
	}

	// Save Mesh
	IGMeshIOService ioService(&blockSet.m_mesh);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("MyBlockSet.vtk");
	
	return 0;
}