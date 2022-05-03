#include <iostream>
#include "gmds/baptiste/RLBlockSet.h"
#include "gmds/io/IGMeshIOService.h"

using namespace gmds;

int main()
{
	std::cout << "Hello World" << "\n";
	RLBlockSet blockSet(MeshModel(DIM3|F|N|F2N|N2F));
	blockSet.setFromFile("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk");

	/*
	double xMin = 0;
	double yMin = 0;
	double xMax = 9;
	double yMax = 9;
	blockSet.setFrame(xMin, yMin, xMax, yMax);

	int nbFaces = blockSet.countBlocks();
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

	nbFaces = blockSet.countBlocks();
	std::cout << "Number of blocks : " << nbFaces << "\n";

	for(auto faceID: blockSet.m_mesh.faces())
	{
		blockSet.editCorner(faceID, false, "y", -3);
	}
	*/

	blockSet.saveMesh("MyBlockSet");
	return 0;
}