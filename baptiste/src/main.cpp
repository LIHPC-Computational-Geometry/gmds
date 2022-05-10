#include <iostream>
#include "gmds/baptiste/RLBlockSet.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/igalgo/VolFracComputation.h"
#include <gmds/io/VTKReader.h>

using namespace gmds;

Mesh readMesh(std::string filename)
{
	Mesh mesh = Mesh(MeshModel(DIM2|F|N|F2N));
	gmds::IGMeshIOService ioService(&mesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(filename);
	return mesh;
}

void testVolFrac()
{
	Mesh AImprintMesh = readMesh("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk");
	RLBlockSet blockSet = RLBlockSet(MeshModel(DIM2|F|N|F2N|N2F));
	blockSet.setFromFile("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk", 8, 4);
	Mesh AMesh = blockSet.m_mesh;
	std::cout << "The files have been read successfully" << "\n";

	for(auto faceID: blockSet.m_mesh.faces())
	{
		std::cout << faceID << "\n";
		if (faceID > 25)
		{
			blockSet.deleteBlock(faceID);
			std::cout << "Face " << faceID << " deleted" << "\n";
		}
	}

	AMesh.newVariable<double, GMDS_FACE>("volFrac");
	volfraccomputation_2d(&AMesh, &AImprintMesh, AMesh.getVariable<double, GMDS_FACE>("volFrac"));
	std::cout << "The function has been called successfully" << "\n";
	Variable<double>* res = AMesh.getVariable<double, GMDS_FACE>("volFrac");
	std::cout << res->value(0) << "\n";
	for (int i=0; i<32; ++i)
	{
		std::cout << "IoU n°" << i << " : " << res->value(i) << "\n";
	}

	for(auto faceID: blockSet.m_mesh.faces())
	{
		blockSet.editCorner(faceID, false, "y", -3);
		blockSet.editCorner(faceID, true, "x", 2);
	}

	volfraccomputation_2d(&AMesh, &AImprintMesh, AMesh.getVariable<double, GMDS_FACE>("volFrac"));

	blockSet.saveMesh("MyBlockSet");
}

int main()
{
	std::cout << "Hello World" << "\n";

	// RLBlockSet blockSet(MeshModel(DIM3|F|N|F2N|N2F));
	// blockSet.setFromFile("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk");

	testVolFrac();

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

	// blockSet.saveMesh("MyBlockSet");
	return 0;
}