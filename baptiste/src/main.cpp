#include <iostream>
#include "gmds/baptiste/RLBlockSet.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/igalgo/VolFracComputation.h"
#include <gmds/io/VTKReader.h>

using namespace gmds;

Mesh readMesh(std::string filename)
{
	Mesh mesh = Mesh(MeshModel(DIM3|F|N));
	gmds::IGMeshIOService ioService(&mesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(filename);
	return mesh;
}

void testVolFrac()
{
	Mesh m1 = readMesh("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk");
	Mesh m2 = readMesh("/home/bonyb/Documents/GitHub/gmds/test_samples/A.vtk");

	m1.newVariable<double, GMDS_FACE>("volFrac");
	volfraccomputation_2d(&m1, &m2, m1.getVariable<double, GMDS_FACE>("volFrac"));
	std::cout << m1.getVariable<double, GMDS_FACE>("volFrac") << "\n";
	// Process finished with exit code 139 (interrupted by signal 11: SIGSEGV)
}

int main()
{
	std::cout << "Hello World" << "\n";
	RLBlockSet blockSet(MeshModel(DIM3|F|N|F2N|N2F));
	blockSet.setFromFile("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk");

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

	blockSet.saveMesh("MyBlockSet");
	return 0;
}