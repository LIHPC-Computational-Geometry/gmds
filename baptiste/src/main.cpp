#include <iostream>
#include "gmds/baptiste/RLBlockSet.h"
#include "gmds/ig/Mesh.h"
#include "gmds/io/IGMeshIOService.h"


//using namespace std;
using namespace gmds;

int main()
{
	std::cout<<"Hello World" << "\n";

	RLBlockSet blockSet(MeshModel(DIM3|F|N|F2N|N2F));

	double xMin = 0.2;
	double yMin = 0.5;
	double xMax = 5.3;
	double yMax = 4.7;

	blockSet.setFrame(xMin, yMin, xMax, yMax);

	int nbFaces = blockSet.getNumberOfBlocks();
	std::cout<<nbFaces;

	return 0;
}