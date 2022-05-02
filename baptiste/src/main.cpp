#include <iostream>
#include "gmds/baptiste/RLBlockSet.h"
#include "gmds/io/IGMeshIOService.h"


// TO DO
// Always displace two corners at the same time to keep the rectangular shape

#include <gmds/io/VTKReader.h>

//using namespace std;
using namespace gmds;

Node findMinNode(Mesh mesh)
{

	std::vector<Node> nodes;
	for (int nodeID : mesh.nodes())
	{
		Node node = mesh.get<Node>(nodeID);
		nodes.push_back(node);
	}
	std::vector<Node>::iterator iter;
	iter = std::min_element(nodes.begin(), nodes.end(),[](Node a, Node b){return a.X() <= b.X() and a.Y() <= b.Y();});
	int nodeID = iter->id();
	Node node = mesh.get<Node>(nodeID);
	return node;
}


Node finMaxNode(Mesh mesh)
{
	std::vector<Node> nodes;
	for (int nodeID : mesh.nodes())
	{
		Node node = mesh.get<Node>(nodeID);
		nodes.push_back(node);
	}
	std::vector<Node>::iterator iter;
	iter = std::max_element(nodes.begin(), nodes.end(),[](Node a, Node b){return a.X() <= b.X() and a.Y() <= b.Y();});
	int nodeID = iter->id();
	Node node = mesh.get<Node>(nodeID);
	return node;
}

int main()
{
	std::cout << "Hello World" << "\n";

	Mesh a_mesh = Mesh(MeshModel(DIM3|F|N));
	gmds::IGMeshIOService ioService(&a_mesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk");

	int nb = a_mesh.getNbFaces();
	std::cout << nb;

	RLBlockSet blockSet(MeshModel(DIM3|F|N|F2N|N2F));

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

	blockSet.saveMesh("MyBlockSet");
	return 0;
}