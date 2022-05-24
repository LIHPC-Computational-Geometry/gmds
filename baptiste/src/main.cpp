#include <iostream>
#include "gmds/baptiste/RLBlockSet.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/igalgo/VolFracComputation.h"
#include <gmds/io/VTKReader.h>
#include "gmds/baptiste/tools.h"
#include <map>
#include <string>

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
	//RLBlockSet blockSet = RLBlockSet(MeshModel(DIM2|F|N|F2N|N2F));
	RLBlockSet blockSet = RLBlockSet(2);
	std::cout << "Constructor Ok" << "\n";
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

double getReward(RLBlockSet blockSet, Mesh targetShape)
{
	blockSet.m_mesh.newVariable<double, GMDS_FACE>("volFrac");
	volfraccomputation_2d(&blockSet.m_mesh, &targetShape, blockSet.m_mesh.getVariable<double, GMDS_FACE>("volFrac"));
	std::cout << "The function has been called successfully" << "\n";
	Variable<double>* res = blockSet.m_mesh.getVariable<double, GMDS_FACE>("volFrac");
	double sum = 0;
	for (int i =0; i < res->getNbValues(); i++)
	{
		sum += res->value(i);
	}
	return sum;
}

void testReward()
{
	Mesh AImprintMesh = readMesh("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk");
	RLBlockSet blockSet = RLBlockSet(2);
	blockSet.setFromFile("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk", 8, 4);
	double r = getReward(blockSet, AImprintMesh);
	std::cout << r << "\n";
}


Mesh copyMesh(RLBlockSet blockSet)
{
	Tools t = Tools();
	std::cout << "tools created";
	blockSet.saveMesh("originalMesh");
	std::cout << "Mesh saved";
	t.readMesh("originalMesh.vtk");
	std::cout << "Mesh read";
	std::cout << "Mesh copied";
	return t.m_mesh;
}

Mesh cloneMesh(const Mesh &mesh)
{
	if (mesh.getDim() == 3)
	{
		throw GMDSException("Dimension must be 2");
	}
	Mesh newMesh = Mesh(MeshModel(DIM3|F|N|F2N));
	for (int nodeID:mesh.nodes())
	{
		newMesh.newNode(mesh.get<Node>(nodeID).point());
	}
	/*
	for (int faceID:mesh.faces())
	{
		newMesh.newFace(mesh.get<Face>(faceID).get<Node>());
	}
	 */
	std::set<TCellID> invalidIndexes;
	for (int faceID = 0; faceID <= mesh.getMaxLocalID(2); faceID++)
	{
		if (mesh.has<Face>(faceID))
		{
			newMesh.newFace(mesh.get<Face>(faceID).get<Node>());
		}
		else
		{
			newMesh.newFace({0, 0, 0, 0});
			invalidIndexes.insert(faceID);
		}
	}
	for (int faceID:invalidIndexes)
	{
		newMesh.deleteFace(faceID);
	}
	std::cout << mesh.getNbNodes() << "\n";
	std::cout << mesh.getNbFaces() << "\n";
	std::cout << newMesh.getNbNodes() << "\n";
	std::cout << newMesh.getNbFaces() << "\n";

	for (int faceID:mesh.faces())
	{
		std::cout << faceID;
	}
	for (int faceID:newMesh.faces())
	{
		std::cout << faceID;
	}
	return newMesh;
}

std::vector<std::map<int , int>> getAllActions()
{
	std::vector<std::map<int , int>> actions;
	for (int type = 0; type <= 1; type++)
	{
		for (int v = 0; v <= 1; v++)
		{
			for (int axis = 0; axis <= 1; axis++)
			{
				std::vector<int> ranges = {-2, -1, 1, 2};
				for (int range:ranges)
				{

					std::map<int, int> action = {
					   { 0, type },
					   { 1, v },
					   { 2, axis },
					   { 3, range}
					};

					actions.push_back(action);
				}
			}
		}
	}
	return actions;
}

void executeAction(RLBlockSet blockSet, std::map<int , int> action, int faceID)
{
	if (action[0] == 0)
	{
		blockSet.deleteBlock(faceID);
	}
	else
	{
		bool v = !!action[1];
		std::string axis;
		if (action[2]==0)
		{
			axis = "x";
		}
		else
		{
			axis = "y";
		}
		int range = action[3];
		blockSet.editCorner(faceID, v, axis, range);
	}
}


void virtualExpert(RLBlockSet blockSet, Mesh targetShape, int nMax)
{
	std::vector<std::map<int , int>> actions = getAllActions();
	for (int step=0; step<nMax; step++)
	{
		std::cout << "Step n°" << step << "\n";
		for(int faceID: blockSet.m_mesh.faces())
		{
			std::cout << "In faceIDs loop n°" << step << "\n";
			std::map<int , int> bestAction = actions[0];
			double maxReward = 0;
			for (std::map<int , int> action:actions)
			{
				std::cout << "In actions loop n°" << step << "\n";
				Mesh meshBis = copyMesh(blockSet);
				std::cout << "Mesh copied cc";
				RLBlockSet blockSetBis = RLBlockSet(2);
				blockSetBis.m_mesh = meshBis;
				executeAction(blockSetBis, action, faceID);
				double reward = getReward(blockSetBis, targetShape);
				if (reward > maxReward)
				{
					maxReward = reward;
					bestAction= action;
				}
				std::cout << "End of n°" << step << "\n";
			}
			executeAction(blockSet, bestAction, faceID);
		}
	}
}

void testDeepCopy()
{
	RLBlockSet blockSet = RLBlockSet(2);
	double xMin = 0;
	double yMin = 0;
	double xMax = 9;
	double yMax = 9;
	blockSet.setFrame(xMin, yMin, xMax, yMax);

	int nbFaces = blockSet.countBlocks();
	std::cout << "Number of blocks : " << nbFaces << "\n";

	//RLBlockSet blockSet2 = RLBlockSet(blockSet.m_mesh);
	RLBlockSet blockSet2 = RLBlockSet(2);
	Mesh mesh2 = cloneMesh(blockSet.m_mesh);

	blockSet2.m_mesh = mesh2;
	std::cout << "New block set created" << "\n";

	nbFaces = blockSet2.countBlocks();
	std::cout << "Number of blocks 2 : " << nbFaces << "\n";

	blockSet.deleteBlock(3);

	nbFaces = blockSet.countBlocks();
	std::cout << "Number of blocks : " << nbFaces << "\n";

	nbFaces = blockSet2.countBlocks();
	std::cout << "Number of blocks 2 : " << nbFaces << "\n";
}





int main()
{
	std::cout << "Hello World" << "\n";
	Mesh targetShape = readMesh("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk");
	RLBlockSet blockSet = RLBlockSet(2);
	blockSet.setFromFile("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk", 8, 4);
	//virtualExpert(blockSet, targetShape, 25);
	//blockSet.saveMesh("virtualExpert");


	//testVolFrac();
	//testReward();
	testDeepCopy();

	/*
	RLBlockSet blockSet(MeshModel(DIM3|F|N|F2N|N2F));
	// blockSet.setFromFile("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk");


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


	// blockSet.saveMesh("MyBlockSet");
	 */
	return 0;
}