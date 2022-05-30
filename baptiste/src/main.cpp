#include <iostream>
#include "gmds/baptiste/RLBlockSet.h"
#include "gmds/io/IGMeshIOService.h"
#include "gmds/igalgo/VolFracComputation.h"
#include "gmds/baptiste/tools.h"
#include <map>
#include <string>

using namespace gmds;

void testVolFrac()
{
	Mesh AImprintMesh = readMesh("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk");
	//RLBlockSet blockSet = RLBlockSet(MeshModel(DIM2|F|N|F2N|N2F));
	RLBlockSet blockSet = RLBlockSet();
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

void executeAction(RLBlockSet &blockSet, std::map<int , int> &action, int faceID)
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

void virtualExpert(RLBlockSet &blockSet, Mesh &targetShape, int nMax)
{
	std::vector<std::map<int , int>> actions = getAllActions();
	for (int step=0; step<nMax; step++)
	{
		std::cout << "Step n°" << step << "\n";
		for(int faceID: blockSet.m_mesh.faces())
		{
			std::cout << "In faceIDs loop with id = " << faceID << "\n";
			std::map<int , int> bestAction = actions[0];
			double maxReward = 0;
			int actionCounter = 0;
			for (auto action:actions)
			{
				if (step <= 0.5 * nMax)
				{
					if (action[0] == 1)
					{
						RLBlockSet blockSetBis = RLBlockSet();
						cloneBlockSet(blockSet, blockSetBis);
						executeAction(blockSetBis, action, faceID);
						double reward = blockSetBis.getReward(targetShape);
						if (reward > maxReward)
						{
							maxReward = reward;
							bestAction = action;
						}
					}
				}
				else
				{
					RLBlockSet blockSetBis = RLBlockSet();
					cloneBlockSet(blockSet, blockSetBis);
					executeAction(blockSetBis, action, faceID);
					double reward = blockSetBis.getReward(targetShape);
					if (reward > maxReward)
					{
						maxReward = reward;
						bestAction = action;
					}
				}
				actionCounter++;
			}
			executeAction(blockSet, bestAction, faceID);
		}
	}
}

void testDeepCopy()
{
	RLBlockSet blockSet = RLBlockSet();
	double xMin = 0;
	double yMin = 0;
	double xMax = 9;
	double yMax = 9;
	blockSet.setFrame(xMin, yMin, xMax, yMax);

	int nbFaces = blockSet.countBlocks();
	std::cout << "Number of blocks : " << nbFaces << "\n";

	RLBlockSet blockSet2 = RLBlockSet();
	cloneBlockSet(blockSet, blockSet2);
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
	//Mesh targetShape = readMesh("/home/bonyb/Documents/GitHub/gmds/test_samples/A.vtk");
	//Mesh targetShape = readMesh("/home/bonyb/Documents/GitHub/gmds/test_samples/HolesInSquare0.vtk");


	RLBlockSet blockSet = RLBlockSet();

	blockSet.setFromFile("/home/bonyb/Documents/GitHub/gmds/test_samples/Curved_Shape1_ref2.vtk", 3, 3);
	//blockSet.setFromFile("/home/bonyb/Documents/GitHub/gmds/test_samples/A.vtk", 3, 3);
	//blockSet.setFromFile("/home/bonyb/Documents/GitHub/gmds/test_samples/HolesInSquare0.vtk", 3, 3);
	//double reward = blockSet.getReward(targetShape);
	//std::cout << reward << "\n";
	virtualExpert(blockSet, targetShape, 5);

	blockSet.saveMesh("/home/bonyb/Documents/virtualExpertCurvedShape");

	//Mesh m = Mesh(MeshModel(F|N|F2N|N2F));

	/*
	RLBlockSet blockSet = RLBlockSet();
	blockSet.setFrame(0, 0, 3, 3, 3, 3);
	//Mesh m2 = cloneMesh(blockSet.m_mesh);
	//RLBlockSet blockSet2 = copyBlockSet(blockSet);


	RLBlockSet blockSetBis = RLBlockSet();
	cloneBlockSet(blockSet, blockSetBis);
	 */

	//testReward();
	//testDeepCopy();

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