//
// Created by Baptiste Bony on 22/04/22.
//

#ifndef GMDS_RL_BLOCK_SET_H
#define GMDS_RL_BLOCK_SET_H
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_BAPTISTE_export.h"
#include <gmds/ig/Mesh.h>
#include "GMDSIgAlgo_export.h"
/*----------------------------------------------------------------------------*/

namespace gmds
{
class RLBlockSet
{
 public:

	RLBlockSet(MeshModel model);

	RLBlockSet();

	~RLBlockSet();

	gmds::Mesh m_mesh;
	double xSize;
	double ySize;

	void setFrame(double xMin, double yMin, double xMax, double yMax, int nX = 3, int nY = 3);

	int countBlocks();

	void deleteBlock(const int faceID);

	void editCorner(const int faceID, bool v, std::string axis, int range);

	Node findCorner(std::vector<Node> nodes, bool v);

	Node findSecondCorner(Node corner, std::vector<Node> nodes, std::string axis);

	void moveCorners(Node corner, Node secondCorner, std::string axis, int range);

	void editCornerBis(const int faceID, bool v, std::string axis, int range);

	void saveMesh(std::string title);

	void setFromFile(std::string filename, int nX = 3, int nY = 3);

	std::vector<int> getAllFaces();

	double getReward(Mesh &targetMesh);

	bool isValid();

	std::string getStateID();

	double overlap();

    double getOverlap();
};

}
#endif     // GMDS_RL_BLOCK_SET_H
