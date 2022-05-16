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

namespace gmds {
class RLBlockSet
{
 public:

	RLBlockSet(MeshModel model);

	RLBlockSet();

	//RLBlockSet(RLBlockSet& original);

	~RLBlockSet();

	gmds::Mesh m_mesh;
	double xSize;
	double ySize;

	std::vector<int> getAllFaces(){
		std::vector<int> v;
		for (int faceID : m_mesh.faces())
		{
			v.push_back(faceID);
		}
		return v;
	};

	void setFrame(double xMin, double yMin, double xMax, double yMax, int nX = 3, int nY = 3);

	int countBlocks();

	void deleteBlock(const int faceID);

	void editCorner(const int faceID, bool v, std::string axis, int range);

	void saveMesh(std::string title);

	void setFromFile(std::string filename, int nX = 3, int nY = 3);
};
}
#endif     // GMDS_RL_BLOCK_SET_H
