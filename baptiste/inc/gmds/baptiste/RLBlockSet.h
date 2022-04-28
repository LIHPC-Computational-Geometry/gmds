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

	virtual ~RLBlockSet();

	gmds::Mesh m_mesh;

	std::vector<double> LinearSpacedArray(double a, double b, std::size_t N);

	void setFrame(double xMin, double yMin, double xMax, double yMax, int nX = 3, int nY = 3);

	int getNumberOfBlocks();

	void deleteBlock(const int faceID);

	void editCorner(const int faceID, bool v, std::string axis, int range);

};
}
#endif     // GMDS_RL_BLOCK_SET_H
