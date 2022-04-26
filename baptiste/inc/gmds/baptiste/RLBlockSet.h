//
// Created by Baptiste Bony on 22/04/22.
//

#ifndef GMDS_RL_BLOCK_SET_H
#define GMDS_RL_BLOCK_SET_H
/*----------------------------------------------------------------------------*/
#include "gmds/ig/Mesh.h"
#include "GMDSIgAlgo_export.h"
/*----------------------------------------------------------------------------*/

namespace gmds {
class RLBlockSet
{
 public:

	gmds::Mesh m_mesh;

	RLBlockSet(Mesh *AMesh);


	
	void setFrame(double xMin, double yMin, double xMax, double yMax, int nX, int nY);

	int getNumberOfBlocks();

	void deleteBlock(const int faceID);

};
}
#endif     // GMDS_RL_BLOCK_SET_H
