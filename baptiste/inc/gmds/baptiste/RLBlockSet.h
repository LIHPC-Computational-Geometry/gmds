//
// Created by Baptiste Bony on 22/04/22.
//

#ifndef GMDS_RL_BLOCK_SET_H
#define GMDS_RL_BLOCK_SET_H
/*----------------------------------------------------------------------------*/
#include "gmds/ig/Mesh.h"
/*----------------------------------------------------------------------------*/

namespace gmds {
class RLBlockSet
{
 public:
	gmds::Mesh m_mesh;
	RLBlockSet(gmds::Mesh m_mesh);

	void createNode(double xA, double yA, double xB, double yB);
	
	void setFrame(double xMin, double yMin, double xMax, double yMax);

	int getNumberOfBlocks();

	void deleteBlock(const int faceID);

};
}
#endif     // GMDS_RL_BLOCK_SET_H
