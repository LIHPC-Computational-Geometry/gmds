//
// Created by Baptiste Bony on 22/04/22.
//

#ifndef BLOCKMESHER_CPP_RLBLOCKSET_H
#define BLOCKMESHER_CPP_RLBLOCKSET_H
/*----------------------------------------------------------------------------*/
#include "gmds/ig/Mesh.h"
/*----------------------------------------------------------------------------*/

namespace gmds {
// Action Delete Face (swap value activate variable)
class RLBlockSet
{
 public:
	RLBlockSet(gmds::Mesh mesh);

	void setFrame(gmds::Node xMin, gmds::Node yMin, gmds::Node xMax, gmds::Node yMax);

	int getNumberOfBlocks();

	void deleteBlock(const int faceID);

};
}
#endif     // BLOCKMESHER_CPP_RLBLOCKSET_H
