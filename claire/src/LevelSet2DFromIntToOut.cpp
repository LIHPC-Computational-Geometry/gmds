//
// Created by rochec on 19/01/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/LevelSet2DFromIntToOut.h>
//#include <limits>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/


LevelSet2DFromIntToOut::LevelSet2DFromIntToOut(Mesh *AMesh, int AmarkFrontNodesInt, int AmarkFrontNodesOut) {
	m_mesh = AMesh;
	m_markFrontNodesInt = AmarkFrontNodesInt;
	m_markFrontNodesOut = AmarkFrontNodesOut;
	m_distance = m_mesh->newVariable<double,GMDS_NODE>("distance");

}