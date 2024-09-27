//
// Created by chenyt on 26/09/24.
//
#include "gmds/medialaxis/NonConformalHalfEdge.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
NonConformalHalfEdge::NonConformalHalfEdge(gmds::TCellID AId, gmds::Face AFace, std::vector<Edge> AEdges)
{
	m_id = AId;
	m_face = AFace;
	m_conformal_edges = AEdges;
	m_next = -1;
}
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
