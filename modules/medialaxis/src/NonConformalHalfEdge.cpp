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
std::vector<Node> NonConformalHalfEdge::getOrderedNodes()
{
	std::vector<Node> ordered_nodes;
	ordered_nodes.push_back(m_first_node);
	Node n;
	for (int i = 0; i < m_conformal_edges.size()-1; i++)
	{
		n = getCommonNode(m_conformal_edges[i],m_conformal_edges[i+1]);
		ordered_nodes.push_back(n);
	}
	ordered_nodes.push_back(m_end_node);
	return ordered_nodes;
}
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/
