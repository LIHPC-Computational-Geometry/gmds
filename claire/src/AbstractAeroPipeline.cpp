//
// Created by rochec on 09/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AbstractAeroPipeline.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

/*
AbstractAeroPipeline::AbstractAeroPipeline(ParamsAero Aparams, gmds::Mesh && Am) :
  m_params(Aparams),
  m_m(std::move(Am))
{
	m_mesh = &m_m;
	m_isOver = false;
}
 */
/*------------------------------------------------------------------------*/
AbstractAeroPipeline::AbstractAeroPipeline(ParamsAero Aparams) :
  m_params(Aparams)
{
	m_isOver = false;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
bool AbstractAeroPipeline::getIsOver(){
	return m_isOver;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
std::vector<Node> AbstractAeroPipeline::AdjacentNodes(Node n){

	std::vector<Edge> adjacent_edges = n.get<Edge>() ;
	std::vector<Node> adjacent_nodes;
	for (auto e:adjacent_edges){
		TCellID ne_id = e.getOppositeNodeId(n);
		Node ne = m_mesh->get<Node>(ne_id);
		adjacent_nodes.push_back(ne);
	}

	return adjacent_nodes;
}
/*------------------------------------------------------------------------*/