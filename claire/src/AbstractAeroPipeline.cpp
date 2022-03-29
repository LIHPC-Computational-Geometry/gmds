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
	m_meshTet = &m_m;
	m_isOver = false;
}
 */
/*------------------------------------------------------------------------*/
AbstractAeroPipeline::AbstractAeroPipeline(ParamsAero Aparams) :
  m_params(Aparams),
  m_meshTet(nullptr),
  m_meshHex(nullptr),
  m_manager(new cad::FACManager()),
  m_linker_TG(new cad::GeomMeshLinker())
{

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
AbstractAeroPipeline::~AbstractAeroPipeline(){
	if(m_manager){
		delete m_manager;
	}
	if(m_linker_TG){
		delete m_linker_TG;
	}
	if(m_meshTet){
		delete m_meshTet;
	}
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
std::vector<Node> AbstractAeroPipeline::AdjacentNodes(Node n){

	std::vector<Edge> adjacent_edges = n.get<Edge>() ;
	std::vector<Node> adjacent_nodes;
	for (auto e:adjacent_edges){
		TCellID ne_id = e.getOppositeNodeId(n);
		Node ne = m_meshTet->get<Node>(ne_id);
		adjacent_nodes.push_back(ne);
	}

	return adjacent_nodes;
}
/*------------------------------------------------------------------------*/