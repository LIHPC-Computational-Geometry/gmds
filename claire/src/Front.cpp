//
// Created by rochec on 20/04/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/Front.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

Front::Front() {

}



/*-------------------------------------------------------------------*/
std::vector<TCellID> Front::getNodes(){
	return m_nodesId;
};
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
std::vector<TCellID> Front::getEdges(){
	return m_edgesId;
};
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void Front::addNodeId(TCellID n_id){
	m_nodesId.push_back(n_id);
}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void Front::addEdgeId(TCellID e_id){
	m_edgesId.push_back(e_id);
}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void Front::initializeFromLayerId(Mesh *m, int layer_id){
	Variable<int>* couche_id = m->getVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	for (auto n_id:m->nodes())
	{
		if(couche_id->value(n_id) == layer_id)
		{
			m_nodesId.push_back(n_id);
		}
	}

	for (auto e_id:m->edges()){
		Edge e = m->get<Edge>(e_id);
		std::vector<Node> nodes = e.get<Node>();
		if ( couche_id->value(nodes[0].id()) == layer_id
		    && couche_id->value(nodes[1].id()) == layer_id )
		{
			m_edgesId.push_back(e_id);
		}
	}
}
/*-------------------------------------------------------------------*/