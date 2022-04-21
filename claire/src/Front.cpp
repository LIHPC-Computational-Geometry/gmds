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
std::vector<TCellID> Front::getNeighbors(TCellID n_id){
	return m_NodeNeighbors[n_id];
}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
int Front::getNodeType(TCellID n_id){
	return m_NodeInfo[n_id].nodeType;
}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
TCellID Front::getNextNode(TCellID n_id, TCellID n_neighbors_id){
	return m_NodeInfo[n_neighbors_id].next_nodes[n_id];
}
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


/*-------------------------------------------------------------------*/
void Front::initializeNodeNeighbors(Mesh *m){
	for (auto e_id:m_edgesId){
		Edge e = m->get<Edge>(e_id);
		std::vector<Node> nodes = e.get<Node>();
		m_NodeNeighbors[nodes[0].id()].push_back(nodes[1].id());
		m_NodeNeighbors[nodes[1].id()].push_back(nodes[0].id());
	}
}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void Front::initializeNodeType(Mesh *m, std::map<TCellID, TCellID> map_idealpositions){

	Variable<int>* couche_id = m->getVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	initializeNodeNeighbors(m);

	for (auto n_id:m_nodesId){
		// Initialisation des positions idéales
		m_NodeInfo[n_id].ideal_node_id = map_idealpositions[n_id];
		// On initialise tous les types à 0 (normal)
		m_NodeInfo[n_id].nodeType = 0;
		for (auto n_neighbor_id:m_NodeNeighbors[n_id]){
			m_NodeInfo[n_id].next_nodes[n_neighbor_id] = map_idealpositions[n_id];
		}
	}

}
/*-------------------------------------------------------------------*/