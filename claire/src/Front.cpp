//
// Created by rochec on 20/04/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/Front.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
Front::Front() {

}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void Front::setFrontID(int layer_id){
	m_FrontID = layer_id;
}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void Front::setMultipleNode(TCellID n_id){
	m_NodeInfo[n_id].nodeType = 1;
}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void Front::setContractedNode(TCellID n_id){
	m_NodeInfo[n_id].nodeType = 2;
}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void Front::setNextNode(TCellID n_id, TCellID n_neighbor_id, TCellID n_new_id){
	m_NodeInfo[n_id].next_nodes[n_neighbor_id] = n_new_id ;
}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void Front::setNonFusionable(TCellID n_id){
	m_NodeInfo[n_id].isFusionable = false;
}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void Front::setNonMultiplicable(TCellID n_id){
	m_NodeInfo[n_id].isMultiplicable = false;
}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
int Front::getFrontID(){
	return m_FrontID;
}
/*-------------------------------------------------------------------*/


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
bool Front::isFusionable(TCellID n_id){
	return m_NodeInfo[n_id].isFusionable;
}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
bool Front::isMultiplicable(TCellID n_id){
	return m_NodeInfo[n_id].isMultiplicable;
}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
TCellID Front::getIdealNode(TCellID n_id){
	return m_NodeInfo[n_id].ideal_node_id;
}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
TCellID Front::getNextNode(TCellID n_id, TCellID n_neighbors_id){
	return m_NodeInfo[n_id].next_nodes[n_neighbors_id];
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

	m_FrontID = layer_id;

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
		// Initialisation de la fusionabilité ou multiplicabilité des noeuds
		m_NodeInfo[n_id].isFusionable = true;
		m_NodeInfo[n_id].isMultiplicable = true;
	}

}
/*-------------------------------------------------------------------*/