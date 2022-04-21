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
		// On initialise tous les types à 0 (normal)
		m_NodeInfo[n_id].nodeType = 0;
		// On cherche les arêtes voisines à n sur le front
		/*
		Node n = m->get<Node>(n_id);
		std::vector<Edge> edges = n.get<Edge>();
		for (auto e:edges){
			std::vector<Node> nodes = e.get<Node>();
			Node n_opp = e.getOppositeNode(n);
			if ( couche_id->value(nodes[0].id()) == couche_id->value(nodes[1].id()) ){
				m_NodeInfo[n_id].next_nodes[n_opp.id()] = map_idealpositions[n_id];
			}
		}
		 */
		for (auto n_neighbor_id:m_NodeNeighbors[n_id]){
			m_NodeInfo[n_id].next_nodes[n_neighbor_id] = map_idealpositions[n_id];
		}
	}

}
/*-------------------------------------------------------------------*/