//
// Created by rochec on 10/02/2022.
//
/*------------------------------------------------------------------------*/
#include <gmds/aero/LevelSetExtended.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

LevelSetExtended::LevelSetExtended(Mesh *AMesh, TInt AmarkFrontNodes, Variable<double>* Adistance) :
  AbstractLevelSet(AMesh,AmarkFrontNodes, Adistance)
{

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
std::vector<Node> LevelSetExtended::getNeighbors(Node n){
	std::vector<Node> vec_Neighbors;
	std::vector<Edge> adjacent_edges = n.get<Edge>() ;
	for(auto const& e:adjacent_edges) {
		TCellID ne_id = e.getOppositeNodeId(n);
		Node ne = m_mesh->get<Node>(ne_id);
		vec_Neighbors.push_back(ne);
		std::vector<Edge> adjacent_edges_2 = ne.get<Edge>() ;
		for (auto const& e2:adjacent_edges_2){
			TCellID ne2_id = e2.getOppositeNodeId(ne);
			Node ne2 = m_mesh->get<Node>(ne2_id);
			vec_Neighbors.push_back(ne2);
		}
	}
	return vec_Neighbors;
}

/*------------------------------------------------------------------------*/