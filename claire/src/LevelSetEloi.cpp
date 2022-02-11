//
// Created by rochec on 10/02/2022.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/LevelSetEloi.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

LevelSetEloi::LevelSetEloi(Mesh *AMesh, int AmarkFrontNodes, Variable<double> *Adistance) :
  AbstractLevelSet(AMesh,AmarkFrontNodes, Adistance)
{

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
std::vector<Node> LevelSetEloi::getNeighbors(Node n){
	std::vector<Node> vec_Neighbors;
	std::vector<Edge> adjacent_edges = n.get<Edge>() ;
	for(auto e:adjacent_edges) {
		TCellID ne_id = e.getOppositeNodeId(n);
		Node ne = m_mesh->get<Node>(ne_id);
		vec_Neighbors.push_back(ne);
	}
	return vec_Neighbors;
}

/*------------------------------------------------------------------------*/