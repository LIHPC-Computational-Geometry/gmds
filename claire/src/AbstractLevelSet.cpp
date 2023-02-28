//
// Created by rochec on 10/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AbstractLevelSet.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AbstractLevelSet::AbstractLevelSet(Mesh *AMesh, int AmarkFrontNodes, Variable<double>* Adistance) {
	m_mesh = AMesh;
	m_markFrontNodes = AmarkFrontNodes;
	m_distance = Adistance;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
AbstractLevelSet::STATUS AbstractLevelSet::execute()
{
	initialisationDistances();

	// Tant qu'il y a des noeuds dans la map
	while(!m_DistanceMap.isEmpty()){
		TCellID n0_id;
		double v0;
		m_DistanceMap.getAndRemoveFirst(v0, n0_id);
		Node n0 = m_mesh->get<Node>(n0_id);
		math::Point p0 = n0.point();
		std::vector<Node> Neighbors_Nodes = getNeighbors(n0);
		for(auto const &n:Neighbors_Nodes){
			TCellID n_id = n.id();
			math::Point p = n.point();
			math::Vector3d vec = p-p0;
			double le = vec.norm();
			double ve ;
			getValue(n_id, ve);
			double ve_new = v0 + le ;
			if (ve_new < ve){
				setValue(n_id, ve_new);
				m_DistanceMap.update(ve, ve_new, n_id);
			}
		}
	}

	return AbstractLevelSet::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void AbstractLevelSet::initialisationDistances(){
	for (auto id:m_mesh->nodes()){
		if(m_mesh->isMarked<Node>(id, m_markFrontNodes)){
			m_DistanceMap.add(0, id);
			m_distance->set(id, 0);
		}
		else{
			m_DistanceMap.add(std::numeric_limits<double>::max(), id);
			m_distance->set(id, std::numeric_limits<double>::max());
		}
	}
}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void AbstractLevelSet::getValue(TCellID n_id, double &v0){
	v0 = m_distance->value(n_id);
}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void AbstractLevelSet::setValue(TCellID n_id, double v0){
	m_distance->value(n_id) = v0 ;
}
/*-------------------------------------------------------------------*/