//
// Created by rochec on 13/01/2022.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/LevelSet2D.h>
#include <limits>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/


LevelSet2D::LevelSet2D(Mesh *AMesh, int AmarkFrontNodes) {
	m_mesh = AMesh;
	m_markFrontNodes = AmarkFrontNodes;
	m_distance = m_mesh->newVariable<double,GMDS_NODE>("distance");

}


/*------------------------------------------------------------------------*/
LevelSet2D::STATUS LevelSet2D::execute()
{
	initialisationDistances();

	// Tant qu'il y a des noeuds dans la map
	while(!m_DistanceMap.isEmpty()){
		TCellID n0_id;
		double v0;
		m_DistanceMap.getAndRemoveFirst(v0, n0_id);
		Node n0 = m_mesh->get<Node>(n0_id);
		std::vector<Edge> adjacent_edges = n0.get<Edge>() ;
		for(auto e:adjacent_edges){
			double le =e.length();
			TCellID ne_id = e.getOppositeNodeId(n0);
			double ve ;
			getValue(ne_id, ve);
			double ve_new = v0 + le ;
			if (ve_new < ve){
				setValue(ne_id, ve_new);
				m_DistanceMap.update(ve, ve_new, ne_id);
			}


			/*
			// --- MODIFICATION NAIVE : On regarde les noeuds à 2 edges d'écart ---
			Node ne = m_mesh->get<Node>(ne_id);
			std::vector<Edge> adjacent_edges_2 = ne.get<Edge>() ;
			for(auto e2:adjacent_edges_2) {
				TCellID ne2_id = e2.getOppositeNodeId(ne);
				Node ne2 = m_mesh->get<Node>(ne2_id);
				math::Point p0 = n0.point();
				math::Point pe2 = ne2.point();
				math::Vector3d Vec = pe2-p0 ;
				double le2 = Vec.norm();
				double ve;
				getValue(ne2_id, ve);
				double ve_new = v0 + le2;
				if (ve_new < ve) {
					setValue(ne2_id, ve_new);
					m_DistanceMap.update(ve, ve_new, ne2_id);
				}
			}
			// --- FIN MODIFICATION NAIVE ---
			*/




		}
	}

	return LevelSet2D::SUCCESS;
}
/*------------------------------------------------------------------------*/




/*-------------------------------------------------------------------*/
void LevelSet2D::initialisationDistances(){
	for (auto id:m_mesh->nodes()){
		if(m_mesh->isMarked<Node>(id, m_markFrontNodes)){
			m_DistanceMap.add(0, id);
			m_distance->set(id, 0);
		}
		else{
			m_DistanceMap.add(std::numeric_limits<double>::max(), id);
			m_distance->set(id, std::numeric_limits<double>::max());
			//m_DistanceMap.add(1000, id);
			//m_distance->set(id, 1000);
		}
	}
};
/*-------------------------------------------------------------------*/





/*-------------------------------------------------------------------*/
void LevelSet2D::getValue(TCellID n_id, double &v0){
	v0 = m_distance->value(n_id);
};
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void LevelSet2D::setValue(TCellID n_id, double v0){
	m_distance->value(n_id) = v0 ;
};
/*-------------------------------------------------------------------*/