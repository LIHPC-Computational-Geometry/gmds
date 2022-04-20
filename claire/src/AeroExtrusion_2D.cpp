//
// Created by rochec on 14/04/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AeroException.h>
#include <gmds/claire/AeroExtrusion_2D.h>
#include <gmds/claire/AdvectedPointRK4_2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AeroExtrusion_2D::AeroExtrusion_2D(Mesh *AMeshT, Mesh *AMeshQ) {
	m_meshT = AMeshT;
	m_meshQ = AMeshQ;
}


/*------------------------------------------------------------------------*/
AeroExtrusion_2D::STATUS
AeroExtrusion_2D::execute()
{
	// Exemple exception
	//if(m_mesh==NULL)
	//	throw AeroException("ERROR: Invalid mesh pointer");

	return AeroExtrusion_2D::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
Front
AeroExtrusion_2D::Compute1stLayer(Variable<double>* A_distance, double dist_cible){
	Front First_Front;

	return First_Front;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
std::map<TCellID, math::Point>
   AeroExtrusion_2D::ComputeIdealPositions(Front AFront, double dist_cible, Variable<double>* A_distance, Variable<math::Vector3d>* A_vectors)
{
	std::map<TCellID, math::Point> map_optnexpoint;
	std::vector<TCellID> front_nodes = AFront.getNodes();

	for (auto n_id:front_nodes){
		Node n = m_meshQ->get<Node>(n_id);
		math::Point M = n.point();
		AdvectedPointRK4_2D advpoint(m_meshT, M, dist_cible, A_distance, A_vectors);
		advpoint.execute();
		//math::Point P = advpoint.getPend();
		map_optnexpoint[n_id] = advpoint.getPend() ;
		}

	return map_optnexpoint;
}
/*------------------------------------------------------------------------*/