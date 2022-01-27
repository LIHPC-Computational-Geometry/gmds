//
// Created by rochec on 27/01/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/PointFollowingVectorField3D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

PointFollowingVectorField3D::PointFollowingVectorField3D(Mesh *AMesh, math::Point A_Pstart, double A_distance, Variable<math::Vector3d>* A_gradient3D) {
	m_mesh = AMesh;
	m_Pstart = A_Pstart;
	m_distance = A_distance;
	m_gradient3D = A_gradient3D;
}


/*------------------------------------------------------------------------*/
PointFollowingVectorField3D::STATUS PointFollowingVectorField3D::execute()
{
	m_Pend = m_Pstart;
	double minLenght = minEdgeLenght();
	std::cout << "min : " << minLenght << std::endl;

	while (m_distance != 0){
		TCellID pointFace_id = inWhichTetra(m_Pend) ;
		math::Vector3d Grad_local = m_gradient3D->value(pointFace_id) ;
		Grad_local.normalize();
		if (m_distance >= minLenght){
			m_Pend = m_Pend + minLenght*Grad_local;
			m_distance = m_distance - minLenght;
		}
		else{
			m_Pend = m_Pend + m_distance*Grad_local;
			m_distance = 0;
		}
		std::cout << "Point intermédiaire : " << m_Pend << std::endl;
	}

	std::cout << "Point final : " << m_Pend << std::endl;

	return PointFollowingVectorField3D::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
bool PointFollowingVectorField3D::sameSide(math::Point v1, math::Point v2, math::Point v3, math::Point v4, math::Point p){

	bool isOnSameSide(false);

	math::Vector3d vec_12 = v2-v1;
	math::Vector3d vec_13 = v3-v1;
	math::Vector3d vec_14 = v4-v1;
	math::Vector3d vec_1p = p-v1;

	math::Vector3d normal = vec_12.cross(vec_13);

	double dotv4 = normal.dot(vec_14);
	double dotP = normal.dot(vec_1p);
	double signe = dotv4*dotP;

	if (signe >= 0){
		isOnSameSide = true;
	}

	return isOnSameSide;
}
/*------------------------------------------------------------------------*/




/*------------------------------------------------------------------------*/
bool PointFollowingVectorField3D::isInTetra(TCellID region_id, math::Point M){

	bool isInRegion(false);

	Region region = m_mesh->get<Region>(region_id);
	std::vector<TCellID> region_nodes_ids = region.getIDs<Node>();
	TCellID vi_id = region_nodes_ids[0];
	TCellID vj_id = region_nodes_ids[1];
	TCellID vk_id = region_nodes_ids[2];
	TCellID vh_id = region_nodes_ids[3];

	Node vi = m_mesh->get<Node>(vi_id);
	Node vj = m_mesh->get<Node>(vj_id);
	Node vk = m_mesh->get<Node>(vk_id);
	Node vh = m_mesh->get<Node>(vh_id);

	math::Point vi_coord = vi.point();
	math::Point vj_coord = vj.point();
	math::Point vk_coord = vk.point();
	math::Point vh_coord = vh.point();

	bool sameSide1 = sameSide(vi_coord, vj_coord, vk_coord, vh_coord, M);
	bool sameSide2 = sameSide(vj_coord, vk_coord, vh_coord, vi_coord, M);
	bool sameSide3 = sameSide(vk_coord, vh_coord, vi_coord, vj_coord, M);
	bool sameSide4 = sameSide(vh_coord, vi_coord, vj_coord, vk_coord, M);

	if (sameSide1 && sameSide2 && sameSide3 && sameSide4) {
		isInRegion = true;
	}

	return isInRegion;
}
/*------------------------------------------------------------------------*/




/*------------------------------------------------------------------------*/
TCellID PointFollowingVectorField3D::inWhichTetra(math::Point M){
	TCellID region_id;
	bool isInRegion(false);
	for (auto reg_id:m_mesh->regions()){
		if(!isInRegion){
			isInRegion = isInTetra(reg_id, M);
			if(isInRegion){
				region_id = reg_id;
			}
		}
	}
	return region_id;
}
/*------------------------------------------------------------------------*/




/*------------------------------------------------------------------------*/
double PointFollowingVectorField3D::minEdgeLenght(){
	// Initialisation avec une arête prise au hasard. Pb : si l'arête 0 a été
	// retirée
	Edge edge_0 = m_mesh->get<Edge>(0);
	double minLenght(edge_0.length());
	for (auto edge_id:m_mesh->edges()){
		Edge edge = m_mesh->get<Edge>(edge_id);
		if(edge.length() < minLenght){
			minLenght = edge.length() ;
			//std::cout << "Edge id : " << edge_id << ", Taille : " << minLenght << std::endl;
		}
	}
	return minLenght;
}
/*------------------------------------------------------------------------*/