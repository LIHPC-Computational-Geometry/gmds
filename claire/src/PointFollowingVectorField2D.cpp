//
// Created by rochec on 26/01/2022.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/PointFollowingVectorField2D.h>
#include <limits>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

PointFollowingVectorField2D::PointFollowingVectorField2D(Mesh *AMesh, math::Point A_Pstart, double A_distance, Variable<math::Vector3d>* A_gradient2D) {
	m_mesh = AMesh;
	m_Pstart = A_Pstart;
	m_distance = A_distance;
	m_gradient2D = A_gradient2D;
	m_discrete_path.push_back(m_Pstart);
}


/*------------------------------------------------------------------------*/
PointFollowingVectorField2D::STATUS PointFollowingVectorField2D::execute()
{
	/*
	math::Point P(0.5, 0.5, 0.0);
	TCellID face_id = inWhichTriangle(P);
	std::cout << "Le point P " << P << " appartient au triangle d'id : " << face_id << std::endl;

	P.setXYZ(0.3, 0.3, 0.0);
	face_id = inWhichTriangle(P);
	std::cout << "Le point P " << P << " appartient au triangle d'id : " << face_id << std::endl;
	 */

	m_Pend = m_Pstart;
	double minLenght = minEdgeLenght();
	std::cout << "min : " << minLenght << std::endl;

	while (m_distance != 0){
		TCellID pointFace_id = inWhichTriangle(m_Pend) ;
		math::Vector3d Grad_local = m_gradient2D->value(pointFace_id) ;
		if (m_distance >= minLenght){
			m_Pend = m_Pend + minLenght*Grad_local;
			//m_distance = m_distance - Grad_local.norm();
			m_distance = m_distance - minLenght;
		}
		else{
			m_Pend = m_Pend + m_distance*Grad_local;
			m_distance = 0;
		}
		//std::cout << "Point intermédiaire : " << m_Pend << std::endl;
		m_discrete_path.push_back(m_Pend);
	}

	std::cout << "Point final : " << m_Pend << std::endl;

	writeDiscretePathInVTK();

	return PointFollowingVectorField2D::SUCCESS;
}
/*------------------------------------------------------------------------*/





/*------------------------------------------------------------------------*/
bool PointFollowingVectorField2D::isInTriangle(TCellID face_id, math::Point M){
	bool isInFace(false);
	Face face = m_mesh->get<Face>(face_id);
	std::vector<TCellID> face_nodes_ids = face.getIDs<Node>();
	TCellID vi_id = face_nodes_ids[0];
	TCellID vj_id = face_nodes_ids[2];
	TCellID vk_id = face_nodes_ids[1];

	Node vi = m_mesh->get<Node>(vi_id);
	Node vj = m_mesh->get<Node>(vj_id);
	Node vk = m_mesh->get<Node>(vk_id);

	math::Point vi_coord = vi.point();
	math::Point vj_coord = vj.point();
	math::Point vk_coord = vk.point();

	math::Vector3d vij = vj_coord-vi_coord ;
	math::Vector3d vjk = vk_coord-vj_coord ;
	math::Vector3d vki = vi_coord-vk_coord ;
	math::Vector3d viM = M-vi_coord ;
	math::Vector3d vjM = M-vj_coord ;
	math::Vector3d vkM = M-vk_coord ;

	double d1 = ( vij.cross(viM) ).dot( viM.cross(-vki) ) ;
	double d2 = ( -vij.cross(vjM) ).dot( vjM.cross(vjk) ) ;
	double d3 = ( vki.cross(vkM) ).dot( vkM.cross(-vjk) ) ;

	if (d1 >= 0 && d2 >= 0 && d3 >= 0) {
		isInFace = true;
	}
	return isInFace;
}
/*------------------------------------------------------------------------*/




/*------------------------------------------------------------------------*/
TCellID PointFollowingVectorField2D::inWhichTriangle(math::Point M){
	TCellID face_id;
	bool isInFace(false);
	for (auto f_id:m_mesh->faces()){
		if(!isInFace){
			isInFace = isInTriangle(f_id, M);
			if(isInFace){
				face_id = f_id;
			}
		}
	}
	return face_id;
}
/*------------------------------------------------------------------------*/




/*------------------------------------------------------------------------*/
double PointFollowingVectorField2D::minEdgeLenght(){
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





/*------------------------------------------------------------------------*/
void PointFollowingVectorField2D::writeDiscretePathInVTK(){

	// First, we create the file where we are going to store the info
	std::ofstream stream= std::ofstream("GMDS_discrete_path.vtk", std::ios::out);
	//set the numerical precision (number of digits)
	stream.precision(15);
	//Header indicating which type of file it is
	//stream << "# Claire PHD Debug Version 1.0\n\n";

	stream << "# vtk DataFile Version 2.0\n";
	stream << "Test moche écriture VTK\n";
	stream << "ASCII\n\n";

	stream << "DATASET UNSTRUCTURED_GRID\n\n";

	stream << "POINTS ";
	stream << m_discrete_path.size() ;
	stream << " float\n";
	for (int i=0; i< m_discrete_path.size(); i++){
		stream << m_discrete_path[i].X() << " " << m_discrete_path[i].Y() << " " << m_discrete_path[i].Z() << "\n";
	}

	stream << "\n";

	stream << "POINT_DATA " << m_discrete_path.size() << "\n" ;
	stream << "SCALARS GMDS_discrete_path float 1\n";
	stream << "LOOKUP_TABLE default\n";
	for (int i=0; i< m_discrete_path.size(); i++){
		stream << 1.0 << "\n";
	}

	stream.close();
}
/*------------------------------------------------------------------------*/