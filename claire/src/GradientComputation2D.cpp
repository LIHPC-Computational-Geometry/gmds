//
// Created by rochec on 21/01/2022.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/GradientComputation2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/


GradientComputation2D::GradientComputation2D(Mesh *AMesh, Variable<double>* Adistance) {
	m_mesh = AMesh;
	m_distance = Adistance;
	m_gradient2D = m_mesh->newVariable<math::Vector3d ,GMDS_FACE>("GMDS_Gradient_2D");
}


/*------------------------------------------------------------------------*/
GradientComputation2D::STATUS GradientComputation2D::execute()
{

	for(auto face_id:m_mesh->faces()) {
		math::Vector3d Gradient = computeGradientOnSimpleFace(face_id);
		m_gradient2D->set(face_id, Gradient) ;
	}

	return GradientComputation2D::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
math::Vector3d GradientComputation2D::getNormalVector(TCellID n0_id, TCellID n1_id){
	Node n0 = m_mesh->get<Node>(n0_id);
	Node n1 = m_mesh->get<Node>(n1_id);
	math::Point p0 = n0.point();
	math::Point p1 = n1.point();

	math::Vector3d Vecteur = p1-p0;

	math::Vector3d Vecteur_Normal;
	Vecteur_Normal.setX( - Vecteur.Y() ) ;
	Vecteur_Normal.setY( Vecteur.X() ) ;
	Vecteur_Normal.setZ(0);

	return Vecteur_Normal ;

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
math::Vector3d GradientComputation2D::computeGradientOnSimpleFace(TCellID face_id){
	math::Vector3d Gradient ;
	Face face = m_mesh->get<Face>(face_id);
	std::vector<TCellID> face_nodes_ids = face.getIDs<Node>();
	TCellID vi_id = face_nodes_ids[0];
	TCellID vj_id = face_nodes_ids[2];
	TCellID vk_id = face_nodes_ids[1];
	double di = m_distance->value(vi_id) ;
	double dj = m_distance->value(vj_id) ;
	double dk = m_distance->value(vk_id) ;
	double At = face.area();

	math::Vector3d vki_ortho = getNormalVector(vk_id, vi_id) ;
	math::Vector3d vij_ortho = getNormalVector(vi_id, vj_id) ;

	Gradient = vki_ortho*(dj-di)/(2.0*At) + vij_ortho*(dk-di)/(2.0*At) ;

	//std::cout << "Gradient : " << Gradient << std::endl;

	return Gradient;
}
/*------------------------------------------------------------------------*/
