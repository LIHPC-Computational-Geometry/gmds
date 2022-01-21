//
// Created by rochec on 21/01/2022.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/GradientComputation2D.h>
#include <limits>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/


GradientComputation2D::GradientComputation2D(Mesh *AMesh, Variable<double> *Adistance) {
	m_mesh = AMesh;
	m_distance = Adistance;
	m_gradient2D = m_mesh->newVariable<math::Vector3d ,GMDS_FACE>("gradient 2D");
}


/*------------------------------------------------------------------------*/
GradientComputation2D::STATUS GradientComputation2D::execute()
{


	return GradientComputation2D::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
math::Vector3d GradientComputation2D::getNormalVector(TCellID n0_id, TCellID n1_id){
	Node n0 = m_mesh->get<Node>(n0_id);
	Node n1 = m_mesh->get<Node>(n1_id);
	math::Point p0 = n0.point();
	math::Point p1 = n1.point();

	math::Vector3d Vecteur;
	Vecteur = p1-p0;

	math::Vector3d Vecteur_Normal;
	Vecteur_Normal.setX( - Vecteur.Y() ) ;
	Vecteur_Normal.setY( Vecteur.X() ) ;
	Vecteur_Normal.setZ(0);

	return Vecteur_Normal ;

}
/*------------------------------------------------------------------------*/