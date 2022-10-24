//
// Created by rochec on 24/10/2022.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/DiffusionEquation2D.h>
#include <Eigen/Sparse>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

DiffusionEquation2D::DiffusionEquation2D(Mesh *AMesh, int AmarkFrontNodes_int, int AmarkFrontNodes_out) {
	m_mesh = AMesh;
	m_markNodes_int = AmarkFrontNodes_int;
	m_markNodes_out = AmarkFrontNodes_out;
}


/*------------------------------------------------------------------------*/
DiffusionEquation2D::STATUS
DiffusionEquation2D::execute()
{

	return DiffusionEquation2D::SUCCESS;
}
/*------------------------------------------------------------------------*/






//
//		GEOMETRY FUNCTIONS
//

/*------------------------------------------------------------------------*/
// Construction des 3 fonctions de base
// (phi_hat_1, phi_hat_2, phi_hat_3)
double DiffusionEquation2D::phi_hat(int hat_i, math::Point X_hat)
{
	double value;
	switch (hat_i)
	{
	case 1:
		value=1-X_hat.X()-X_hat.Y();
		break;
	case 2:
		value=X_hat.X();
		break;
	case 3:
		value=X_hat.Y();
		break;
	default:
		std::cout << "Please choose a valid number of hat_i for the reference triangle (1, 2 or 3)" << std::endl;
		abort();
	}
	return value;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
// Construction des 3 gradients des fonctions de base
// (gradphihat1, gradphihat2, gradphihat3)
// (indépendant du vecteur X_hat)
Eigen::Vector2d DiffusionEquation2D::grad_phi_hat(int hat_i)
{
	Eigen::Vector2d gradphi;
	switch (hat_i)
	{
	case 1:
		gradphi(0) = -1; gradphi(1) = -1;
		break;
	case 2:
		gradphi(0) = 1; gradphi(1) = 0;
		break;
	case 3:
		gradphi(0) = 0; gradphi(1) = 1;
		break;
	default:
		std::cout << "Please choose a valid number of hati for the reference triangle (1, 2 or 3)" << std::endl;
		abort();
	}
	return gradphi;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
// Construction de la fonction de transformation
// du triangle de référence K_hat en K
Eigen::Vector2d DiffusionEquation2D::FK(TCellID triK_id, math::Point X_hat)
{
	Face triK = m_mesh->get<Face>(triK_id);
	std::vector<Node> tri_nodes = triK.get<Node>() ;

	Eigen::Vector2d X;
	X(0) = tri_nodes[0].X()* phi_hat(1, X_hat) + tri_nodes[1].X()* phi_hat(2, X_hat) + tri_nodes[2].X()* phi_hat(3, X_hat);
	X(1) = tri_nodes[0].Y()* phi_hat(1, X_hat) + tri_nodes[1].Y()* phi_hat(2, X_hat) + tri_nodes[2].Y()* phi_hat(3, X_hat);

	return X;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
// Construction de la jacobienne de la fonction de transformation
// du triangle de référence K_hat en K
// (indépendant du vecteur X_hat)
Eigen::Matrix2d DiffusionEquation2D::JFK(TCellID triK_id)
{
	Face triK = m_mesh->get<Face>(triK_id);
	std::vector<Node> tri_nodes = triK.get<Node>() ;

	Eigen::Matrix2d X;
	X(0,0) = tri_nodes[1].X() - tri_nodes[0].X();
	X(0,1) = tri_nodes[2].X() - tri_nodes[0].X();
	X(1,0) = tri_nodes[1].Y() - tri_nodes[0].Y();
	X(1,1) = tri_nodes[2].Y() - tri_nodes[0].Y();

	return X;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
// Construction de la valeur absolue du déterminant de la jacobienne
// de la fonction de transformation
// du triangle de référence hatK en K
// (indépendant du vecteur hatX)
double DiffusionEquation2D::absdetJFK(TCellID triK_id)
{
	return fabs((JFK(triK_id)).determinant());
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
// Construction du gradient
// (gradphi_i) avec i qui est le hati-ème sommet de K
// appliqué en F(hatX)
Eigen::Vector2d DiffusionEquation2D::grad_phi(int i_hat, TCellID triK_id)
{
	return (((JFK(triK_id)).inverse()).transpose())*grad_phi_hat(i_hat);
}
/*------------------------------------------------------------------------*/