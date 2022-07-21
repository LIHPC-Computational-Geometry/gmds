//
// Created by rochec on 20/07/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AbstractSmoothLineSweeping_2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AbstractSmoothLineSweeping_2D::AbstractSmoothLineSweeping_2D(Blocking2D::Block* B, int nb_max_it) {
	m_B = B;
	m_nb_max_iterations = nb_max_it;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
AbstractSmoothLineSweeping_2D::STATUS AbstractSmoothLineSweeping_2D::execute()
{
	//std::map<TCellID, math::Point> old_coords;
	Array2D<TCoord>* old_coord_x;
	Array2D<TCoord>* old_coord_y;
	int Nx = m_B->getNbDiscretizationI();
	int Ny = m_B->getNbDiscretizationJ();

	for (int i=0; i <= Nx; i++)
	{
		for (int j=0; j <= Ny; j++)
		{
			Node n = (*m_B)(i, j);
			(*old_coord_x)(i,j) = n.X() ;
			(*old_coord_y)(i,j) = n.Y() ;
		}
	}

	return AbstractSmoothLineSweeping_2D::SUCCESS;
}
/*------------------------------------------------------------------------*/







/*------------------------------------------------------------------------*/
// Fonction FindMidBranche : Si on considère une branche composée de 3 points, cette fonction retourne le point
//                            positionné au milieu de cette branche
// En entrée : A, B, C -> 3 points
// En sortie : le point milieu
math::Point AbstractSmoothLineSweeping_2D::FindMidBranche(const math::Point A, const math::Point B, const math::Point C) {
	math::Point Point_Milieu ;
	math::Vector3d Vec_AB = B-A ;
	math::Vector3d Vec_BC = C-B ;
	double norme_1 = Vec_AB.norm() ;
	double norme_2 = Vec_BC.norm() ;
	double norme_branche = norme_1 + norme_2 ;
	double norme_milieu = norme_branche / 2.0 ;

	if (norme_milieu <= norme_1){
		Vec_AB.normalize();
		Point_Milieu = A + norme_milieu*Vec_AB ;
	}
	else if (norme_milieu > norme_1){
		math::Vector3d Vec_CB = - Vec_BC ;
		Vec_CB.normalize();
		Point_Milieu = C + norme_milieu*Vec_CB ;
	}

	return Point_Milieu;
}
/*------------------------------------------------------------------------*/