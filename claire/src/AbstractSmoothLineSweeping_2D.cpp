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
	m_Nx = m_B->getNbDiscretizationI();
	m_Ny = m_B->getNbDiscretizationJ();
	m_tol = pow(10,-6);
	m_theta = 0.0;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
AbstractSmoothLineSweeping_2D::STATUS AbstractSmoothLineSweeping_2D::execute()
{
	Array2D<TCoord>* old_coord_x;
	Array2D<TCoord>* old_coord_y;

	double err = pow(10,6);
	int iteration(0);

	while ( err > m_tol && iteration < m_nb_max_iterations )
	{
		One_Step_Smoothing();

		err = L2_norm_relative_error(old_coord_x, old_coord_y);
		iteration++;
		Update_old_coords(old_coord_x, old_coord_y);
	}

	return AbstractSmoothLineSweeping_2D::SUCCESS;
}
/*------------------------------------------------------------------------*/




/*------------------------------------------------------------------------*/
void AbstractSmoothLineSweeping_2D::One_Step_Smoothing(){

	for (int i=0; i < m_Nx; i++)
	{
		for (int j=0; j < m_Ny; j++)
		{
			Node n = (*m_B)(i, j);
			math::Point new_ideal_pos = ComputeNewPosition(i,j) ;

			// Damping
			n.setPoint( m_theta*n.point() + (1.0-m_theta)*new_ideal_pos );
		}
	}

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
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


/*------------------------------------------------------------------------*/
double AbstractSmoothLineSweeping_2D::L2_norm_relative_error(Array2D<TCoord>* old_coord_x, Array2D<TCoord>* old_coord_y)
{
	double num(0.0);
	double denom(0.0);

	int Nx = m_B->getNbDiscretizationI();
	int Ny = m_B->getNbDiscretizationJ();
	for (int i=0; i <= Nx; i++) {
		for (int j = 0; j <= Ny; j++) {
			Node n = (*m_B)(i, j);
			(*old_coord_x)(i, j) = n.X();
			(*old_coord_y)(i, j) = n.Y();

			double err_loc = sqrt( pow((*old_coord_x)(i, j) - n.X() , 2) + pow((*old_coord_y)(i, j) - n.Y(),2) );

			num += pow(err_loc, 2);
			denom += pow((*old_coord_x)(i, j), 2) + pow((*old_coord_y)(i, j), 2) ;
		}
	}

	return sqrt(num/denom);
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AbstractSmoothLineSweeping_2D::Update_old_coords(Array2D<TCoord>* old_coord_x, Array2D<TCoord>* old_coord_y)
{
	for (int i=0; i < m_Nx; i++)
	{
		for (int j=0; j < m_Ny; j++)
		{
			Node n = (*m_B)(i, j);
			(*old_coord_x)(i,j) = n.X() ;
			(*old_coord_y)(i,j) = n.Y() ;
		}
	}
}
/*------------------------------------------------------------------------*/