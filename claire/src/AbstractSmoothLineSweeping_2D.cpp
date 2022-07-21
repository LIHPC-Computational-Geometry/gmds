//
// Created by rochec on 20/07/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AbstractSmoothLineSweeping_2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AbstractSmoothLineSweeping_2D::AbstractSmoothLineSweeping_2D(Blocking2D::Block* B, int nb_max_it):
  	m_B(B),
  	m_nb_max_iterations(nb_max_it),
  	m_Nx(m_B->getNbDiscretizationI()-1),
	m_Ny(m_B->getNbDiscretizationJ()-1),
  	m_old_coord_x(m_Nx,m_Ny),
	m_old_coord_y(m_Nx,m_Ny)
{
	m_tol = pow(10,-6);
	m_theta = 0.0;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
AbstractSmoothLineSweeping_2D::STATUS AbstractSmoothLineSweeping_2D::execute()
{
	double err = pow(10,6);
	int iteration(0);

	while ( err > m_tol && iteration < m_nb_max_iterations )
	{
		One_Step_Smoothing();
		BoundarySlipping();

		err = L2_norm_relative_error(&m_old_coord_x, &m_old_coord_y);
		iteration++;
		Update_old_coords(&m_old_coord_x, &m_old_coord_y);
	}

	return AbstractSmoothLineSweeping_2D::SUCCESS;
}
/*------------------------------------------------------------------------*/




/*------------------------------------------------------------------------*/
void AbstractSmoothLineSweeping_2D::One_Step_Smoothing(){

	for (int i=1; i < m_Nx; i++)
	{
		for (int j=1; j < m_Ny; j++)
		{
			math::Point new_ideal_pos = ComputeNewPosition(i,j) ;

			// Damping
			Node n = (*m_B)(i, j);
			n.setPoint( m_theta*n.point() + (1.0-m_theta)*new_ideal_pos );
		}
	}

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AbstractSmoothLineSweeping_2D::BoundarySlipping()
{
	for (int i=1;i<m_Nx;i++)
	{
		math::Point Mid = FindMidBranche((*m_B)(i-1,0).point(), (*m_B)(i,0).point(), (*m_B)(i+1,0).point());
		(*m_B)(i, 0).setPoint(m_theta*(*m_B)(i, 0).point() + (1.0-m_theta)*Mid );

		Mid = FindMidBranche((*m_B)(i-1,m_Ny).point(), (*m_B)(i,m_Ny).point(), (*m_B)(i+1,m_Ny).point());
		(*m_B)(i, m_Ny).setPoint(m_theta*(*m_B)(i, m_Ny).point() + (1.0-m_theta)*Mid );
	}


	for (int j=1;j<m_Ny;j++)
	{
		math::Point Mid = FindMidBranche((*m_B)(0,j-1).point(), (*m_B)(0,j).point(), (*m_B)(0,j+1).point());
		(*m_B)(0, j).setPoint(m_theta*(*m_B)(0, j).point() + (1.0-m_theta)*Mid);

		Mid = FindMidBranche((*m_B)(m_Nx,j-1).point(), (*m_B)(m_Nx,j).point(), (*m_B)(m_Nx,j+1).point());
		(*m_B)(m_Nx, j).setPoint(m_theta*(*m_B)(m_Nx, j).point() + (1.0-m_theta)*Mid);
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

	for (int i=0; i < m_Nx; i++) {
		for (int j = 0; j < m_Ny; j++) {
			Node n = (*m_B)(i, j);

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