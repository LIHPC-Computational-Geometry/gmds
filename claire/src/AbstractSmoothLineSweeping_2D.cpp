//
// Created by rochec on 20/07/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AbstractSmoothLineSweeping_2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AbstractSmoothLineSweeping_2D::AbstractSmoothLineSweeping_2D(Blocking2D::Block* AB, int Anb_max_it, double Atheta):
  	m_B(AB),
  	m_nb_max_iterations(Anb_max_it),
  	m_theta(Atheta),
  	m_Nx(m_B->getNbDiscretizationI()-1),
	m_Ny(m_B->getNbDiscretizationJ()-1),
  	m_P_new(m_Nx+1, m_Ny+1)
{
	m_tol = pow(10,-3);
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
AbstractSmoothLineSweeping_2D::STATUS AbstractSmoothLineSweeping_2D::execute()
{
	double err = pow(10,6);
	int iteration(0);

	// Block points initialization. The only points fixed during the whole smoothing.
	m_P_new(0,0) = (*m_B)(0,0).point() ;
	m_P_new(m_Nx,0) = (*m_B)(m_Nx,0).point() ;
	m_P_new(0,m_Ny) = (*m_B)(0,m_Ny).point() ;
	m_P_new(m_Nx,m_Ny) = (*m_B)(m_Nx,m_Ny).point() ;

	while ( err > m_tol && iteration < m_nb_max_iterations )
	{
		One_Step_Smoothing();
		BoundarySlipping();

		err = L2_norm_relative_error();
		Update_new_coords();
		iteration++;
	}

	std::cout << "Nombre d'itérations lissage : " << iteration << std::endl;

	return AbstractSmoothLineSweeping_2D::SUCCESS;
}
/*------------------------------------------------------------------------*/




/*------------------------------------------------------------------------*/
void AbstractSmoothLineSweeping_2D::One_Step_Smoothing(){

	for (int i=1; i < m_Nx; i++)
	{
		for (int j=1; j < m_Ny; j++)
		{
			math::Point P_new_ideal = ComputeNewPosition(i,j) ;

			// Damping
			Node n = (*m_B)(i, j);
			m_P_new(i,j) = m_theta*n.point() + (1.0-m_theta)*P_new_ideal;
		}
	}

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AbstractSmoothLineSweeping_2D::BoundarySlipping()
{
	for (int i=1;i<m_Nx;i++)
	{
		math::Point Mid = WeightedPointOnBranch((*m_B)(i-1,0).point(), (*m_B)(i,0).point(), (*m_B)(i+1,0).point(), 0.5);
		m_P_new(i,0) =( m_theta*(*m_B)(i, 0).point() + (1.0-m_theta)*Mid ) ;

		Mid = WeightedPointOnBranch((*m_B)(i-1,m_Ny).point(), (*m_B)(i,m_Ny).point(), (*m_B)(i+1,m_Ny).point(), 0.5);
		m_P_new(i,m_Ny) = ( m_theta*(*m_B)(i, m_Ny).point() + (1.0-m_theta)*Mid ) ;
	}

	for (int j=1;j<m_Ny;j++)
	{
		math::Point Mid = WeightedPointOnBranch((*m_B)(0,j-1).point(), (*m_B)(0,j).point(), (*m_B)(0,j+1).point(), 0.5);
		m_P_new(0,j) = ( m_theta*(*m_B)(0, j).point() + (1.0-m_theta)*Mid ) ;

		Mid = WeightedPointOnBranch((*m_B)(m_Nx,j-1).point(), (*m_B)(m_Nx,j).point(), (*m_B)(m_Nx,j+1).point(), 0.5);
		m_P_new(m_Nx,j) = ( m_theta*(*m_B)(m_Nx, j).point() + (1.0-m_theta)*Mid ) ;
	}

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
math::Point AbstractSmoothLineSweeping_2D::WeightedPointOnBranch(const math::Point A, const math::Point B, const math::Point C, double alpha) {
	math::Point P_Weighted;
	math::Vector3d Vec_AB = B-A ;
	math::Vector3d Vec_BC = C-B ;
	double norme_1 = Vec_AB.norm() ;
	double norme_2 = Vec_BC.norm() ;
	double norme_branche = norme_1 + norme_2 ;
	double norme_cible = alpha*norme_branche ;

	if (norme_cible <= norme_1){
		Vec_AB.normalize();
		P_Weighted = A + norme_cible*Vec_AB ;
	}
	else if (norme_cible > norme_1){
		math::Vector3d Vec_CB = - Vec_BC ;
		Vec_CB.normalize();
		P_Weighted = C + (norme_branche-norme_cible)*Vec_CB ;
	}
	else
	{
		P_Weighted = B ;
	}

	return P_Weighted;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double AbstractSmoothLineSweeping_2D::L2_norm_relative_error()
{
	double num(0.0);
	double denom(0.0);

	for (int i=0; i < m_Nx; i++) {
		for (int j = 0; j < m_Ny; j++)
		{
			Node n = (*m_B)(i, j);

			double err_loc = sqrt( pow(m_P_new(i,j).X() - n.X() , 2) + pow(m_P_new(i,j).Y() - n.Y(),2) );

			num += pow(err_loc, 2);
			denom += pow(m_P_new(i,j).X(), 2) + pow(m_P_new(i,j).Y(), 2) ;
		}
	}

	return sqrt(num/denom);
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AbstractSmoothLineSweeping_2D::Update_new_coords()
{
	for (int i=0; i <= m_Nx; i++)
	{
		for (int j=0; j <= m_Ny; j++)
		{
			Node n = (*m_B)(i, j);
			n.setPoint(m_P_new(i,j));
		}
	}
}
/*------------------------------------------------------------------------*/