//
// Created by rochec on 20/07/2022.
//

/*------------------------------------------------------------------------*/
#include<gmds/math/Line.h>
#include <gmds/claire/AbstractSmoothLineSweeping_2D.h>
#include <gmds/claire/Utils.h>
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

	for (int i=0; i <= m_Nx; i++)
	{
		for (int j=0; j <= m_Ny; j++)
		{
			Node n = (*m_B)(i, j);
			m_P_new(i,j) = n.point() ;
		}
	}
}
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

	//std::cout << "Nombre d'itérations lissage : " << iteration << std::endl;

	return AbstractSmoothLineSweeping_2D::SUCCESS;
}
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
void AbstractSmoothLineSweeping_2D::BoundarySlipping()
{

	for (int i=1;i<m_Nx;i++)
	{
		//math::Point Mid = math::Utils::WeightedPointOnBranch((*m_B)(i-1,0).point(), (*m_B)(i,0).point(), (*m_B)(i+1,0).point(), 0.5);
		//m_P_new(i,0) =( m_theta*(*m_B)(i, 0).point() + (1.0-m_theta)*Mid ) ;

		math::Point Mid = math::Utils::WeightedPointOnBranch((*m_B)(i-1,m_Ny).point(), (*m_B)(i,m_Ny).point(), (*m_B)(i+1,m_Ny).point(), 0.5);
		m_P_new(i,m_Ny) =( m_theta*(*m_B)(i, m_Ny).point() + (1.0-m_theta)*Mid ) ;

		/*
		math::Vector3d v = ( (*m_B)(i,m_Ny-1).point() - (*m_B)(i,m_Ny-2).point() ).normalize() ;
		math::Segment tan( (*m_B)(i,m_Ny-2).point(), (*m_B)(i,m_Ny-2).point() + pow(10,6)*v ) ;
		math::Segment b1( (*m_B)(i-1,m_Ny).point(), (*m_B)(i,m_Ny).point() );
		math::Segment b2( (*m_B)(i,m_Ny).point(), (*m_B)(i+1,m_Ny).point() );
		math::Point P_intersection;
		bool int_found = tan.intersect2D(b1, P_intersection);
		if (!int_found)
		{
			tan.intersect2D(b2, P_intersection);
		}
		if (!int_found)
		{
			P_intersection = math::Utils::WeightedPointOnBranch((*m_B)(i-1,m_Ny).point(), (*m_B)(i,m_Ny).point(), (*m_B)(i+1,m_Ny).point(), 0.5);
		}

		m_P_new(i,m_Ny) = ( m_theta*(*m_B)(i, m_Ny).point() + (1.0-m_theta)*P_intersection ) ;
		 */

	}

	for (int j=1;j<m_Ny;j++)
	{
		math::Point Mid = math::Utils::WeightedPointOnBranch((*m_B)(0,j-1).point(), (*m_B)(0,j).point(), (*m_B)(0,j+1).point(), 0.5);
		m_P_new(0,j) = ( m_theta*(*m_B)(0, j).point() + (1.0-m_theta)*Mid ) ;

		Mid = math::Utils::WeightedPointOnBranch((*m_B)(m_Nx,j-1).point(), (*m_B)(m_Nx,j).point(), (*m_B)(m_Nx,j+1).point(), 0.5);
		m_P_new(m_Nx,j) = ( m_theta*(*m_B)(m_Nx, j).point() + (1.0-m_theta)*Mid ) ;
	}

}
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