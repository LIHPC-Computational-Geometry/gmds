//
// Created by rochec on 18/01/24.
//
/*------------------------------------------------------------------------*/
#include <gmds/aero/ComputeBezierCtrlPtstoInterpolateSurface.h>
#include <gmds/math/BezierSurface.h>
#include <gmds/math/TransfiniteInterpolation.h>
#include <gmds/aero/Utils.h>
#include <gmds/ig/Mesh.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

ComputeBezierCtrlPtstoInterpolateSurface::ComputeBezierCtrlPtstoInterpolateSurface(cad::GeomSurface* ASurface,
                                                                                   Array2D<math::Point>* AInitCtrlPts) :
  m_Surface(ASurface),
  m_InitCtrlPts(AInitCtrlPts),
  m_degree_m(AInitCtrlPts->nbLines()-1),
  m_degree_n(AInitCtrlPts->nbColumns()-1)
{
	m_CtrlPts = new Array2D<math::Point>(m_degree_m+1,m_degree_n+1);
	m_InterpolatePts = new Array2D<math::Point>(m_degree_m+1,m_degree_n+1);
	for (int i=0;i<=m_degree_m;i++)
	{
		for (int j=0;j<=m_degree_n;j++)
		{
			(*m_CtrlPts)(i,j) = (*m_InitCtrlPts)(i,j);
			(*m_InterpolatePts)(i,j) = (*m_InitCtrlPts)(i,j);
		}
	}
}
/*------------------------------------------------------------------------*/
ComputeBezierCtrlPtstoInterpolateSurface::STATUS
ComputeBezierCtrlPtstoInterpolateSurface::execute()
{
	InitInterpolatePoints();
	ComputeControlPoints();

	return ComputeBezierCtrlPtstoInterpolateSurface::SUCCESS;
}
/*------------------------------------------------------------------------*/
Array2D<math::Point>
ComputeBezierCtrlPtstoInterpolateSurface::getCtrlPts()
{
	return (*m_CtrlPts);
}
/*------------------------------------------------------------------------*/
void
ComputeBezierCtrlPtstoInterpolateSurface::InitInterpolatePoints()
{
	// Init the grid of the control points with the four nodes of the face
	//(*m_InterpolatePts)(0,0) = f_nodes[0].point();
	//(*m_InterpolatePts)(m_degree_m,0) = f_nodes[1].point();
	//(*m_InterpolatePts)(m_degree_m,m_degree_n) = f_nodes[2].point();
	//(*m_InterpolatePts)(0,m_degree_n) = f_nodes[3].point();

	// We consider here that the four corner points are on the surface
	// Init the boundary nodes of the points to interpolate
	for (int i=1;i<m_degree_m;i++)
	{
		(*m_InterpolatePts)(i,0) = (*m_InterpolatePts)(0,0) + i*((*m_InterpolatePts)(m_degree_m,0)-(*m_InterpolatePts)(0,0))/(float(m_degree_m)) ;
		(*m_InterpolatePts)(i,m_degree_n) = (*m_InterpolatePts)(0,m_degree_n) + i*((*m_InterpolatePts)(m_degree_m,m_degree_n)-(*m_InterpolatePts)(0,m_degree_n))/(float(m_degree_m)) ;
	}
	for (int j=1;j<m_degree_n;j++)
	{
		(*m_InterpolatePts)(0,j) = (*m_InterpolatePts)(0,0) + j*((*m_InterpolatePts)(0,m_degree_n)-(*m_InterpolatePts)(0,0))/(float(m_degree_n)) ;
		(*m_InterpolatePts)(m_degree_n,j) = (*m_InterpolatePts)(m_degree_m,0) + j*((*m_InterpolatePts)(m_degree_m,m_degree_n)-(*m_InterpolatePts)(m_degree_m,0))/(float(m_degree_n)) ;
	}

	// Init the inner nodes of the grid by transfinite interpolation
	math::TransfiniteInterpolation::computeQuad((*m_InterpolatePts));

	// Project the nodes onto the surface to interpolate
	for (int i=0;i<=m_degree_m;i++)
	{
		for (int j=0;j<=m_degree_n;j++)
		{
			math::Point p = (*m_InterpolatePts)(i,j);
			m_Surface->project(p);
			(*m_InterpolatePts)(i,j) = p;
		}
	}

}
/*------------------------------------------------------------------------*/
void
ComputeBezierCtrlPtstoInterpolateSurface::ComputeControlPoints()
{
	Eigen::MatrixXd mat_B((m_degree_m+1)*(m_degree_n+1), (m_degree_m+1)*(m_degree_n+1));
	Eigen::VectorXd ctrl_points_x((m_degree_m+1)*(m_degree_n+1));
	Eigen::VectorXd ctrl_points_y((m_degree_m+1)*(m_degree_n+1));
	Eigen::VectorXd ctrl_points_z((m_degree_m+1)*(m_degree_n+1));
	Eigen::VectorXd interp_points_x((m_degree_m+1)*(m_degree_n+1));
	Eigen::VectorXd interp_points_y((m_degree_m+1)*(m_degree_n+1));
	Eigen::VectorXd interp_points_z((m_degree_m+1)*(m_degree_n+1));

	for (int i=0;i<=m_degree_m;i++)
	{
		for (int j=0;j<=m_degree_n;j++)
		{
			interp_points_x[i+j*(m_degree_m+1)] = (*m_InterpolatePts)(i,j).X();
			interp_points_y[i+j*(m_degree_m+1)] = (*m_InterpolatePts)(i,j).Y();
			interp_points_z[i+j*(m_degree_m+1)] = (*m_InterpolatePts)(i,j).Z();
		}
	}

	// Matrix Assembly
	for (int i=0;i<=m_degree_m;i++)
	{
		for (int j=0;j<=m_degree_n;j++)
		{
			double u = 1.0*i/m_degree_m ;
			double v = 1.0*j/m_degree_n ;
			for (int i2=0;i2<=m_degree_m;i2++)
			{
				for (int j2=0;j2<=m_degree_n;j2++)
				{
					double Bi2n = math::Utils::BernsteinPolynomial(m_degree_m, i2, u) ;
					double Bj2m = math::Utils::BernsteinPolynomial(m_degree_n, j2, v) ;
					mat_B(i+j*(m_degree_m+1), i2+j2*(m_degree_m+1)) = Bi2n*Bj2m ;
				}
			}
		}
	}

	// Solve the system
	Eigen::MatrixXd mat_B_inv = mat_B.inverse();
	ctrl_points_x = mat_B_inv*interp_points_x;
	ctrl_points_y = mat_B_inv*interp_points_y;
	ctrl_points_z = mat_B_inv*interp_points_z;

	// Put the positions of the control points in the array
	for (int i=0;i<=m_degree_m;i++)
	{
		for (int j=0;j<=m_degree_n;j++)
		{
			math::Point p;
			p.setX( ctrl_points_x(i+j*(m_degree_m+1)) );
			p.setY( ctrl_points_y(i+j*(m_degree_m+1)) );
			p.setZ( ctrl_points_z(i+j*(m_degree_m+1)) );
			(*m_CtrlPts)(i,j) = p;
		}
	}

}
/*------------------------------------------------------------------------*/