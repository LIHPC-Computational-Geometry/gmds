//
// Created by rochec on 13/02/24.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/ComputeBezierCurveCtrlPtstoInterpolateCurve.h>
#include <gmds/math/BezierSurface.h>
#include <gmds/math/TransfiniteInterpolation.h>
#include <gmds/claire/Utils.h>
#include <gmds/ig/Mesh.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

ComputeBezierCurveCtrlPtstoInterpolateCurve::ComputeBezierCurveCtrlPtstoInterpolateCurve(cad::GeomCurve* ACurve,
                                                                                         std::vector<math::Point>* AInitCtrlPts) :
  m_Curve(ACurve),
  m_InitCtrlPts(AInitCtrlPts),
  m_degree_n(AInitCtrlPts->size()-1)
{
	m_CtrlPts = new std::vector<math::Point>(m_degree_n+1);
	m_InterpolatePts = new std::vector<math::Point>(m_degree_n+1);
	for (int i=0;i<=m_degree_n;i++)
	{
		(*m_CtrlPts)[i] = (*m_InitCtrlPts)[i];
		(*m_InterpolatePts)[i] = (*m_InitCtrlPts)[i];
	}
}
/*------------------------------------------------------------------------*/
ComputeBezierCurveCtrlPtstoInterpolateCurve::STATUS
ComputeBezierCurveCtrlPtstoInterpolateCurve::execute()
{
	InitInterpolatePoints();
	ComputeControlPoints();

	return ComputeBezierCurveCtrlPtstoInterpolateCurve::SUCCESS;
}
/*------------------------------------------------------------------------*/
std::vector<math::Point>
ComputeBezierCurveCtrlPtstoInterpolateCurve::getCtrlPts()
{
	return (*m_CtrlPts);
}
/*------------------------------------------------------------------------*/
void
ComputeBezierCurveCtrlPtstoInterpolateCurve::InitInterpolatePoints()
{
	// We consider here that the two corner points are on the curve
	// Init the boundary nodes of the points to interpolate

	// Project the points onto the surface to interpolate
	for (int i=0;i<=m_degree_n;i++)
	{
		math::Point p = (*m_InterpolatePts)[i];
		m_Curve->project(p);
		(*m_InterpolatePts)[i] = p;
	}

}
/*------------------------------------------------------------------------*/
void
ComputeBezierCurveCtrlPtstoInterpolateCurve::ComputeControlPoints()
{
	Eigen::MatrixXd mat_B(m_degree_n+1, m_degree_n+1);
	Eigen::VectorXd ctrl_points_x(m_degree_n+1);
	Eigen::VectorXd ctrl_points_y(m_degree_n+1);
	Eigen::VectorXd ctrl_points_z(m_degree_n+1);
	Eigen::VectorXd interp_points_x(m_degree_n+1);
	Eigen::VectorXd interp_points_y(m_degree_n+1);
	Eigen::VectorXd interp_points_z(m_degree_n+1);

	for (int i=0;i<=m_degree_n;i++)
	{
		interp_points_x[i] = (*m_InterpolatePts)[i].X();
		interp_points_y[i] = (*m_InterpolatePts)[i].Y();
		interp_points_z[i] = (*m_InterpolatePts)[i].Z();
	}

	// Matrix Assembly
	for (int i=0;i<=m_degree_n;i++)
	{
		double u = 1.0*i/m_degree_n ;
			for (int i2=0;i2<=m_degree_n;i2++)
			{
			   double Bi2n = math::Utils::BernsteinPolynomial(m_degree_n, i2, u) ;
			   mat_B(i, i2) = Bi2n;
			}
	}

	// Solve the system
	Eigen::MatrixXd mat_B_inv = mat_B.inverse();
	ctrl_points_x = mat_B_inv*interp_points_x;
	ctrl_points_y = mat_B_inv*interp_points_y;
	ctrl_points_z = mat_B_inv*interp_points_z;

	// Put the positions of the control points in the array
	for (int i=0;i<=m_degree_n;i++)
	{
			math::Point p;
			p.setX( ctrl_points_x(i) );
			p.setY( ctrl_points_y(i) );
			p.setZ( ctrl_points_z(i) );
			(*m_CtrlPts)[i] = p;
	}

}
/*------------------------------------------------------------------------*/