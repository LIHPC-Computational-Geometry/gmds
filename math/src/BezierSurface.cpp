//
// Created by rochec on 28/10/2022.
//
/*----------------------------------------------------------------------------*/
#include <gmds/math/BezierSurface.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*--------------------------------------------------------------------------*/
namespace math{
/*------------------------------------------------------------------------*/
BezierSurface::BezierSurface(const Array2D<math::Point>& APts) :
  m_control_points(APts)
{
	m_degree_n = m_control_points.nbLines()-1;
	m_degree_m = m_control_points.nbColumns()-1;
}
/*------------------------------------------------------------------------*/
math::Point BezierSurface::operator()(const double& u, const double& v) const {

	Array2D<math::Point> current_control = m_control_points;
	math::Point M({0,0,0});

	for (int i=0; i<=m_degree_n; i++)
	{
		for (int j=0; j<=m_degree_m; j++)
		{
			double Bnu = BernsteinPolynomial(m_degree_n, i, u);
			double Bnv = BernsteinPolynomial(m_degree_m, j, v);
			M = M + Bnu*Bnv*m_control_points(i,j) ;
		}
	}

	return M;
}
/*------------------------------------------------------------------------*/
double BezierSurface::BinomialCoefficient(int n, int k) const
{
	double factorial_n(1);
	double factorial_k(1);
	double factorial_n_k(1);

	for (int i=2; i <= n; i++)
	{
		factorial_n = factorial_n*i;
		if (i <= k)
		{
			factorial_k = factorial_k*i;
		}
		if (i <= n-k)
		{
			factorial_n_k = factorial_n_k*i;
		}
	}

	return factorial_n/(factorial_k*factorial_n_k);

}
/*------------------------------------------------------------------------*/
double BezierSurface::BernsteinPolynomial(int n, int i, double u) const
{
	double C_n_i = BinomialCoefficient(n, i);
	return C_n_i* pow(u, i)* pow(1-u, n-i);
}
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/