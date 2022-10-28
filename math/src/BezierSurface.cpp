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
	m_degree_m = m_control_points.nbLines()-1;
	m_degree_n = m_control_points.nbColumns()-1;
}
/*------------------------------------------------------------------------*/
math::Point BezierSurface::operator()(const double& Au, const double& Av) const {

	Array2D<math::Point> current_control = m_control_points;
	math::Point M({0,0,0});

	for (int i=0; i<=m_degree_m; i++)
	{
		for (int j=0; j<=m_degree_n; j++)
		{
			double Bmu = BernsteinPolynomial(m_degree_m, i, Au);
			double Bnv = BernsteinPolynomial(m_degree_n, j, Av);
			M = M + Bmu*Bnv*m_control_points(i,j) ;
		}
	}

	return M;
}
/*------------------------------------------------------------------------*/
Array2D<math::Point> BezierSurface:: getDiscretization(const int ANb_m, const int ANb_n) const {
	if(ANb_m<1 || ANb_n<1)
		throw GMDSException("BezierSurface discretization impossible with this parameters");

	Array2D<math::Point> points(ANb_m+1, ANb_n+1);
	double step_m = 1.0/ANb_m;
	double step_n = 1.0/ANb_n;
	for(int i=0; i<=ANb_m; i++)
	{
		for (int j=0; j<=ANb_n; j++)
		{
			points(i,j) = this->operator()(i*step_m, j*step_n);
		}
	}
	return points;
}
/*----------------------------------------------------------------------------*/
double BezierSurface::BinomialCoefficient(const int An, const int Ak) const
{
	double factorial_n(1);
	double factorial_k(1);
	double factorial_n_k(1);

	for (int i=2; i <= An; i++)
	{
		factorial_n = factorial_n*i;
		if (i <= Ak)
		{
			factorial_k = factorial_k*i;
		}
		if (i <= An-Ak)
		{
			factorial_n_k = factorial_n_k*i;
		}
	}

	return factorial_n/(factorial_k*factorial_n_k);

}
/*------------------------------------------------------------------------*/
double BezierSurface::BernsteinPolynomial(const int An, const int Ai, const double Au) const
{
	double C_n_i = BinomialCoefficient(An, Ai);
	return C_n_i* pow(Au, Ai)* pow(1-Au, An-Ai);
}
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/