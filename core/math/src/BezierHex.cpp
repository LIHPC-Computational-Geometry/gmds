//
// Created by rochec on 05/12/23.
//
/*----------------------------------------------------------------------------*/
#include <gmds/math/BezierHex.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*--------------------------------------------------------------------------*/
namespace math{
/*------------------------------------------------------------------------*/
BezierHex::BezierHex(const Array3D<math::Point>& APts) :
  m_control_points(APts)
{
	m_degree_l = m_control_points.nbElements(0)-1;
	m_degree_m = m_control_points.nbElements(1)-1;
	m_degree_n = m_control_points.nbElements(2)-1;
}
/*------------------------------------------------------------------------*/
math::Point BezierHex::operator()(const double& Au, const double& Av, const double& Aw) const {

	math::Point M(0.,0.,0.);

	for (int i=0; i<=m_degree_l; i++)
	{
		for (int j=0; j<=m_degree_m; j++)
		{
			for (int k=0; k<=m_degree_n; k++)
			{
				double Blu = BernsteinPolynomial(m_degree_l, i, Au);
				double Bmv = BernsteinPolynomial(m_degree_m, j, Av);
				double Bnw = BernsteinPolynomial(m_degree_n, k, Aw);
				M = M + Blu * Bmv * Bnw * m_control_points(i, j, k);
			}
		}
	}

	return M;
}
/*------------------------------------------------------------------------*/
Array3D<math::Point> BezierHex:: getDiscretization(const int ANb_l, const int ANb_m, const int ANb_n) const {
	if(ANb_l<1 || ANb_m<1 || ANb_n<1)
		throw GMDSException("BezierHex discretization impossible with this parameters");

	Array3D<math::Point> points(ANb_l+1,ANb_m+1, ANb_n+1);
	double step_l = 1.0/ANb_l;
	double step_m = 1.0/ANb_m;
	double step_n = 1.0/ANb_n;
	for(int i=0; i<=ANb_l; i++)
	{
		for (int j=0; j<=ANb_m; j++)
		{
			for (int k=0; k<=ANb_n; k++)
			{
				points(i, j, k) = this->operator()(i * step_l, j * step_m, k*step_n);
			}
		}
	}
	return points;
}
/*----------------------------------------------------------------------------*/
double BezierHex::BinomialCoefficient(const int An, const int Ak) const
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
double BezierHex::BernsteinPolynomial(const int An, const int Ai, const double Au) const
{
	double C_n_i = BinomialCoefficient(An, Ai);
	return C_n_i* pow(Au, Ai)* pow(1-Au, An-Ai);
}
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
      /*----------------------------------------------------------------------------*/