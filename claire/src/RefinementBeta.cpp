//
// Created by rochec on 27/07/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/RefinementBeta.h>
//#include <math.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

RefinementBeta::RefinementBeta(std::vector<math::Point> APoints, double Asize_first_edge)
{
	m_Points = APoints;
	m_size_first_edge = Asize_first_edge;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
RefinementBeta::STATUS
RefinementBeta::execute()
{
	int Nbr_Points = m_Points.size() ;
	double sum_size_edges = 0.0;
	std::vector<double> longueurs;
	longueurs.push_back(0.0);		// The first node is at the lenght 0

	for (int i=0; i < Nbr_Points-1; i++)
	{
		math::Vector3d v = m_Points[i+1]-m_Points[i] ;
		sum_size_edges += v.norm() ;
		longueurs.push_back(sum_size_edges);
	}

	double Beta = ComputeBeta(m_size_first_edge, sum_size_edges, Nbr_Points) ;

	m_PointsRefined.push_back(m_Points[0]);	// The first position is fixed

	// Compute the position of the new points
	for (int i=1; i < Nbr_Points-1; i++)
	{
		double n = (i-1.0)/(Nbr_Points-1.0);
		double r = (Beta+1.0)/(Beta-1.0);
		double z = log(r);
		double p = z*(1.0-n);
		double fn = 1.0+Beta*(1.0-exp(p))/(1.0+exp(p));
		double l = fn*sum_size_edges;

		for (int j=1; j < Nbr_Points-1; j++)
		{
			if ( l >= longueurs[j] && l <= longueurs[j+1] )
			{
				l = l - longueurs[j];
				math::Vector3d v = ( m_Points[j+1]-m_Points[j] ).normalize() ;
				m_PointsRefined.push_back( m_Points[j] + l*v );
			}
		}

	}

	m_PointsRefined.push_back(m_Points[Nbr_Points-1]);		// The last position is fixed

	return RefinementBeta::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double RefinementBeta::ComputeBeta(double first_edge_size, double sum_edge_sizes, int Nbr_Points)
{

	double n = 1.0/(Nbr_Points-1.0);
	double tol = pow(10,-6); 	// Tolerance for the Newton algorithm
	double Beta = 1.000000001;			// Initialization of the Beta parameter
	double F = 1.0;						// Initialization of the function F
	int iter = 0;							// Initialization iterations of the Newton

	while ( abs(F) >= tol && iter < 10000 )
	{
		double r = (Beta+1.0)/(Beta-1.0) ;
		double z = log(r);
		double p = z*(1.0-n);
		double fn = 1.0+Beta*(1.0-exp(p))/(1.0+exp(p));

		F = - first_edge_size + fn*sum_edge_sizes; 	// Compute F

		// Compute F'
		double u = pow(Beta-1.0, 1.0-n) - pow(Beta+1.0, 1.0-n) ;
		double up = (1.0-n)*pow(Beta-1.0, -n) - (1.0-n)*pow(Beta+1.0, -n) ;
		double v = pow(Beta-1.0, 1.0-n) + pow(Beta+1.0, 1.0-n) ;
		double vp = (1.0-n)*pow(Beta-1.0, -n) + (1.0-n)*pow(Beta+1.0, -n) ;
		double Fp = u/v + Beta*(up*v-u*vp)/(pow(v,2.0));
		Fp = Fp*sum_edge_sizes;

		if (Fp != 0.0)
		{
			Beta = Beta - F/Fp;
			iter ++;
		}
		else
		{
			std::cout << "ERREUR : Division par 0 dans RefinementBeta." << std::endl;
			F = 0;
			iter = pow(10,8);
		}

	}

	return Beta;

}
/*------------------------------------------------------------------------*/