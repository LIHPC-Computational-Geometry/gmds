//
// Created by rochec on 22/07/2022.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/SmoothLineSweepingOrtho.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

SmoothLineSweepingOrtho::SmoothLineSweepingOrtho(Blocking2D::Block *AB, int Anb_max_it, double Atheta) :
  AbstractSmoothLineSweeping_2D(AB, Anb_max_it, Atheta)
{

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
math::Point SmoothLineSweepingOrtho::ComputeNewPosition(int i, int j)
{
	double alpha = pow( (j-1.0)/( 6.0*(m_Ny-1.0) ) , 0.01) ;

	// Compute the 6 points for the Yao Smoother
	math::Point V1 = WeightedPointOnBranch((*m_B)(i-1,j-1).point(), (*m_B)(i,j-1).point(), (*m_B)(i+1,j-1).point(), 0.5);
	math::Point V2 = WeightedPointOnBranch((*m_B)(i-1,j).point(), (*m_B)(i,j).point(), (*m_B)(i+1,j).point(), 0.5);
	math::Point V3 = WeightedPointOnBranch((*m_B)(i-1,j+1).point(), (*m_B)(i,j+1).point(), (*m_B)(i+1,j+1).point(), 0.5);

	math::Point H1 = WeightedPointOnBranch((*m_B)(i-1,j-1).point(), (*m_B)(i-1,j).point(), (*m_B)(i-1,j+1).point(), 0.5);
	math::Point H2 = WeightedPointOnBranch((*m_B)(i,j-1).point(), (*m_B)(i,j).point(), (*m_B)(i,j+1).point(), 0.5);
	math::Point H3 = WeightedPointOnBranch((*m_B)(i+1,j-1).point(), (*m_B)(i+1,j).point(), (*m_B)(i+1,j+1).point(), 0.5);

	// Finding the intersection between the 4 segments
	math::Point X2;
	math::Segment Seg_Vert_1(V1, V2);
	math::Segment Seg_Vert_2(V2, V3);
	math::Segment Seg_Hori_1(H1, H2);
	math::Segment Seg_Hori_2(H2, H3);

	bool intersection_trouvee = Seg_Vert_1.intersect2D(Seg_Hori_1, X2);
	if (!intersection_trouvee) {
		intersection_trouvee = Seg_Vert_1.intersect2D(Seg_Hori_2, X2);
	}
	if (!intersection_trouvee) {
		intersection_trouvee = Seg_Vert_2.intersect2D(Seg_Hori_1, X2);
	}
	if (!intersection_trouvee) {
		intersection_trouvee = Seg_Vert_2.intersect2D(Seg_Hori_2, X2);
	}

	if (!intersection_trouvee) {
		X2 = (*m_B)(i,j).point();
	}

	math::Point X1 = ComputeOrtho(i, j);

	math::Point P_new;

	if ( Seg_Hori_1.isIn(X1) && Seg_Hori_1.isIn(X2) )
	{
		P_new = alpha*X2 + (1.0-alpha)*X1 ;
	}
	else if ( Seg_Hori_2.isIn(X1) && Seg_Hori_2.isIn(X2) )
	{
		P_new = alpha*X2 + (1.0-alpha)*X1 ;
	}
	else if ( Seg_Hori_1.isIn(X1) && Seg_Hori_2.isIn(X2) )
	{
		math::Vector3d v_s1 = H2-X1;
		math::Vector3d v_s2 = X2-H2;
		double l = alpha*( v_s1.norm() + v_s2.norm() );
		if ( l <= v_s1.norm() )
		{
			math::Vector3d u = H2-X1;
			u.normalize();
			P_new = X1 + l*u;
		}
		else
		{
			math::Vector3d u = H2-X2;
			u.normalize();
			P_new = X2 + ( v_s1.norm() + v_s2.norm() - l )*u;
		}
	}
	else if ( Seg_Hori_2.isIn(X1) && Seg_Hori_1.isIn(X2) )
	{
		math::Vector3d v_s1 = H2-X2;
		math::Vector3d v_s2 = X1-H2;
		double l = (1.0-alpha)*( v_s1.norm() + v_s2.norm() );
		if ( l <= v_s1.norm() )
		{
			math::Vector3d u = H2-X2;
			u.normalize();
			P_new = X2 + l*u;
		}
		else
		{
			math::Vector3d u = H2-X1;
			u.normalize();
			P_new = X1 + ( v_s1.norm() + v_s2.norm() - l )*u;
		}
	}
	else
	{
		P_new = X2 ;
	}

	return P_new;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
math::Point SmoothLineSweepingOrtho::ComputeOrtho(int i, int j)
{
	math::Point X1;

	math::Vector3d v1 = (*m_B)(i, 0).point() - (*m_B)(i-1, 0).point() ;
	math::Vector3d v2 = (*m_B)(i+1, 0).point() - (*m_B)(i, 0).point() ;

	math::Vector3d n1({v1.Y(), -v1.X(), 0});
	math::Vector3d n2({v2.Y(), -v2.X(), 0});

	n1.normalize();
	n2.normalize();

	math::Vector3d n = n1+n2;

	// Finding the intersection between the 4 segments
	math::Point H1 = WeightedPointOnBranch((*m_B)(i-1,j-1).point(), (*m_B)(i-1,j).point(), (*m_B)(i-1,j+1).point(), 0.5);
	math::Point H2 = WeightedPointOnBranch((*m_B)(i,j-1).point(), (*m_B)(i,j).point(), (*m_B)(i,j+1).point(), 0.5);
	math::Point H3 = WeightedPointOnBranch((*m_B)(i+1,j-1).point(), (*m_B)(i+1,j).point(), (*m_B)(i+1,j+1).point(), 0.5);

	math::Segment Seg_Hori_1(H1, H2);
	math::Segment Seg_Hori_2(H2, H3);

	//math::Segment Seg_Ortho((*m_B)(i,j-1).point()+pow(10,6)*(-n), (*m_B)(i,j-1).point()+pow(10,6)*n);
	math::Point Q = (*m_B)(i,j-1).point() ;
	math::Point P1_s, P2_s;
	P1_s.setXYZ( Q.X() - pow(10,1)*n.X(),  Q.Y() - pow(10,1)*n.Y(), 0);
	P2_s.setXYZ( Q.X() + pow(10,1)*n.X(),  Q.Y() + pow(10,1)*n.Y(), 0);
	math::Segment Seg_Ortho(P1_s, P2_s);

	bool intersection_trouvee = Seg_Ortho.intersect2D(Seg_Hori_1, X1);
	if (!intersection_trouvee) {
		intersection_trouvee = Seg_Ortho.intersect2D(Seg_Hori_2, X1);
	}

	if (!intersection_trouvee) {
		X1 = (*m_B)(i,j).point();
	}

	/*
	if (i==5 && j==5)
	{
		std::cout << v1 << std::endl;
		std::cout << v2 << std::endl;
		std::cout << " " << std::endl;
		std::cout << n1 << std::endl;
		std::cout << n2 << std::endl;
		std::cout << n << std::endl;
		std::cout << " " << std::endl;
		std::cout << H1 << std::endl;
		std::cout << H2 << std::endl;
		std::cout << H3 << std::endl;
		std::cout << " " << std::endl;
		std::cout << P1_s << std::endl;
		std::cout << P2_s << std::endl;
		std::cout << " " << std::endl;
		std::cout << "Intersection trouvée " << intersection_trouvee << std::endl;
		std::cout << (*m_B)(i,j).point() << std::endl;
		std::cout << X1 << std::endl;
	}
	 */

	return X1;
}
/*------------------------------------------------------------------------*/