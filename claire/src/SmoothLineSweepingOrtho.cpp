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
	double alpha = pow( (j-1.0)/( 6.0*(m_Ny-2.0) ) , 0.1) ;

	// Compute the 6 points for the Yao Smoother
	math::Point V1 = WeightedPointOnBranch((*m_B)(i-1,j-1).point(), (*m_B)(i,j-1).point(), (*m_B)(i+1,j-1).point(), 0.5);
	math::Point V2 = WeightedPointOnBranch((*m_B)(i-1,j).point(), (*m_B)(i,j).point(), (*m_B)(i+1,j).point(), 0.5);
	math::Point V3 = WeightedPointOnBranch((*m_B)(i-1,j+1).point(), (*m_B)(i,j+1).point(), (*m_B)(i+1,j+1).point(), 0.5);

	math::Point H1 = WeightedPointOnBranch((*m_B)(i-1,j-1).point(), (*m_B)(i-1,j).point(), (*m_B)(i-1,j+1).point(), 0.5);
	math::Point H2 = WeightedPointOnBranch((*m_B)(i,j-1).point(), (*m_B)(i,j).point(), (*m_B)(i,j+1).point(), 0.5);
	math::Point H3 = WeightedPointOnBranch((*m_B)(i+1,j-1).point(), (*m_B)(i+1,j).point(), (*m_B)(i+1,j+1).point(), 0.5);

	// Finding the intersection between the 4 segments
	bool intersection_trouvee(false);
	math::Point X2;
	math::Segment Seg_Vert_1(V1, V2);
	math::Segment Seg_Vert_2(V2, V3);
	math::Segment Seg_Hori_1(H1, H2);
	math::Segment Seg_Hori_2(H2, H3);

	intersection_trouvee = Seg_Vert_1.intersect2D(Seg_Hori_1, X2);
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

	math::Vector3d v1 = (*m_B)(i, j-1).point() - (*m_B)(i-1, j-1).point() ;
	math::Vector3d v2 = (*m_B)(i, j-1).point() - (*m_B)(i+1, j-1).point() ;

	math::Vector3d n1 = v1.getOneOrtho();
	math::Vector3d n2 = v2.getOneOrtho();

	math::Vector3d n = n1+n2;

	// Finding the intersection between the 4 segments
	math::Point H1 = WeightedPointOnBranch((*m_B)(i-1,j-1).point(), (*m_B)(i-1,j).point(), (*m_B)(i-1,j+1).point(), 0.5);
	math::Point H2 = WeightedPointOnBranch((*m_B)(i,j-1).point(), (*m_B)(i,j).point(), (*m_B)(i,j+1).point(), 0.5);
	math::Point H3 = WeightedPointOnBranch((*m_B)(i+1,j-1).point(), (*m_B)(i+1,j).point(), (*m_B)(i+1,j+1).point(), 0.5);

	math::Segment Seg_Hori_1(H1, H2);
	math::Segment Seg_Hori_2(H2, H3);

	math::Segment Seg_Ortho((*m_B)(i,j-1).point()+pow(10,6)*(-n), (*m_B)(i,j-1).point()+pow(10,6)*n);

	bool intersection_trouvee(false);
	if (!intersection_trouvee) {
		intersection_trouvee = Seg_Ortho.intersect2D(Seg_Hori_1, X1);
	}
	if (!intersection_trouvee) {
		intersection_trouvee = Seg_Ortho.intersect2D(Seg_Hori_2, X1);
	}

	if (!intersection_trouvee) {
		X1 = (*m_B)(i,j).point();
	}

	return X1;
}
/*------------------------------------------------------------------------*/