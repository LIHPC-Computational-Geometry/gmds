//
// Created by rochec on 14/04/2022.
//

/*----------------------------------------------------------------------------*/
#include <gmds/claire/AeroMeshQuality.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace math {


/*------------------------------------------------------------------------*/
double AeroMeshQuality::minlenghtedge(math::Point p0, math::Point p1, math::Point p2, math::Point p3){

	Vector3d v1 = p1-p0;
	Vector3d v2 = p1-p2;
	Vector3d v3 = p2-p3;
	Vector3d v4 = p3-p0;

	return std::min(std::min(v1.norm(), v2.norm()), std::min(v3.norm(), v4.norm()));
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double AeroMeshQuality::AngleOuverture(const math::Point& p0, math::Point p1, math::Point p2, math::Point p3)
{
	Vector3d v1 = p1-p0;
	Vector3d v2 = p2-p0;
	Vector3d v3 = p3-p0;
	v1.normalize();
	v2.normalize();
	v3.normalize();

	double angle1 = acos(v1.dot(v2));
	double angle2 = acos(v1.dot(v3));

	return angle1+angle2;

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double AeroMeshQuality::AspectRatioQUAD(math::Point p0, math::Point p1, math::Point p2, math::Point p3){

	Vector3d a = p2-p0;
	Vector3d b = p1-p3;

	return std::min(a.norm(), b.norm())/std::max(a.norm(), b.norm());

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double AeroMeshQuality::InternalAngleDeviationQUAD(math::Point p0, math::Point p1, math::Point p2, math::Point p3){

	Vector3d e1 = p1-p0;
	Vector3d e2 = p2-p1;
	Vector3d e3 = p3-p2;
	Vector3d e4 = p0-p3;

	e1.normalize();
	e2.normalize();
   e3.normalize();
	e4.normalize();

	double teta1 = acos(e1.dot(e2)) * 180.0 / M_PI;
	double teta2 = acos(e2.dot(e3)) * 180.0 / M_PI;
	double teta3 = acos(e3.dot(e4)) * 180.0 / M_PI;
	double teta4 = acos(e4.dot(e1)) * 180.0 / M_PI;

	double teta_min = std::min( std::min(teta1, teta2), std::min(teta3, teta4));
   double teta_max = std::max( std::max(teta1, teta2), std::max(teta3, teta4));

	return std::max( abs(90.0-teta_min), abs(teta_max-90.0) ) ;

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double AeroMeshQuality::EquiAngleSkewnessQUAD(math::Point p0, math::Point p1, math::Point p2, math::Point p3)
{
	Vector3d e1 = p1-p0;
	Vector3d e2 = p2-p1;
	Vector3d e3 = p3-p2;
	Vector3d e4 = p0-p3;

	e1.normalize();
	e2.normalize();
	e3.normalize();
	e4.normalize();

	double teta1 = acos(e1.dot(e2)) * 180.0 / M_PI;
	double teta2 = acos(e2.dot(e3)) * 180.0 / M_PI;
	double teta3 = acos(e3.dot(e4)) * 180.0 / M_PI;
	double teta4 = acos(e4.dot(e1)) * 180.0 / M_PI;

	double Q_min = std::min( std::min(teta1, teta2), std::min(teta3, teta4));
	double Q_max = std::max( std::max(teta1, teta2), std::max(teta3, teta4));
	double Q_eq = 90.0;

	return std::max( (Q_max-Q_eq)/(180.0-Q_eq), (Q_eq-Q_min)/Q_eq ) ;

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double AeroMeshQuality::ConditionQUAD(math::Point p0, math::Point p1, math::Point p2, math::Point p3)
{
	Vector3d L0 = p1-p0;
	Vector3d L1 = p2-p1;
	Vector3d L2 = p3-p2;
	Vector3d L3 = p0-p3;

	double l0 = L0.norm() ;
	double l1 = L1.norm() ;
	double l2 = L2.norm() ;
	double l3 = L3.norm() ;

	Vector3d N0 = L3.cross(L0) ;
	Vector3d N1 = L0.cross(L1) ;
	Vector3d N2 = L1.cross(L2) ;
	Vector3d N3 = L2.cross(L3) ;

	Vector3d X1 = (p1-p0) + (p2-p3) ;
	Vector3d X2 = (p2-p1) + (p3-p0) ;

	Vector3d Nc = X1.cross(X2) ;
	Vector3d nc = Nc /Nc.norm() ;

	double a0 = nc.dot(N0) ;
	double a1 = nc.dot(N1) ;
	double a2 = nc.dot(N2) ;
	double a3 = nc.dot(N3) ;

	double max1 = std::max( (pow(l0,2) + pow(l3, 2))/a0, (pow(l0,2) + pow(l1, 2))/a1 ) ;
	double max2 = std::max( (pow(l1,2) + pow(l2, 2))/a2, (pow(l2,2) + pow(l3, 2))/a3 ) ;

	return (1.0/2.0)*std::max( max1, max2  );
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double AeroMeshQuality::EdgeRatioQUAD(math::Point p0, math::Point p1, math::Point p2, math::Point p3)
{
	Vector3d L0 = p1-p0;
	Vector3d L1 = p2-p1;
	Vector3d L2 = p3-p2;
	Vector3d L3 = p0-p3;

	double l0 = L0.norm() ;
	double l1 = L1.norm() ;
	double l2 = L2.norm() ;
	double l3 = L3.norm() ;

	double Lmax = std::max( std::max(l0, l1), std::max(l2, l3));
	double Lmin = std::min( std::min(l0, l1), std::min(l2, l3));

	return Lmax/Lmin;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double AeroMeshQuality::JacobianQUAD(math::Point p0, math::Point p1, math::Point p2, math::Point p3)
{
	Vector3d L0 = p1-p0;
	Vector3d L1 = p2-p1;
	Vector3d L2 = p3-p2;
	Vector3d L3 = p0-p3;

	Vector3d N0 = L3.cross(L0) ;
	Vector3d N1 = L0.cross(L1) ;
	Vector3d N2 = L1.cross(L2) ;
	Vector3d N3 = L2.cross(L3) ;

	Vector3d X1 = (p1-p0) + (p2-p3) ;
	Vector3d X2 = (p2-p1) + (p3-p0) ;

	Vector3d Nc = X1.cross(X2) ;
	Vector3d nc = Nc / Nc.norm() ;

	double a0 = nc.dot(N0) ;
	double a1 = nc.dot(N1) ;
	double a2 = nc.dot(N2) ;
	double a3 = nc.dot(N3) ;

	return std::min( std::min(a0, a1), std::min(a2, a3) );
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double AeroMeshQuality::ScaledJacobianQUAD(math::Point p0, math::Point p1, math::Point p2, math::Point p3)
{
	Vector3d L0 = p1-p0;
	Vector3d L1 = p2-p1;
	Vector3d L2 = p3-p2;
	Vector3d L3 = p0-p3;

	double l0 = L0.norm() ;
	double l1 = L1.norm() ;
	double l2 = L2.norm() ;
	double l3 = L3.norm() ;

	Vector3d N0 = L3.cross(L0) ;
	Vector3d N1 = L0.cross(L1) ;
	Vector3d N2 = L1.cross(L2) ;
	Vector3d N3 = L2.cross(L3) ;

	Vector3d X1 = (p1-p0) + (p2-p3) ;
	Vector3d X2 = (p2-p1) + (p3-p0) ;

	Vector3d Nc = X1.cross(X2) ;
	Vector3d nc = Nc /Nc.norm() ;

	double a0 = nc.dot(N0) ;
	double a1 = nc.dot(N1) ;
	double a2 = nc.dot(N2) ;
	double a3 = nc.dot(N3) ;

	double sj0 = a0/(l0*l3);
	double sj1 = a1/(l0*l1);
	double sj2 = a2/(l1*l2);
	double sj3 = a3/(l2*l3);

	return std::min( std::min(sj0, sj1), std::min(sj2, sj3) );
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double AeroMeshQuality::ShapeQUAD(math::Point p0, math::Point p1, math::Point p2, math::Point p3)
{
	Vector3d L0 = p1-p0;
	Vector3d L1 = p2-p1;
	Vector3d L2 = p3-p2;
	Vector3d L3 = p0-p3;

	double l0 = L0.norm() ;
	double l1 = L1.norm() ;
	double l2 = L2.norm() ;
	double l3 = L3.norm() ;

	Vector3d N0 = L3.cross(L0) ;
	Vector3d N1 = L0.cross(L1) ;
	Vector3d N2 = L1.cross(L2) ;
	Vector3d N3 = L2.cross(L3) ;

	Vector3d X1 = (p1-p0) + (p2-p3) ;
	Vector3d X2 = (p2-p1) + (p3-p0) ;

	Vector3d Nc = X1.cross(X2) ;
	Vector3d nc = Nc /Nc.norm() ;

	double a0 = nc.dot(N0) ;
	double a1 = nc.dot(N1) ;
	double a2 = nc.dot(N2) ;
	double a3 = nc.dot(N3) ;

	double s0 = a0/(pow(l0,2)+pow(l3,2));
	double s1 = a1/(pow(l0,2)+pow(l1,2));
	double s2 = a2/(pow(l1,2)+pow(l2,2));
	double s3 = a3/(pow(l2,2)+pow(l3,2));

	return 2.0*std::min( std::min(s0, s1), std::min(s2, s3) );
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double AeroMeshQuality::SkewQUAD(math::Point p0, math::Point p1, math::Point p2, math::Point p3)
{
	Vector3d X1 = (p1-p0) + (p2-p3) ;
	Vector3d X2 = (p2-p1) + (p3-p0) ;

	Vector3d x1 = X1 / X1.norm() ;
	Vector3d x2 = X2 / X2.norm() ;

	return abs(x1.dot(x2));
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double AeroMeshQuality::StretchQUAD(math::Point p0, math::Point p1, math::Point p2, math::Point p3)
{
	Vector3d L0 = p1-p0;
	Vector3d L1 = p2-p1;
	Vector3d L2 = p3-p2;
	Vector3d L3 = p0-p3;

	double l0 = L0.norm() ;
	double l1 = L1.norm() ;
	double l2 = L2.norm() ;
	double l3 = L3.norm() ;

	Vector3d D0 = p2-p0;
	Vector3d D1 = p3-p1;

	double Dmax = std::max(D0.norm(), D1.norm());
	double lmin = std::min(std::min(l0, l1), std::min(l2, l3));

	return sqrt(2.0)*lmin/Dmax ;
}
/*------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
}  // namespace math
/*----------------------------------------------------------------------------*/
}  // namespace gmds
      /*----------------------------------------------------------------------------*/