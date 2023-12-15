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
double AeroMeshQuality::minlenghtedge(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3){

	Vector3d v1 = p1-p0;
	Vector3d v2 = p1-p2;
	Vector3d v3 = p2-p3;
	Vector3d v4 = p3-p0;

	return std::min(std::min(v1.norm(), v2.norm()), std::min(v3.norm(), v4.norm()));
}
/*------------------------------------------------------------------------*/
double AeroMeshQuality::AngleOuverture(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3)
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
double AeroMeshQuality::AspectRatioQUAD(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3){

	Vector3d a = p2-p0;
	Vector3d b = p1-p3;

	return std::min(a.norm(), b.norm())/std::max(a.norm(), b.norm());

}
/*------------------------------------------------------------------------*/
double AeroMeshQuality::InternalAngleDeviationQUAD(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3){

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
double AeroMeshQuality::EquiAngleSkewnessQUAD(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3)
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
double AeroMeshQuality::ConditionQUAD(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3)
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
double AeroMeshQuality::EdgeRatioQUAD(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3)
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
double AeroMeshQuality::JacobianQUAD(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3)
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
double AeroMeshQuality::ScaledJacobianQUAD(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3)
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
double AeroMeshQuality::ShapeQUAD(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3)
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
double AeroMeshQuality::SkewQUAD(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3)
{
	Vector3d X1 = (p1-p0) + (p2-p3) ;
	Vector3d X2 = (p2-p1) + (p3-p0) ;

	Vector3d x1 = X1 / X1.norm() ;
	Vector3d x2 = X2 / X2.norm() ;

	return abs(x1.dot(x2));
}
/*------------------------------------------------------------------------*/
double AeroMeshQuality::StretchQUAD(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3)
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
double
AeroMeshQuality::ScaledJacobianHEX(const math::Point& p0, const math::Point& p1, const math::Point& p2, const math::Point& p3,
                  const math::Point& p4, const math::Point& p5, const math::Point& p6, const math::Point& p7)
{
	Vector3d L0 = (p1-p0).normalize() ;
	Vector3d L1 = (p2-p1).normalize() ;
	Vector3d L2 = (p3-p2).normalize() ;
	Vector3d L3 = (p3-p0).normalize() ;
	Vector3d L4 = (p4-p0).normalize() ;
	Vector3d L5 = (p5-p1).normalize() ;
	Vector3d L6 = (p6-p2).normalize() ;
	Vector3d L7 = (p7-p3).normalize() ;
	Vector3d L8 = (p5-p4).normalize() ;
	Vector3d L9 = (p6-p5).normalize() ;
	Vector3d L10 = (p7-p6).normalize() ;
	Vector3d L11 = (p7-p4).normalize() ;

	Vector3d X1 = ((p1-p0) + (p2-p3) + (p5-p4) + (p6-p7)).normalize() ;
	Vector3d X2 = ((p3-p0) + (p2-p1) + (p7-p4) + (p6-p5)).normalize() ;
	Vector3d X3 = ((p4-p0) + (p5-p1) + (p6-p2) + (p7-p3)).normalize() ;

	Eigen::Matrix3d A0;
	Eigen::Matrix3d A1;
	Eigen::Matrix3d A2;
	Eigen::Matrix3d A3;
	Eigen::Matrix3d A4;
	Eigen::Matrix3d A5;
	Eigen::Matrix3d A6;
	Eigen::Matrix3d A7;
	Eigen::Matrix3d A8;

	// Fill the 3x3 Jacobian matrices
	for (int i=0;i<3;i++)
	{
		A0(i,0) = L0[i];
		A0(i,1) = L3[i];
		A0(i,2) = L4[i];
		A1(i,0) = L1[i];
		A1(i,1) = -L0[i];
		A1(i,2) = L5[i];
		A2(i,0) = L2[i];
		A2(i,1) = -L1[i];
		A2(i,2) = L6[i];
		A3(i,0) = -L3[i];
		A3(i,1) = -L2[i];
		A3(i,2) = L7[i];
		A4(i,0) = L11[i];
		A4(i,1) = L8[i];
		A4(i,2) = -L4[i];
		A5(i,0) = -L8[i];
		A5(i,1) = L9[i];
		A5(i,2) = -L5[i];
		A6(i,0) = -L9[i];
		A6(i,1) = L10[i];
		A6(i,2) = -L6[i];
		A7(i,0) = -L10[i];
		A7(i,1) = -L11[i];
		A7(i,2) = -L7[i];
		A8(i,0) = X1[i];
		A8(i,1) = X2[i];
		A8(i,2) = X3[i];
	}

	return std::min({A0.determinant(), A1.determinant(), A2.determinant(), A3.determinant(), A4.determinant(),
	                       A5.determinant(), A6.determinant(), A7.determinant(), A8.determinant()}) ;

}
/*------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
}  // namespace math
/*----------------------------------------------------------------------------*/
}  // namespace gmds
      /*----------------------------------------------------------------------------*/