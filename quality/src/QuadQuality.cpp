/*---------------------------------------------------------------------------*/
#include <gmds/quality/QuadQuality.h>
#include <iostream>
/*---------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::quality;
/*---------------------------------------------------------------------------*/
QuadQuality
QuadQuality::build(const math::Point& AP0, const math::Point& AP1,
								 const math::Point& AP2, const math::Point& AP3){
	return QuadQuality({math::vec(AP0),
	                 math::vec(AP1),
	                 math::vec(AP2),
	                 math::vec(AP3)});
}
/*---------------------------------------------------------------------------*/
double
QuadQuality::signedArea() const {
	return 0.25*a0()+0.25*a1()+0.25*a2()+0.25*a3();
}
/*---------------------------------------------------------------------------*/
double
QuadQuality::angleDeviation() const
{
		math::Vector3d e1 = L01().normalize();
	   math::Vector3d e2 = L12().normalize();
	   math::Vector3d e3 = L23().normalize();
	   math::Vector3d e4 = L30().normalize();

		auto teta1 = acos(e1.dot(e2)) * 180.0 / M_PI;
		auto teta2 = acos(e2.dot(e3)) * 180.0 / M_PI;
		auto teta3 = acos(e3.dot(e4)) * 180.0 / M_PI;
		auto teta4 = acos(e4.dot(e1)) * 180.0 / M_PI;

		auto teta_min = std::min( std::min(teta1, teta2), std::min(teta3, teta4));
		auto teta_max = std::max( std::max(teta1, teta2), std::max(teta3, teta4));

		return std::max( fabs(90.0-teta_min), fabs(teta_max-90.0) ) ;
}
/*---------------------------------------------------------------------------*/
double
QuadQuality::aspectRatio() const {
	return (lmax()*(l0()+l1()+l2()+l3()))/4*signedArea();
}
/*---------------------------------------------------------------------------*/
double
QuadQuality::condition() const {
	double x0 = (L0().norm2()+L3().norm2())/a0();
	double x1 = (L1().norm2()+L0().norm2())/a1();
	double x2 = (L2().norm2()+L1().norm2())/a2();
	double x3 = (L3().norm2()+L2().norm2())/a3();
	return 0.5*std::max(x0,std::max(x1,std::max(x2,x3)));
}
/*---------------------------------------------------------------------------*/
double
QuadQuality::jacobian() const {
	return std::min(a0(),std::min(a1(),std::min(a2(),a3())));
}
/*---------------------------------------------------------------------------*/
double
QuadQuality::scaledJacobian() const {

	double x0 =a0()/(l0()*l3());
	double x1 =a1()/(l1()*l0());
	double x2 =a2()/(l2()*l1());
	double x3 =a3()/(l3()*l2());

	return std::min(x0,std::min(x1,std::min(x2,x3)));

}
/*---------------------------------------------------------------------------*/
double
QuadQuality::minAngle() const {
	throw GMDSException("Not yet implemented!");
}
/*---------------------------------------------------------------------------*/
double
QuadQuality::maxAngle() const {
	throw GMDSException("Not yet implemented!");
}
/*---------------------------------------------------------------------------*/
