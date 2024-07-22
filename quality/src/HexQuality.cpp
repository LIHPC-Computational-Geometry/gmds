/*---------------------------------------------------------------------------*/
#include <gmds/quality/HexQuality.h>
#include <algorithm>
#include <limits>
/*---------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::quality;
/*---------------------------------------------------------------------------*/
HexQuality
HexQuality::build(const math::Point &AP0, const math::Point &AP1,
						const math::Point &AP2, const math::Point &AP3,
						const math::Point &AP4, const math::Point &AP5,
						const math::Point &AP6, const math::Point &AP7) {
	return HexQuality({math::vec(AP0),
							 math::vec(AP1),
							 math::vec(AP2),
							 math::vec(AP3),
							 math::vec(AP4),
							 math::vec(AP5),
							 math::vec(AP6),
							 math::vec(AP7)});
}
/*---------------------------------------------------------------------------*/
double HexQuality::volume() const{
	return A8().det()/64.0;
}
/*---------------------------------------------------------------------------*/
double HexQuality::diagonal() const {
	return dmin()/dmax();
}

/*---------------------------------------------------------------------------*/
double HexQuality::edgeRatio() const {
	return lmax()/lmin();
}
/*---------------------------------------------------------------------------*/
double HexQuality::maximumEdgeRatio() const
{
	double x1 = X1().norm();
	double x2 = X2().norm();
	double x3 = X3().norm();
	double a12 = std::max(x1/x2,x2/x1);
	double a13 = std::max(x1/x3,x3/x1);
	double a23 = std::max(x2/x3,x3/x2);
	return std::max(a12,std::max(a13,a23));
}
/*---------------------------------------------------------------------------*/
double HexQuality::jacobian() const
{
	double j[9]={A0().det(),A1().det(),A2().det(),
	              A3().det(),A4().det(),A5().det(),
	              A6().det(),A7().det(),A8().det()/64.0};
	double* ptr_min = std::min_element(j,j+9);
	return *ptr_min;

}
/*---------------------------------------------------------------------------*/
double HexQuality::scaledJacobian() const {
	double sj[9]={A0unit().det(),A1unit().det(),A2unit().det(),
					 A3unit().det(),A4unit().det(),A5unit().det(),
					 A6unit().det(),A7unit().det(),A8unit().det()};
	double* ptr_min = std::min_element(sj,sj+9);
	return *ptr_min;

}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
double HexQuality::lmin() const {
	double l[]={l0(),l1(),l2(),l3(),l4(),l5(),l6(),l7(),l8(),l9(),l10(),l11()};
	return *std::min_element(l,l+11);
}
/*---------------------------------------------------------------------------*/
double HexQuality::lmax() const {
	double l[]={l0(),l1(),l2(),l3(),l4(),l5(),l6(),l7(),l8(),l9(),l10(),l11()};
	return *std::max_element(l,l+11);
}
/*---------------------------------------------------------------------------*/
double HexQuality::maximumAspectFrobenius() const {

	math::Matrix33 A[8]= {A0(),A1(),A2(),A3(),A4(),A5(),A6(),A7()};
	double frob[8];
	for(auto i=0; i<8; ++i)
	{
		double det = A[i].det();
		if(det <= std::numeric_limits<double>::min())
			return std::numeric_limits<double>::max();


		double t1 = A[i][0].dot(A[i][0]) + A[i][1].dot(A[i][1]) + A[i][2].dot(A[i][2]);
		float t2 = (A[i][0].cross(A[i][1])).dot(A[i][0].cross(A[i][1])) +
		           (A[i][1].cross(A[i][2])).dot(A[i][1].cross(A[i][2])) +
		           (A[i][2].cross(A[i][0])).dot(A[i][2].cross(A[i][0]));
		frob[i]  = sqrt(t1*t2)/(det*3);
	}
	return *std::max_element(frob, frob+8);
}
/*---------------------------------------------------------------------------*/
double HexQuality::meanAspectFrobenius() const {

	math::Matrix33 A[8]= {A0(),A1(),A2(),A3(),A4(),A5(),A6(),A7()};
	double frob=0;
	for(auto i=0; i<8; ++i)
	{
		double det = A[i].det();
		if(det <= std::numeric_limits<double>::min())
			return std::numeric_limits<double>::max();


		double t1 = A[i][0].dot(A[i][0]) + A[i][1].dot(A[i][1]) + A[i][2].dot(A[i][2]);
		float t2 = (A[i][0].cross(A[i][1])).dot(A[i][0].cross(A[i][1])) +
		           (A[i][1].cross(A[i][2])).dot(A[i][1].cross(A[i][2])) +
		           (A[i][2].cross(A[i][0])).dot(A[i][2].cross(A[i][0]));
		frob+= sqrt(t1*t2)/(det*3);
	}
	return frob;
}
/*---------------------------------------------------------------------------*/
double HexQuality::oddy() const {
	static double four_over_three = 4.0/3.0;
	math::Matrix33 A[9]= {A0(),A1(),A2(),A3(),A4(),
	                       A5(),A6(),A7(),A8()};

	double oddy[9];
	for(auto i=0; i<9; ++i)
	{
		double det = A[i].det();
		if(det > std::numeric_limits<double>::min())
		{
			double a11 = A[i][0].dot(A[i][0]);
			double a12 = A[i][0].dot(A[i][1]);
			double a13 = A[i][0].dot(A[i][2]);
			double a22 = A[i][1].dot(A[i][1]);
			double a23 = A[i][1].dot(A[i][2]);
			double a33 = A[i][2].dot(A[i][2]);

			double AtA_sqrd = a11*a11 + 2.0*a12*a12 + 2.0*a13*a13 + a22*a22 + 2.0*a23*a23 +a33*a33;
			double A_sqrd   = a11 + a22 + a33;

			oddy[i] = (AtA_sqrd - A_sqrd*A_sqrd/3.0) / pow(det,four_over_three);
		}
		else return std::numeric_limits<double>::max();
	}
	return *std::max_element(oddy,oddy+9);
}
/*---------------------------------------------------------------------------*/
double HexQuality::shape() const {
	static double two_over_three = 2.0/3.0;
	math::Matrix33 A[9]= {A0(),A1(),A2(),A3(),A4(),
	                       A5(),A6(),A7(),A8()};
	double sh[9];
	for(auto i=0; i<9; ++i){
		float det = A[i].det();
		if(det<=std::numeric_limits<double>::min())
			return 0;
		double n = pow(det, two_over_three);
		double d = A[i][0].dot(A[i][0]) + A[i][1].dot(A[i][1]) + A[i][2].dot(A[i][2]);
		if(d<=std::numeric_limits<double>::min())
			return 0;
		sh[i] = 3.0 * n/d;
	}
	return *std::min_element(sh, sh+9);
}

/*---------------------------------------------------------------------------*/
double HexQuality::shear() const {
	math::Matrix33 A[9]= {A0(),A1(),A2(),A3(),A4(),
	                       A5(),A6(),A7(),A8()};

	double sh[9];
	for(auto i=0; i<9; ++i)
	{
		sh[i] = A[i].det();
		if(sh[i]<=std::numeric_limits<double>::min())
		              return 0;
	}
	return *std::min_element(sh, sh+9);
}
/*---------------------------------------------------------------------------*/
double HexQuality::skew() const {

	if(X1().norm() <= std::numeric_limits<double>::min())
		return 0;
	if(X2().norm() <= std::numeric_limits<double>::min())
		return 0;
	if(X3().norm() <= std::numeric_limits<double>::min())
		return 0;

	double sk[3] ={
	      std::fabs(X1().normalize().dot(X2().normalize())),
	      std::fabs(X1().normalize().dot(X3().normalize())),
	      std::fabs(X2().normalize().dot(X3().normalize()))};

	return *std::max_element(sk,sk+3);
}

/*---------------------------------------------------------------------------*/
double HexQuality::stretch() const {

	double  l[12]={l0(),l1(),l2(),l3(),l4(),l5(),l6(),l7(),l8(),l9(),l10(),l1()};
	double  d[4]={D0().norm(),D1().norm(),D2().norm(),D3().norm()};

	return 1.732050807568877f * (*std::min_element(l, l+12)) /
			 (*std::max_element(d, d+4));
}

/*---------------------------------------------------------------------------*/
double HexQuality::taper() const{
	double     x[3]= {X1().norm(),X2().norm(),X3().norm()};
	double    cd[3]={X12().norm(),X13().norm(),X23().norm()};

	if(x[0] <= std::numeric_limits<double>::min())
		return std::numeric_limits<double>::max();
	if(x[1] <= std::numeric_limits<double>::min())
		return std::numeric_limits<double>::max();
	if(x[2] <= std::numeric_limits<double>::min())
		return std::numeric_limits<double>::max();

	double taper[3] ={
	      cd[0] / std::min(x[0], x[1]),
	      cd[1] / std::min(x[0], x[2]),
	      cd[2] / std::min(x[1], x[2]),
	   };

	return *std::max_element(taper, taper+3);
}
