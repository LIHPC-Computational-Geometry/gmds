
/*----------------------------------------------------------------------------*/
/*
 * Prims3.cpp
 *
 *  Created on: 27 nov 2014
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#include "gmds/math/Prism3.h"
/*----------------------------------------------------------------------------*/
#include "gmds/math/Vector.h"
#include "gmds/math/Triangle.h"
#include "gmds/math/Tetrahedron.h"
#include "gmds/math/Pyramid.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
Prism3::Prism3()
{
	m_pnts[0] = Point(0,0,0);
	m_pnts[1] = Point(1,0,0);
	m_pnts[2] = Point(0,1,0);
	m_pnts[3] = Point(0,0,1);
	m_pnts[4] = Point(1,0,1);
	m_pnts[5] = Point(0,1,1);
}
/*----------------------------------------------------------------------------*/
Prism3::Prism3(const Point& AP0, const Point& AP1, const Point& AP2,
		const Point& AP3, const Point& AP4, const Point& AP5)
{
	m_pnts[0] = AP0;
	m_pnts[1] = AP1;
	m_pnts[2] = AP2;
	m_pnts[3] = AP3;
	m_pnts[4] = AP4;
	m_pnts[5] = AP5;
}
/*----------------------------------------------------------------------------*/
Prism3::Prism3(Point APoints[6])
{
        m_pnts[0] = APoints[0];
        m_pnts[1] = APoints[1];
        m_pnts[2] = APoints[2];
        m_pnts[3] = APoints[3];
        m_pnts[4] = APoints[4];
        m_pnts[5] = APoints[5];
}
/*----------------------------------------------------------------------------*/
Prism3::Prism3(const std::vector<Point>& APoints)
{
        if(APoints.size() != 6) {
                throw GMDSException("Prism3::Prism3 APoints not of size 6.");
        }

        m_pnts[0] = APoints[0];
        m_pnts[1] = APoints[1];
        m_pnts[2] = APoints[2];
        m_pnts[3] = APoints[3];
        m_pnts[4] = APoints[4];
	m_pnts[5] = APoints[5];
}
/*----------------------------------------------------------------------------*/
Prism3::Prism3(const Prism3& APrsm)
{
	m_pnts[0] = APrsm.m_pnts[0];
	m_pnts[1] = APrsm.m_pnts[1];
	m_pnts[2] = APrsm.m_pnts[2];
	m_pnts[3] = APrsm.m_pnts[3];
	m_pnts[4] = APrsm.m_pnts[4];
	m_pnts[5] = APrsm.m_pnts[5];
}
/*----------------------------------------------------------------------------*/
Prism3::~Prism3(){}
/*----------------------------------------------------------------------------*/
const Point& Prism3::getPoint(const TInt& AIndex) const
{
	return m_pnts[AIndex];
}
/*----------------------------------------------------------------------------*/
Point Prism3::getCenter() const
{
	TCoord coordX = 0.;
	TCoord coordY = 0.;
	TCoord coordZ = 0.;
	
	for(const auto & m_pnt : m_pnts) {
		coordX += m_pnt.X();
		coordY += m_pnt.Y();
		coordZ += m_pnt.Z();
	}
	coordX /= 6.;
    coordY /= 6.;
    coordZ /= 6.;

    return Point(coordX,coordY,coordZ);
}
/*----------------------------------------------------------------------------*/
    double
    Prism3::getVolume() const
    {

        Point center(getCenter());
        Tetrahedron t0(m_pnts[0],m_pnts[1],m_pnts[2],center);
        Tetrahedron t1(m_pnts[5],m_pnts[4],m_pnts[3],center);
        Pyramid p0(m_pnts[0],m_pnts[3],m_pnts[4],m_pnts[1],center);
        Pyramid p1(m_pnts[1],m_pnts[4],m_pnts[5],m_pnts[2],center);
        Pyramid p2(m_pnts[0],m_pnts[2],m_pnts[5],m_pnts[3],center);

        double vol = t0.getVolume() + t1.getVolume() + p0.getVolume() + p1.getVolume() + p2.getVolume();
        return vol;
    }
/*----------------------------------------------------------------------------*/
    double
    Prism3::computeScaledJacobian() const
    {
        int neighbors[6][3] =
                {
                        {1,2,3},
                        {2,0,4},
                        {0,1,5},
                        {0,5,4},
                        {1,3,5},
                        {2,4,3}
                };

        double scaledJ[6];
        for (int iVertex = 0; iVertex < 6; iVertex++) {

            Matrix<3,3,double> A = this->jacobian(iVertex);

            int i0    = neighbors[iVertex][0];
            int i1    = neighbors[iVertex][1];
            int i2    = neighbors[iVertex][2];
            double l0 = (m_pnts[iVertex]-m_pnts[i0]).norm();
            double l1 = (m_pnts[iVertex]-m_pnts[i1]).norm();
            double l2 = (m_pnts[iVertex]-m_pnts[i2]).norm();
            double l012 = l0*l1*l2;
            double det = A.det();
            if (l012 == 0 || det == 0)
                scaledJ[iVertex] = 0.;
            else
                scaledJ[iVertex] = det/l012;
        }
        double scaledJmin = HUGE_VALF;
        for(double iVertex : scaledJ) {
            if(iVertex < scaledJmin) {
                scaledJmin = iVertex;
            }
        }

        return scaledJmin;
    }
/*----------------------------------------------------------------------------*/
    double
Prism3::computeNormalizedScaledJacobian() const
{
        double scaledJ = this->computeScaledJacobian(); 
	scaledJ *= 2./sqrt(3.);
	
	if(scaledJ >  1.) scaledJ =  1.;
	if(scaledJ < -1.) scaledJ = -1.;
	
	return scaledJ;
}
/*----------------------------------------------------------------------------*/
double
Prism3::computeMeanRatio()
{
	throw GMDSException("Prism3::computeMeanRatio not available yet.");

    /*    const int neighbors[6][3] =
        {
                {1,2,3},
                {2,0,4},
                {0,1,5},
                {5,4,0},
		{3,5,1},
		{3,4,2}
        };
 
	TCoord matValues[3][3];
	matValues[0][0] = 1.;
        matValues[1][0] = 0.;
        matValues[2][0] = 0.;
        matValues[0][1] = 1./2.;
        matValues[1][1] = std::sqrt(3.)/2.;
        matValues[2][1] = 0.;
        matValues[0][2] = 0.;
        matValues[1][2] = 0.;
        matValues[2][2] = 1.;
	const math::Matrix<3,3,double> inverseW((math::Matrix<3,3,double> (matValues)).inverse());

	double meanRatio = 0.;

        for (int iVertex = 0; iVertex < 6; iVertex++) {                

                int i0    = neighbors[iVertex][0];
                int i1    = neighbors[iVertex][1];
                int i2    = neighbors[iVertex][2];
		
		math::Vector v0(m_pnts[i0] - m_pnts[iVertex]);
		math::Vector v1(m_pnts[i1] - m_pnts[iVertex]);
		math::Vector v2(m_pnts[i2] - m_pnts[iVertex]);

		TCoord matValues_dtk[3][3];
		matValues_dtk[0][0] = v0.X();
		matValues_dtk[1][0] = v0.Y();
		matValues_dtk[2][0] = v0.Z();
		matValues_dtk[0][1] = v1.X();
		matValues_dtk[1][1] = v1.Y();
		matValues_dtk[2][1] = v1.Z();
		matValues_dtk[0][2] = v2.X();
		matValues_dtk[1][2] = v2.Y();
		matValues_dtk[2][2] = v2.Z();
		const math::Matrix<3,3,double> dtk (matValues_dtk);
		const math::Matrix<3,3,double> sk (dtk * inverseW);
		
		double det = sk.det();
		double frob = sk.frobeniusNorm2();
		
		if(frob == 0.) {
			throw GMDSException("Prism3::computeMeanRatio frob is zero.");
		}

        	meanRatio += (3. * std::cbrt(std::pow(det,2))) / frob;
        }

	return (1./6.) * meanRatio;*/
}
/*----------------------------------------------------------------------------*/
double
Prism3::computeMeanEdgeLength() const
{
  double sumLength = 0;
  sumLength += m_pnts[0].distance(m_pnts[1]);
  sumLength += m_pnts[1].distance(m_pnts[2]);
  sumLength += m_pnts[2].distance(m_pnts[0]);
  sumLength += m_pnts[3].distance(m_pnts[4]);
  sumLength += m_pnts[4].distance(m_pnts[5]);
  sumLength += m_pnts[5].distance(m_pnts[3]);
  sumLength += m_pnts[0].distance(m_pnts[3]);
  sumLength += m_pnts[1].distance(m_pnts[4]);
  sumLength += m_pnts[2].distance(m_pnts[5]);

  sumLength /= 9.;  
  return sumLength;
}
/*---------------------------------------------------------------------------*/
math::Matrix<3,3,double>
Prism3::jacobian(const int iVertex) const
{
	math::Matrix<3,3,double> mat;

	switch(iVertex) {
	case 0:
		mat= {	 m_pnts[1]- m_pnts[0],
					 m_pnts[2]- m_pnts[0],
					 m_pnts[3]- m_pnts[0]};
		break;
	case 1:

		mat= {	 m_pnts[1]- m_pnts[0],
					 m_pnts[2]- m_pnts[0],
					 m_pnts[4]- m_pnts[1]};
		break;
	case 2:

		mat= {	 m_pnts[1]- m_pnts[0],
					 m_pnts[2]- m_pnts[0],
					 m_pnts[5]- m_pnts[2]};
		break;
	case 3:

		mat= {	 m_pnts[4]- m_pnts[3],
					 m_pnts[5]- m_pnts[3],
					 m_pnts[3]- m_pnts[0]};
		break;
	case 4:
		mat= {	 m_pnts[4]- m_pnts[3],
					 m_pnts[5]- m_pnts[3],
					 m_pnts[4]- m_pnts[1]};
		break;
	case 5:
		mat= {	 m_pnts[4]- m_pnts[3],
					 m_pnts[5]- m_pnts[3],
					 m_pnts[5]- m_pnts[2]};
		break;
	default:
		throw GMDSException("Prism3::jacobian invalid vertex number.");
		break;
	}

	return mat;
}
/*----------------------------------------------------------------------------*/
bool
Prism3::intersect(const Triangle& ATri, const bool AProper) const
{
	Triangle T1(m_pnts[0], m_pnts[1], m_pnts[2]);
	bool inter = ATri.intersect(T1, AProper);
	if (inter) {
		return true;
	}

	Triangle T2(m_pnts[3], m_pnts[4], m_pnts[5]);
	inter = ATri.intersect(T2, AProper);
	if (inter) {
		return true;
	}

	Triangle T3(m_pnts[0], m_pnts[1], m_pnts[4]);
	inter = ATri.intersect(T3, AProper);
	if (inter) {
		return true;
	}

	Triangle T4(m_pnts[0], m_pnts[4], m_pnts[3]);
	inter = ATri.intersect(T4, AProper);
	if (inter) {
		return true;
	}

        Triangle T5(m_pnts[1], m_pnts[2], m_pnts[5]);
        inter = ATri.intersect(T5, AProper);
        if (inter) {
                return true;
        }

        Triangle T6(m_pnts[1], m_pnts[5], m_pnts[4]);
        inter = ATri.intersect(T6, AProper);
        if (inter) {
                return true;
        }

        Triangle T7(m_pnts[2], m_pnts[0], m_pnts[3]);
        inter = ATri.intersect(T7, AProper);
        if (inter) {
                return true;
        }

        Triangle T8(m_pnts[2], m_pnts[3], m_pnts[5]);
        inter = ATri.intersect(T8, AProper);
        if (inter) {
                return true;
        }

	return false;	
}
/*---------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& AStr, const Prism3& APrsm){
	AStr<<"Prism3"<<std::endl;
	for(int iPoint=0; iPoint<6; iPoint++) {
		AStr<<std::endl<<APrsm.getPoint(iPoint);
	}
        return AStr;
}
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/




