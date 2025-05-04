
/*----------------------------------------------------------------------------*/
/*
 * Pyramid.cpp
 *
 *  Created on: 27 nov 2014
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#include <gmds/math/Pyramid.h>
/*----------------------------------------------------------------------------*/
#include <gmds/math/Tetrahedron.h>
#include <gmds/math/Triangle.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
Pyramid::Pyramid()
{
	m_pnts[0] = Point(0,0,0);
	m_pnts[1] = Point(1,0,0);
	m_pnts[2] = Point(1,1,0);
	m_pnts[3] = Point(0,1,0);
	m_pnts[4] = Point(0,0,1);
}
/*----------------------------------------------------------------------------*/
Pyramid::Pyramid(const Point& AP0, const Point& AP1, const Point& AP2, const Point& AP3, const Point& AP4)
{
	m_pnts[0] = AP0;
	m_pnts[1] = AP1;
	m_pnts[2] = AP2;
	m_pnts[3] = AP3;
	m_pnts[4] = AP4;
}
/*----------------------------------------------------------------------------*/
Pyramid::Pyramid(Point APoints[5])
{
        m_pnts[0] = APoints[0];
        m_pnts[1] = APoints[1];
        m_pnts[2] = APoints[2];
        m_pnts[3] = APoints[3];
        m_pnts[4] = APoints[4];
}
/*----------------------------------------------------------------------------*/
Pyramid::Pyramid(const std::vector<Point>& APoints)
{
	if(APoints.size() != 5) {
		throw GMDSException("Pyramid::Pyramid APoints not of size 5.");
	}

        m_pnts[0] = APoints[0];
        m_pnts[1] = APoints[1];
        m_pnts[2] = APoints[2];
        m_pnts[3] = APoints[3];
        m_pnts[4] = APoints[4];
}
/*----------------------------------------------------------------------------*/
Pyramid::Pyramid(const Pyramid& APyr)
{
	m_pnts[0] = APyr.m_pnts[0];
	m_pnts[1] = APyr.m_pnts[1];
	m_pnts[2] = APyr.m_pnts[2];
	m_pnts[3] = APyr.m_pnts[3];
	m_pnts[4] = APyr.m_pnts[4];
}
/*----------------------------------------------------------------------------*/
Pyramid::~Pyramid(){}
/*----------------------------------------------------------------------------*/
const Point& Pyramid::getPoint(const TInt& AIndex) const
{
	return m_pnts[AIndex];
}
/*----------------------------------------------------------------------------*/
Point Pyramid::getCenter() const
{
	TCoord coordX = 0.;
	TCoord coordY = 0.;
	TCoord coordZ = 0.;
	
	for(const auto & m_pnt : m_pnts) {
		coordX += m_pnt.X();
		coordY += m_pnt.Y();
		coordZ += m_pnt.Z();
	}
	coordX /= 5.;
	coordY /= 5.;
	coordZ /= 5.;	

        return Point(coordX,coordY,coordZ);
}
/*----------------------------------------------------------------------------*/
 double Pyramid::getVolume() const
{       
        Point baseCenter((1./4.) * (m_pnts[0]+m_pnts[1]+m_pnts[2]+m_pnts[3]));
	Tetrahedron t0(baseCenter,m_pnts[0],m_pnts[1],m_pnts[4]);
	Tetrahedron t1(baseCenter,m_pnts[1],m_pnts[2],m_pnts[4]);
	Tetrahedron t2(baseCenter,m_pnts[2],m_pnts[3],m_pnts[4]);
	Tetrahedron t3(baseCenter,m_pnts[3],m_pnts[0],m_pnts[4]);

	double vol = t0.getVolume() + t1.getVolume() + t2.getVolume() + t3.getVolume();
	return vol;
}
/*----------------------------------------------------------------------------*/
double
Pyramid::computeScaledJacobian() const	
{
	double scaledJacobianMin = HUGE_VAL;

	// compute the jacobian for the 4 tetrahedra formed by the top vertex 
        // and the four base vertices
	Tetrahedron t0(m_pnts[0],m_pnts[1],m_pnts[2],m_pnts[4]);
	Tetrahedron t1(m_pnts[1],m_pnts[2],m_pnts[3],m_pnts[4]);
	Tetrahedron t2(m_pnts[2],m_pnts[3],m_pnts[0],m_pnts[4]);
	Tetrahedron t3(m_pnts[3],m_pnts[0],m_pnts[1],m_pnts[4]);

	if(t0.computeScaledJacobian() < scaledJacobianMin) {
		scaledJacobianMin = t0.computeScaledJacobian();
	}
	if(t1.computeScaledJacobian() < scaledJacobianMin) {
                scaledJacobianMin = t1.computeScaledJacobian();
        }
	if(t2.computeScaledJacobian() < scaledJacobianMin) {
                scaledJacobianMin = t2.computeScaledJacobian();
        }
	if(t3.computeScaledJacobian() < scaledJacobianMin) {
                scaledJacobianMin = t3.computeScaledJacobian();
        }

	return scaledJacobianMin;
}
/*----------------------------------------------------------------------------*/
double
Pyramid::computeNormalizedScaledJacobian() const	
{
        double scaledJ = this->computeScaledJacobian();
	scaledJ *= 2.;
	
	if(scaledJ >  1.) scaledJ =  1.;
	if(scaledJ < -1.) scaledJ = -1.;
	
	return scaledJ;
}
/*----------------------------------------------------------------------------*/
double
Pyramid::computeMeanRatio()
{
	throw GMDSException("Pyramid::computeMeanRatio not available yet.");

      /*  // We use only the four vertices of the base of the pyramid
        const int neighbors[4][3] =
        {
                {1,3,4},
                {2,0,4},
                {3,1,4},
                {0,2,4},
        };
 
	TCoord matValues[3][3];
	matValues[0][0] = 1.;
        matValues[1][0] = 0.;
        matValues[2][0] = 0.;
        matValues[0][1] = 0.;
        matValues[1][1] = 1.;
        matValues[2][1] = 0.;
        matValues[0][2] = 1./2.;
        matValues[1][2] = 1./2.;
        matValues[2][2] = std::sqrt(2.)/2.;
	const math::Matrix<3,3,double> inverseW((math::Matrix<3,3,double> (matValues)).inverse());

	double meanRatio = 0.;

        for (int iVertex = 0; iVertex < 4; iVertex++) {                

                int i0    = neighbors[iVertex][0];
                int i1    = neighbors[iVertex][1];
                int i2    = neighbors[iVertex][2];
		
		Vector3d v0(m_pnts[i0],m_pnts[iVertex]);
		Vector3d v1(m_pnts[i1],m_pnts[iVertex]);
		Vector3d v2(m_pnts[i2],m_pnts[iVertex]);

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
			throw GMDSException("Pyramid::computeMeanRatio frob is zero.");
		}

        	meanRatio += (3. * std::cbrt(std::pow(det,2))) / frob;
        }

	return (1./4.) * meanRatio;	*/
}
/*----------------------------------------------------------------------------*/
double
Pyramid::computeMeanEdgeLength() const
{
  double sumLength = 0;
  sumLength += m_pnts[0].distance(m_pnts[1]);
  sumLength += m_pnts[1].distance(m_pnts[2]);
  sumLength += m_pnts[2].distance(m_pnts[3]);
  sumLength += m_pnts[3].distance(m_pnts[0]);
  sumLength += m_pnts[0].distance(m_pnts[4]);
  sumLength += m_pnts[1].distance(m_pnts[4]);
  sumLength += m_pnts[2].distance(m_pnts[4]);
  sumLength += m_pnts[3].distance(m_pnts[4]);

  sumLength /= 8.;  
  return sumLength;
}
/*----------------------------------------------------------------------------*/
bool
Pyramid::intersect(const Triangle& ATri, const bool AProper) const
{
	Triangle T1(m_pnts[0], m_pnts[1], m_pnts[2]);
	bool inter = ATri.intersect(T1, AProper);
	if (inter) {
		return true;
	}

	Triangle T2(m_pnts[0], m_pnts[2], m_pnts[3]);
	inter = ATri.intersect(T2, AProper);
	if (inter) {
		return true;
	}

	Triangle T3(m_pnts[0], m_pnts[1], m_pnts[4]);
	inter = ATri.intersect(T3, AProper);
	if (inter) {
		return true;
	}

	Triangle T4(m_pnts[1], m_pnts[2], m_pnts[4]);
	inter = ATri.intersect(T4, AProper);
	if (inter) {
		return true;
	}

        Triangle T5(m_pnts[2], m_pnts[3], m_pnts[4]);
        inter = ATri.intersect(T5, AProper);
        if (inter) {
                return true;
        }

        Triangle T6(m_pnts[3], m_pnts[0], m_pnts[4]);
        inter = ATri.intersect(T6, AProper);
        if (inter) {
                return true;
        }

	return false;	
}
/*---------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& AStr, const Pyramid& APyr){
	AStr<<"Pyramid"<<std::endl;
	for(int iPoint=0; iPoint<5; iPoint++) {
		AStr<<std::endl<<APyr.getPoint(iPoint);
	}
        return AStr;
}
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/




