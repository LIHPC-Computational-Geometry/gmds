
/*----------------------------------------------------------------------------*/
/*
 * Tetrahedron.cpp
 *
 *  Created on: 27 nov 2014
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#include <gmds/math/Tetrahedron.h>
/*----------------------------------------------------------------------------*/
#include <gmds/math/Numerics.h>
#include <gmds/math/Vector.h>
#include <gmds/math/Triangle.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
Tetrahedron::Tetrahedron()
{
	m_pnts[0] = Point(0,0,0);
	m_pnts[1] = Point(1,0,0);
	m_pnts[2] = Point(1,1,0);
	m_pnts[3] = Point(1,0,1);
}
/*----------------------------------------------------------------------------*/
Tetrahedron::Tetrahedron(const Point& AP0, const Point& AP1, const Point& AP2, const Point& AP3)
{
	m_pnts[0] = AP0;
	m_pnts[1] = AP1;
	m_pnts[2] = AP2;
	m_pnts[3] = AP3;
}
/*----------------------------------------------------------------------------*/
Tetrahedron::Tetrahedron(Point APoints[4])
{
        m_pnts[0] = APoints[0];
        m_pnts[1] = APoints[1];
        m_pnts[2] = APoints[2];
        m_pnts[3] = APoints[3];
}
/*----------------------------------------------------------------------------*/
Tetrahedron::Tetrahedron(const std::vector<Point>& APoints)
{
	if(APoints.size() != 4) {
                throw GMDSException("Tetrahedron::Tetrahedron APoints not of size 4.");
        }

        m_pnts[0] = APoints[0];
        m_pnts[1] = APoints[1];
        m_pnts[2] = APoints[2];
        m_pnts[3] = APoints[3];
}
/*----------------------------------------------------------------------------*/
Tetrahedron::Tetrahedron(const Tetrahedron& ATet)
{
	m_pnts[0] = ATet.m_pnts[0];
	m_pnts[1] = ATet.m_pnts[1];
	m_pnts[2] = ATet.m_pnts[2];
	m_pnts[3] = ATet.m_pnts[3];
}
/*----------------------------------------------------------------------------*/
Tetrahedron::~Tetrahedron(){}
/*----------------------------------------------------------------------------*/
const Point& Tetrahedron::getPoint(const TInt& AIndex) const
{
	return m_pnts[AIndex];
}
    /*----------------------------------------------------------------------------*/
    const Point Tetrahedron::getCenter() const
    {
        TCoord coordX = 0.;
        TCoord coordY = 0.;
        TCoord coordZ = 0.;
        
        for(int iPoint=0; iPoint<4; iPoint++) {
            coordX += m_pnts[iPoint].X();
            coordY += m_pnts[iPoint].Y();
            coordZ += m_pnts[iPoint].Z();
        }
        coordX /= 4.;
        coordY /= 4.;
        coordZ /= 4.;	
        
        return Point(coordX,coordY,coordZ);
    }
    /*----------------------------------------------------------------------------*/
     TCoord Tetrahedron::getVolume() const
    {
        math::Vector3d v01=m_pnts[1]-m_pnts[0];
        math::Vector3d v02=m_pnts[2]-m_pnts[0];
        math::Vector3d v03=m_pnts[3]-m_pnts[0];
        return (v03.dot(v01.cross(v02)))/ 6.0;
    }
/*----------------------------------------------------------------------------*/
//double
//Tetrahedron::computeMinJacobian() const
//{
//        Matrix<3,3,double> A = this->jacobian();
//
//        return A.det();
//}
/*----------------------------------------------------------------------------*/
double
Tetrahedron::computeScaledJacobian() const
{
        int neighbors[4][3] =
        {
                {1,2,3},
                {2,0,3},
                {0,1,3},
                {0,1,2},
        };

	Matrix<3,3,double> A = this->jacobian();

	double maxLength = -HUGE_VALF;
	for (int iVertex = 0; iVertex < 4; iVertex++) {

		int i0    = neighbors[iVertex][0];
		int i1    = neighbors[iVertex][1];
		int i2    = neighbors[iVertex][2];
		double l0 = (m_pnts[iVertex]-m_pnts[i0]).norm();
		double l1 = (m_pnts[iVertex]-m_pnts[i1]).norm();
		double l2 = (m_pnts[iVertex]-m_pnts[i2]).norm();
		double length = l0*l1*l2;

		if(length > maxLength) {
			maxLength = length;
		}
	}

	if(maxLength == 0.) {
		return 0.;
	} else {

		double scaledJmin = A.det() / maxLength;

		return scaledJmin;
	}
}
/*----------------------------------------------------------------------------*/
		double
		Tetrahedron::computeNormalizedScaledJacobian() const
		{
        double scaledJ = sqrt(2.)*this->computeScaledJacobian();
	if(scaledJ >  1.) scaledJ =  1.;
	if(scaledJ < -1.) scaledJ = -1.;

	return scaledJ;
}
/*----------------------------------------------------------------------------*/
double
Tetrahedron::computeMeanRatio() const
{
	throw GMDSException("Tetrahedron::computeMeanRatio not available yet.");
/*
        const int neighbors[4][3] =
        {
                {1,2,3},
                {2,0,3},
                {0,1,3},
                {0,1,2},
        };
 
	TCoord matValues[3][3];
	matValues[0][0] = 1.;
        matValues[1][0] = 0.;
        matValues[2][0] = 0.;
        matValues[0][1] = 1./2.;
        matValues[1][1] = std::sqrt(3.)/2.;
        matValues[2][1] = 0.;
        matValues[0][2] = 1./2.;
        matValues[1][2] = std::sqrt(3.)/6.;
        matValues[2][2] = std::sqrt(2./3.);
	const math::Matrix<3,3,double> inverseW((math::Matrix<3,3,double> (matValues)).inverse());

	double meanRatio = 0.;

        for (int iVertex = 0; iVertex < 4; iVertex++) {                

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
			throw GMDSException("Tetrahedron::computeMeanRatio frob is zero.");
		}

        	meanRatio += (3. * std::cbrt(std::pow(det,2))) / frob;
        }

	return (1./4.) * meanRatio;	*/
}
/*----------------------------------------------------------------------------*/
double
Tetrahedron::computeMeanEdgeLength() const
{
  double sumLength = 0;
  sumLength += m_pnts[0].distance(m_pnts[1]);
  sumLength += m_pnts[1].distance(m_pnts[2]);
  sumLength += m_pnts[2].distance(m_pnts[0]);
  sumLength += m_pnts[0].distance(m_pnts[3]);
  sumLength += m_pnts[1].distance(m_pnts[3]);
  sumLength += m_pnts[2].distance(m_pnts[3]);

  sumLength /= 6.;  
  return sumLength;
}
/*----------------------------------------------------------------------------*/
void
Tetrahedron::computeBoundingBox(TCoord AMinXYZ[3], TCoord AMaxXYZ[3]) const
{
	AMinXYZ[0] = min4(m_pnts[0].X(),m_pnts[1].X(),m_pnts[2].X(),m_pnts[3].X());
	AMinXYZ[1] = min4(m_pnts[0].Y(),m_pnts[1].Y(),m_pnts[2].Y(),m_pnts[3].Y());
	AMinXYZ[2] = min4(m_pnts[0].Z(),m_pnts[1].Z(),m_pnts[2].Z(),m_pnts[3].Z());
	AMaxXYZ[0] = max4(m_pnts[0].X(),m_pnts[1].X(),m_pnts[2].X(),m_pnts[3].X());
	AMaxXYZ[1] = max4(m_pnts[0].Y(),m_pnts[1].Y(),m_pnts[2].Y(),m_pnts[3].Y());
	AMaxXYZ[2] = max4(m_pnts[0].Z(),m_pnts[1].Z(),m_pnts[2].Z(),m_pnts[3].Z());
}
/*----------------------------------------------------------------------------*/
math::Matrix<3,3,double>
Tetrahedron::jacobian() const
{
	return{m_pnts[1]-m_pnts[0],m_pnts[2]-m_pnts[0],m_pnts[3]-m_pnts[0]};
}
/*----------------------------------------------------------------------------*/
bool
Tetrahedron::intersect(const Triangle& ATri, const bool AProper) const
{
	Triangle T1(m_pnts[0], m_pnts[1], m_pnts[2]);
	bool inter = ATri.intersect(T1, AProper);
	if (inter) {
		return true;
	}

	Triangle T2(m_pnts[0], m_pnts[1], m_pnts[3]);
	inter = ATri.intersect(T2, AProper);
	if (inter) {
		return true;
	}

	Triangle T3(m_pnts[1], m_pnts[2], m_pnts[3]);
	inter = ATri.intersect(T3, AProper);
	if (inter) {
		return true;
	}

	Triangle T4(m_pnts[2], m_pnts[0], m_pnts[3]);
	inter = ATri.intersect(T4, AProper);
	if (inter) {
		return true;
	}

	return false;	
}
/*---------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& AStr, const Tetrahedron& ATet){
	AStr<<"Tetahedron"<<std::endl;
	for(int iPoint=0; iPoint<4; iPoint++) {
		AStr<<std::endl<<ATet.getPoint(iPoint);
	}
        return AStr;
}
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/




