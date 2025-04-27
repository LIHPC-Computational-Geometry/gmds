
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
    const Point Tetrahedron::getCircumcenter() const
    {
	    // Get the tetra vertices
	    math::Point A = m_pnts[0];
	    math::Point B = m_pnts[1];
	    math::Point C = m_pnts[2];
	    math::Point D = m_pnts[3];
	    // Get the coordinates
	    double xA = A.X();
	    double yA = A.Y();
	    double zA = A.Z();
	    double xB = B.X();
	    double yB = B.Y();
	    double zB = B.Z();
	    double xC = C.X();
	    double yC = C.Y();
	    double zC = C.Z();
	    double xD = D.X();
	    double yD = D.Y();
	    double zD = D.Z();
	    math::Matrix33 M;
	    M(0,0) = xB-xA;
	    M(0,1) = yB-yA;
	    M(0,2) = zB-zA;
	    M(1,0) = xC-xB;
	    M(1,1) = yC-yB;
	    M(1,2) = zC-zB;
	    M(2,0) = xD-xA;
	    M(2,1) = yD-yA
	       ;
	    M(2,2) = zD-zA;
	    double xV = (xB-xA)*(1./2.)*(xA+xB)+(yB-yA)*(1./2.)*(yA+yB)+(zB-zA)*(1./2.)*(zA+zB);
	    double yV = (xC-xB)*(1./2.)*(xB+xC)+(yC-yB)*(1./2.)*(yB+yC)+(zC-zB)*(1./2.)*(zB+zC);
	    double zV = (xD-xA)*(1./2.)*(xA+xD)+(yD-yA)*(1./2.)*(yA+yD)+(zD-zA)*(1./2.)*(zA+zD);
	    math::Vector3d V;
	    V.set(0, xV);
	    V.set(1, yV);
	    V.set(2, zV);
	    math::Vector3d X = M.solve(V);
	    math::Point Ctr(X.X(), X.Y(), X.Z());
	    return Ctr;
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
Tetrahedron::computeQuality() const
{
	const math::Point p0 = m_pnts[0];
	const math::Point p1 = m_pnts[1];
	const math::Point p2 = m_pnts[2];
	const math::Point p3 = m_pnts[3];

	const math::Vector3d v01 = p1-p0;
	const math::Vector3d v21 = p2-p1;
	const math::Vector3d v02 = p0-p2;
	const math::Vector3d v30 = p3-p0;
	const math::Vector3d v31 = p3-p1;
	const math::Vector3d v32 = p3-p2;

	const double l01 = v01.norm();
	const double l21 = v21.norm();
	const double l02 = v02.norm();
	const double l30 = v30.norm();
	const double l31 = v31.norm();
	const double l32 = v32.norm();
  double volume = std::abs(getVolume());

	return (36.0 / std::pow(3.0, 1.0/3.0)) * std::pow(volume, 2.0/3.0) / (l01*l01 + l21*l21 +
					l02*l02 + l30*l30 + l31*l31 + l32*l32);
}
double
Tetrahedron::computeQualityWithMetric(const Eigen::Matrix3d& m0, const Eigen::Matrix3d& m1,
					const Eigen::Matrix3d& m2, const Eigen::Matrix3d& m3) const
{
	const math::Point p0 = m_pnts[0];
	const math::Point p1 = m_pnts[1];
	const math::Point p2 = m_pnts[2];
	const math::Point p3 = m_pnts[3];

	const math::Vector3d v01_ = p1-p0;
	const math::Vector3d v21_ = p2-p1;
	const math::Vector3d v02_ = p2-p0;
	const math::Vector3d v30_ = p3-p0;
	const math::Vector3d v31_ = p3-p1;
	const math::Vector3d v32_ = p3-p2;

	const Eigen::Vector3d v01 = Eigen::Vector3d(v01_.X(), v01_.Y(), v01_.Z());
	const Eigen::Vector3d v21 = Eigen::Vector3d(v21_.X(), v21_.Y(), v21_.Z());
	const Eigen::Vector3d v02 = Eigen::Vector3d(v02_.X(), v02_.Y(), v02_.Z());
	const Eigen::Vector3d v30 = Eigen::Vector3d(v30_.X(), v30_.Y(), v30_.Z());
	const Eigen::Vector3d v31 = Eigen::Vector3d(v31_.X(), v31_.Y(), v31_.Z());
	const Eigen::Vector3d v32 = Eigen::Vector3d(v32_.X(), v32_.Y(), v32_.Z());

	const double l01_m = 0.5*sqrt(v01.transpose() * (m0*v01)) + 0.5*sqrt(v01.transpose() * (m1*v01));
	const double l21_m = 0.5*sqrt(v21.transpose() * (m1*v21)) + 0.5*sqrt(v21.transpose() * (m2*v21));
	const double l02_m = 0.5*sqrt(v02.transpose() * (m0*v02)) + 0.5*sqrt(v02.transpose() * (m2*v02));
	const double l30_m = 0.5*sqrt(v30.transpose() * (m0*v30)) + 0.5*sqrt(v30.transpose() * (m3*v30));
	const double l31_m = 0.5*sqrt(v31.transpose() * (m1*v31)) + 0.5*sqrt(v31.transpose() * (m3*v31));
	const double l32_m = 0.5*sqrt(v32.transpose() * (m2*v32)) + 0.5*sqrt(v32.transpose() * (m3*v32));

	//I use uniform interpolation to compute the volume in the metric see aulauzet summer school
	double volume = std::abs(getVolume());
	double vol_m = volume * (sqrt(m0.determinant()) + sqrt(m1.determinant()) + sqrt(m2.determinant()) + sqrt(m3.determinant())) * 0.25;

	return (36.0 / std::pow(3.0, 1.0/3.0)) * std::pow(vol_m, 2.0/3.0) / (l01_m*l01_m+ l21_m*l21_m +
					l02_m*l02_m + l30_m*l30_m + l31_m*l31_m + l32_m*l32_m);
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
