
/*----------------------------------------------------------------------------*/
/*
 * Quadrilateral.cpp
 *
 *  Created on: 26 oct. 2015
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#include <gmds/math/Quadrilateral.h>
/*----------------------------------------------------------------------------*/
#include <gmds/math/Vector.h>
#include <gmds/math/Triangle.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
Quadrilateral::Quadrilateral()
{
	m_pnts[0] = Point(0,0,0);
	m_pnts[1] = Point(1,0,0);
	m_pnts[2] = Point(1,1,0);
	m_pnts[3] = Point(0,1,0);
}
/*----------------------------------------------------------------------------*/
Quadrilateral::Quadrilateral(const Point& AP1, const Point& AP2, const Point& AP3, const Point& AP4)
{
	m_pnts[0] = AP1;
	m_pnts[1] = AP2;
	m_pnts[2] = AP3;
	m_pnts[3] = AP4;
}
/*----------------------------------------------------------------------------*/
Quadrilateral::Quadrilateral(const Quadrilateral& AQ)
{
	m_pnts[0] = AQ.m_pnts[0];
	m_pnts[1] = AQ.m_pnts[1];
	m_pnts[2] = AQ.m_pnts[2];
	m_pnts[3] = AQ.m_pnts[3];
}
/*----------------------------------------------------------------------------*/
Quadrilateral::~Quadrilateral(){}
/*----------------------------------------------------------------------------*/
const Point& Quadrilateral::getPoint(const TInt& AIndex) const
{
	return m_pnts[AIndex];
}
/*----------------------------------------------------------------------------*/
int Quadrilateral::getNbPoints()
{
        return 4;
}
/*----------------------------------------------------------------------------*/
Point
Quadrilateral::getCenter() const {
	Point pt((m_pnts[0]+m_pnts[1]+m_pnts[2]+m_pnts[3])*(1./4.));
	return pt;
}
/*----------------------------------------------------------------------------*/
Vector3d
Quadrilateral::getNormal() const {
        Point center(this->getCenter());

        Triangle tri0(center, m_pnts[0], m_pnts[1]);
        Triangle tri1(center, m_pnts[1], m_pnts[0]);
        Triangle tri2(center, m_pnts[2], m_pnts[3]);
        Triangle tri3(center, m_pnts[0], m_pnts[2]);

        Vector3d normal = tri0.getNormal() + tri1.getNormal() + tri2.getNormal() + tri3.getNormal();

        return normal;
    }
/*----------------------------------------------------------------------------*/
double
Quadrilateral::area() const
    {
        Point center(this->getCenter());

        Triangle tri0(center, m_pnts[0], m_pnts[1]);
        Triangle tri1(center, m_pnts[1], m_pnts[2]);
        Triangle tri2(center, m_pnts[2], m_pnts[3]);
        Triangle tri3(center, m_pnts[3], m_pnts[0]);

        double area = tri0.area() + tri1.area() + tri2.area() + tri3.area();

        return area;
    }
/*----------------------------------------------------------------------------*/
double
Quadrilateral::computeScaledJacobian2D() const
{        
	int neighbors[4][2] =
        {       
                {1,3},
                {2,0},
                {3,1},
                {0,2}
        };
        
        double scaledJ[4]; 
        for (int iVertex = 0; iVertex < 4; iVertex++) {
                
                Matrix<2,2,double> A = this->jacobian2D(iVertex);
                
                int i0    = neighbors[iVertex][0];
                int i1    = neighbors[iVertex][1];
                double l0 = (m_pnts[iVertex]-m_pnts[i0]).norm2();
                double l1 = (m_pnts[iVertex]-m_pnts[i1]).norm2();
                double l01 = sqrt(l0*l1);
                double det = A.det();
                if (l01 == 0 || det == 0) 
                        scaledJ[iVertex] = 0.;
                else    
                        scaledJ[iVertex] = det/l01;
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
Quadrilateral::computeScaledJacobian2DAt(int AVert) const
    {
        int neighbors[4][2] =
                {
                        {1,3},
                        {2,0},
                        {3,1},
                        {0,2}
                };

        double scaledJmin;

        Matrix<2,2,double> A = this->jacobian2D(AVert);

        int i0    = neighbors[AVert][0];
        int i1    = neighbors[AVert][1];
        double l0 = (m_pnts[AVert]-m_pnts[i0]).norm2();
        double l1 = (m_pnts[AVert]-m_pnts[i1]).norm2();
        double l01 = sqrt(l0*l1);
        double det = A.det();
        if (l01 == 0 || det == 0)
            scaledJmin = 0.;
        else
            scaledJmin = det/l01;

        return scaledJmin;
    }
/*----------------------------------------------------------------------------*/
double
Quadrilateral::computeNormalizedScaledJacobian2D() const
{
        double scaledJ = this->computeScaledJacobian2D();

	if(scaledJ >  1.) scaledJ =  1.;
	if(scaledJ < -1.) scaledJ = -1.;

	return scaledJ;
}
/*----------------------------------------------------------------------------*/
double
Quadrilateral::computeMeanEdgeLength() const
{
  double sumLength = 0;
  sumLength += m_pnts[0].distance(m_pnts[1]);
  sumLength += m_pnts[1].distance(m_pnts[2]);
  sumLength += m_pnts[2].distance(m_pnts[3]);
  sumLength += m_pnts[3].distance(m_pnts[0]);

  sumLength /= 4.;  
  return sumLength;
}
/*---------------------------------------------------------------------------*/
math::Matrix<2,2,double>
Quadrilateral::jacobian2D(const int iVertex) const
{
	math::Matrix<2,2,double> mat;

        switch(iVertex) {
        case 0:
		     mat={m_pnts[1].X()-m_pnts[0].X(),
					 m_pnts[1].Y()-m_pnts[0].Y(),
                m_pnts[3].X()-m_pnts[0].X(),
                m_pnts[3].Y()-m_pnts[0].Y()};
                break;
        case 1:
		     mat={m_pnts[1].X()-m_pnts[0].X(),
                m_pnts[1].Y()-m_pnts[0].Y(),
                m_pnts[2].X()-m_pnts[1].X(),
                m_pnts[2].Y()-m_pnts[1].Y()};
                break;
        case 2:
		     mat={ m_pnts[2].X()-m_pnts[3].X(),
                 m_pnts[2].Y()-m_pnts[3].Y(),
                m_pnts[2].X()-m_pnts[1].X(),
                 m_pnts[2].Y()-m_pnts[1].Y()};
                break;
        case 3:
		     mat={m_pnts[2].X()-m_pnts[3].X(),
                m_pnts[2].Y()-m_pnts[3].Y(),
                m_pnts[3].X()-m_pnts[0].X(),
                 m_pnts[3].Y()-m_pnts[0].Y()};
                break;
	default:
                throw GMDSException("Quadrilateral::jacobian invalid vertex number.");
                break;
        }

        return mat;
}
/*----------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& AStr, const Quadrilateral& AQuad){
        AStr<<"Quadrilateral"<<std::endl;
        for(int iPoint=0; iPoint<4; iPoint++) {
                AStr<<std::endl<<AQuad.getPoint(iPoint);
        }
        return AStr;
}
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/

