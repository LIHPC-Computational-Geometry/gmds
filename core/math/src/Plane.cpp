/*----------------------------------------------------------------------------*/
/*
 * Plane.cpp
 *
 */
/*----------------------------------------------------------------------------*/
#include <gmds/math/Plane.h>
#include <gmds/math/Numerics.h>
/*----------------------------------------------------------------------------*/
#include <gmds/math/Segment.h>
#include <gmds/math/Triangle.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*----------------------------------------------------------------------------*/
    namespace math{
        /*----------------------------------------------------------------------------*/
        Plane::Plane()
        = default;
        /*----------------------------------------------------------------------------*/
        Plane::Plane(const Point&AP, const Vector3d& AN)
        : m_pnt(AP),m_normal(AN)
        {
            m_normal.normalize();
        }
        /*----------------------------------------------------------------------------*/
        Plane::Plane(const Point&AP1, const Point &AP2, const Point &AP3)
        : m_pnt(AP1)
        {
            Vector3d v1= AP2-AP1;
            Vector3d  v2=AP3-AP1;
            m_normal = v1.cross(v2);
            m_normal.normalize();
        }
        /*----------------------------------------------------------------------------*/
        Plane::Plane(const Triangle& AT)
        {
            
            m_pnt = AT.getPoint(0);
            Point p1 = AT.getPoint(1);
            Point p2 = AT.getPoint(2);
            Vector3d v1=p1-m_pnt;
            Vector3d v2=p2-m_pnt;
            m_normal = v1.cross(v2);
            m_normal.normalize();
        }
        /*----------------------------------------------------------------------------*/
        const Point& Plane::getPoint() const
        {
            return m_pnt;
        }
        /*----------------------------------------------------------------------------*/
        const Vector3d& Plane::getNormal() const {
            return m_normal;
        }
        /*----------------------------------------------------------------------------*/
        Plane& Plane::operator=(const Plane& APlane) {
            if (APlane == *this) {
                return *this;
            }
            set(APlane.getPoint(), APlane.getNormal());
            return *this;
        }
        /*----------------------------------------------------------------------------*/
        bool Plane::operator==(const Plane& AP) const {
            if (&AP == this) {
                return true;
            }
            return (AP.getPoint()  == getPoint() &&
                    AP.getNormal() == getNormal() );
        }
        /*----------------------------------------------------------------------------*/
        bool Plane::operator!=(const Plane& AP) const {
            if (&AP == this) {
                return false;
            }
            return AP.getPoint() != getPoint() ||
            (AP.getNormal() != getNormal());
        }
        /*----------------------------------------------------------------------------*/
        bool  Plane::isIn(const Point& AP) const
        {
            return isZero(m_normal.dot(AP-m_pnt));
        }
        /*----------------------------------------------------------------------------*/
        bool Plane::isStrictlyOnLeft(const Point& AP) const
        {
            return (m_normal.dot(AP-m_pnt) < 0.0);
        }
        /*----------------------------------------------------------------------------*/
        TCoord Plane::distance(const Point& AP) const
        {
            Vector3d v= AP-m_pnt;
            return fabs(m_normal.dot(v));
        }
        
        /*----------------------------------------------------------------------------*/
        Point Plane::project(const Point& AP) const
        {
	        Vector3d v= AP-m_pnt;
            TCoord alpha = -m_normal.dot(v);		  
		  math::Point thePoint = AP;		  
		  
		  thePoint.X() = thePoint.X() + (alpha * m_normal)[0];
		  thePoint.Y() = thePoint.Y() + (alpha * m_normal)[1];
		  thePoint.Z() = thePoint.Z() + (alpha * m_normal)[2];		 
            //return AP+ (alpha * m_normal);
		  return thePoint;
        }
        
        /*----------------------------------------------------------------------------*/
        bool Plane::intersect(const Segment& AS, const bool AProper) const
        {
            if (isIn(AS.getPoint(0)) && isIn(AS.getPoint(1)))
                return !(AProper);
            
            if((this->isStrictlyOnLeft(AS.getPoint(0)) && !this->isStrictlyOnLeft(AS.getPoint(1))) ||
               (!this->isStrictlyOnLeft(AS.getPoint(0)) && this->isStrictlyOnLeft(AS.getPoint(1))))
                return true;

            return false;
        }
        /*----------------------------------------------------------------------------*/
        Plane::IntersectionType
        Plane::intersect(const Segment& AS, Point &PI,
                         double& AW0,
                         double& AW1,const bool AProper) const
        {
            Plane::IntersectionType return_value = SEGMENT_MIDDLE;
            Vector3d OP1 ({AS.getPoint(0).X(), AS.getPoint(0).Y(), AS.getPoint(0).Z()});
            Vector3d P1P2= AS.getPoint(1)-AS.getPoint(0);

            TCoord a, b, c, d;
            getEquationCoeffs(a, b, c, d);

            TCoord t = (d - OP1.dot(m_normal)) / (P1P2.dot(m_normal));

            if (t < 0.0 || t > 1.0) {
                return_value=NO_INTERSECTION;
                AW0=-1;
                AW1=-1;
            } else if (isZero(t)) {
                PI = AS.getPoint(0);
                AW0 = 1;
                AW1 = 0;
                return_value=SEGMENT_FIRST_PNT;
            } else if (isZero(t-1.0)) {
                PI = AS.getPoint(1);
                AW0=0;
                AW1=1;
                return_value=SEGMENT_SECOND_PNT;
            } else {
                PI = AS.getPoint(0) + t * P1P2;
                AW0 = 1-t;
                AW1 = t;
            }
            return return_value;
        }
        /*----------------------------------------------------------------------------*/
        bool Plane::intersect(const Plane& AP, const bool AProper) const
        {
            //check parallelism
            if (isZero(m_normal.cross(AP.m_normal).norm2()))
            {
                //we have parallel planes
                if (isZero(m_normal.dot(AP.m_pnt-m_pnt)))
                    return !(AProper); // the same plane !)
                else
                    return false; //different planes
            }
            else // now planes intersect each others
                return true;
        }
        /*----------------------------------------------------------------------------*/
    } // namespace math
    /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
