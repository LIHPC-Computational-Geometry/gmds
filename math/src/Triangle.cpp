
/*----------------------------------------------------------------------------*/
/*
 * Triangle.cpp
 *
 *  Created on: 3 juil. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <gmds/math/Triangle.h>
/*----------------------------------------------------------------------------*/
#include <gmds/math/Vector.h>
#include <gmds/math/Segment.h>
#include <gmds/math/Plane.h>
#include <gmds/math/Ray.h>
#include <gmds/math/Numerics.h>
#include <gmds/math/Constants.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
Triangle::Triangle()
{
	m_pnts[0] = Point(0,0,0);
	m_pnts[1] = Point(1,0,0);
	m_pnts[2] = Point(0,1,0);
}
/*----------------------------------------------------------------------------*/
Triangle::Triangle(const Point& AP1, const Point& AP2, const Point& AP3)
{
	m_pnts[0] = AP1;
	m_pnts[1] = AP2;
	m_pnts[2] = AP3;
}
/*----------------------------------------------------------------------------*/
Triangle::Triangle(const Triangle& AT)
{
	m_pnts[0] = AT.m_pnts[0];
	m_pnts[1] = AT.m_pnts[1];
	m_pnts[2] = AT.m_pnts[2];
}
/*----------------------------------------------------------------------------*/
Triangle::~Triangle(){}
/*----------------------------------------------------------------------------*/
    const Point& Triangle::getPoint(const TInt& AIndex) const
    {
        return m_pnts[AIndex];
    }
/*----------------------------------------------------------------------------*/
    void Triangle::setPoint(const TInt& AIndex, Point& APnt)
    {
        m_pnts[AIndex]= APnt;
    }
    /*----------------------------------------------------------------------------*/
    int Triangle::getNbPoints() const
    {
        return 3;
    }
    
    /*----------------------------------------------------------------------------*/
    double Triangle::area() const
    {
        Vector3d v1 = m_pnts[1] - m_pnts[0];
        Vector3d v2 = m_pnts[2] - m_pnts[0];
        return 0.5 * (v1.cross(v2)).norm();
    }
    /*----------------------------------------------------------------------------*/
    double Triangle::signedArea() const
    {
        Vector3d v1 = m_pnts[1] - m_pnts[0];
        Vector3d v2 = m_pnts[2] - m_pnts[0];
        Vector3d cross_prod = v1.cross(v2);
        if (cross_prod.Z() >= 0)
          return 0.5 * cross_prod.norm();
        else
          return -0.5 * cross_prod.norm();
    }
    /*----------------------------------------------------------------------------*/
    double Triangle::angle() const
    {
	    Vector3d v1 = m_pnts[1] - m_pnts[0];
	    Vector3d v2 = m_pnts[2] - m_pnts[0];
        v1.normalize();
        v2.normalize();

        return std::acos(v1.dot(v2));
    }
/*----------------------------------------------------------------------------*/
bool Triangle::isGood() const {
	Vector3d v1 = m_pnts[1] - m_pnts[0];
	Vector3d v2 = m_pnts[2] - m_pnts[0];
	return !v1.isColinear(v2);
}
/*----------------------------------------------------------------------------*/
    Vector3d Triangle::getNormal() const {
	    Vector3d v1 = m_pnts[1] - m_pnts[0];
	    Vector3d v2 = m_pnts[2] - m_pnts[0];
        return v1.cross(v2);
}
/*----------------------------------------------------------------------------*/
Point
Triangle::getCenter() const {
	Point pt((m_pnts[0]+m_pnts[1]+m_pnts[2])*(1./3.));
	return pt;
}
/*----------------------------------------------------------------------------*/
Point Triangle::getCircumcenter() const
{
	math::Point A = m_pnts[0];
	math::Point B = m_pnts[1];
	math::Point C = m_pnts[2];
	double xA = A.X();
	double yA = A.Y();
	double xB = B.X();
	double yB = B.Y();
	double xC = C.X();
	double yC = C.Y();
	double S = signedArea();
	double xO = ((xA*xA+yA*yA)*(yB-yC)-(xB*xB+yB*yB)*(yA-yC)+(xC*xC+yC*yC)*(yA-yB))*(1./(4.*S));
	double yO = -((xA*xA+yA*yA)*(xB-xC)-(xB*xB+yB*yB)*(xA-xC)+(xC*xC+yC*yC)*(xA-xB))*(1./(4.*S));
	math::Point Ctr(xO, yO, 0.);
	return Ctr;
}
/*----------------------------------------------------------------------------*/
double
Triangle::computeScaledJacobian2D() const
{
  throw GMDSException("Triangle::computeScaledJacobian2D not implemented yet.");

  math::Vector3d crossProduct = (m_pnts[1]-m_pnts[0]).cross(m_pnts[2]-m_pnts[0]);
  
  
}
/*----------------------------------------------------------------------------*/
double
Triangle::computeNormalizedScaledJacobian2D() const
{
  throw GMDSException("Triangle::computeNormalizedScaledJacobian2D not implemented yet.");

  math::Vector3d crossProduct = (m_pnts[1]-m_pnts[0]).cross(m_pnts[2]-m_pnts[0]);
  
  
}
/*----------------------------------------------------------------------------*/
double
Triangle::computeMeanEdgeLength() const
{
  double sumLength = 0;
  sumLength += m_pnts[0].distance(m_pnts[1]);
  sumLength += m_pnts[1].distance(m_pnts[2]);
  sumLength += m_pnts[2].distance(m_pnts[0]);

  sumLength /= 3.;  
  return sumLength;
}
/*----------------------------------------------------------------------------*/
void
Triangle::computeBoundingBox(TCoord AMinXYZ[3], TCoord AMaxXYZ[3]) const
{
	AMinXYZ[0] = min3(m_pnts[0].X(),m_pnts[1].X(),m_pnts[2].X());
	AMinXYZ[1] = min3(m_pnts[0].Y(),m_pnts[1].Y(),m_pnts[2].Y());
	AMinXYZ[2] = min3(m_pnts[0].Z(),m_pnts[1].Z(),m_pnts[2].Z());
	AMaxXYZ[0] = max3(m_pnts[0].X(),m_pnts[1].X(),m_pnts[2].X());
        AMaxXYZ[1] = max3(m_pnts[0].Y(),m_pnts[1].Y(),m_pnts[2].Y());
        AMaxXYZ[2] = max3(m_pnts[0].Z(),m_pnts[1].Z(),m_pnts[2].Z());
}
/*----------------------------------------------------------------------------*/
bool
Triangle::intersect(const Triangle& ATri, const bool AProper) const
{
  Point P1[3] = {this->getPoint(0), this->getPoint(1), this->getPoint(2)};
  Point P2[3] = {ATri.getPoint(0), ATri.getPoint(1), ATri.getPoint(2)};
  
  Vector3d N2 = ATri.getNormal();
  TCoord D2 = -N2.dot(Vector3d({P2[0].X(),P2[0].Y(),P2[0].Z()}));

  TCoord dist_P10_to_Plane2 = N2.dot(Vector3d({P1[0].X(),P1[0].Y(),P1[0].Z()}))+D2;
  TCoord dist_P11_to_Plane2 = N2.dot(Vector3d({P1[1].X(),P1[1].Y(),P1[1].Z()}))+D2;
  TCoord dist_P12_to_Plane2 = N2.dot(Vector3d({P1[2].X(),P1[2].Y(),P1[2].Z()}))+D2;
  
  // Triangle T1 does not intersect the plane containing T2
  if ( dist_P10_to_Plane2>0.0 && dist_P11_to_Plane2>0.0 && dist_P12_to_Plane2>0.0) {
    return false;
  }
  if ( dist_P10_to_Plane2<0.0 && dist_P11_to_Plane2<0.0 && dist_P12_to_Plane2<0.0) {
    return false;
  }
  
  // We check if the two triangles are coplanar
  if ( dist_P10_to_Plane2==0.0 && dist_P11_to_Plane2==0.0 && dist_P12_to_Plane2==0.0) { //YES THEY ARE
    Vector3d normal = this->getNormal();
    int maxIndex = normal.getMaxAbsComponentIndex();
    if(maxIndex==2) {
      // we can project on plane Oxy
      Point p1(this->getPoint(0).X(),this->getPoint(0).Y());
      Point p2(this->getPoint(1).X(),this->getPoint(1).Y());
      Point p3(this->getPoint(2).X(),this->getPoint(2).Y());
      
      Point p4(ATri.getPoint(0).X(),ATri.getPoint(0).Y());
      Point p5(ATri.getPoint(1).X(),ATri.getPoint(1).Y());
      Point p6(ATri.getPoint(2).X(),ATri.getPoint(2).Y());

      Triangle t1(p1,p2,p3);
      Triangle t2(p4,p5,p6);
      return t1.intersect2D(t2,AProper);
    } else if(maxIndex==1) {
      // we can project on plane Oxz
      Point p1(this->getPoint(0).X(),this->getPoint(0).Z());
      Point p2(this->getPoint(1).X(),this->getPoint(1).Z());
      Point p3(this->getPoint(2).X(),this->getPoint(2).Z());
      
      Point p4(ATri.getPoint(0).X(),ATri.getPoint(0).Z());
      Point p5(ATri.getPoint(1).X(),ATri.getPoint(1).Z());
      Point p6(ATri.getPoint(2).X(),ATri.getPoint(2).Z());

      Triangle t1(p1,p2,p3);
      Triangle t2(p4,p5,p6);
      return t1.intersect2D(t2,AProper);
    } else {
      // we can project on plane Oyz
      Point p1(this->getPoint(0).Y(),this->getPoint(0).Z());
      Point p2(this->getPoint(1).Y(),this->getPoint(1).Z());
      Point p3(this->getPoint(2).Y(),this->getPoint(2).Z());
      
      Point p4(ATri.getPoint(0).Y(),ATri.getPoint(0).Z());
      Point p5(ATri.getPoint(1).Y(),ATri.getPoint(1).Z());
      Point p6(ATri.getPoint(2).Y(),ATri.getPoint(2).Z());

      Triangle t1(p1,p2,p3);
      Triangle t2(p4,p5,p6);
      return t1.intersect2D(t2,AProper);
    }

  }

  // check if T1 is in contact with the plane containing T2 but does not cross the plane
  // only one point in contact
  if ( dist_P10_to_Plane2==0.0 && ((dist_P11_to_Plane2>0.0 && dist_P12_to_Plane2>0.0) || (dist_P11_to_Plane2<0.0 && dist_P12_to_Plane2<0.0))) {
    if(!ATri.isIn(P1[0])) {
      return false;
    } else {
	return true;
    }
  }
  if ( dist_P11_to_Plane2==0.0 && ((dist_P10_to_Plane2>0.0 && dist_P12_to_Plane2>0.0) || (dist_P10_to_Plane2<0.0 && dist_P12_to_Plane2<0.0))) {
    if(!ATri.isIn(P1[1])) {
      return false;
    } else {
      return true;
    }
  }
  if ( dist_P12_to_Plane2==0.0 && ((dist_P11_to_Plane2>0.0 && dist_P10_to_Plane2>0.0) || (dist_P11_to_Plane2<0.0 && dist_P10_to_Plane2<0.0))) {
    if(!ATri.isIn(P1[2])) {
      return false;
    } else {
      return true;
    }
  }
  // segment in contact
  if ( dist_P10_to_Plane2==0.0 && dist_P11_to_Plane2==0.0) {
    Segment seg(P1[0],P1[1]);
    return ATri.intersect(seg,AProper);
  }
  if ( dist_P11_to_Plane2==0.0 && dist_P12_to_Plane2==0.0) {
    Segment seg(P1[1],P1[2]);
    return ATri.intersect(seg,AProper);
  }
  if ( dist_P10_to_Plane2==0.0 && dist_P12_to_Plane2==0.0) {
    Segment seg(P1[0],P1[2]);
    return ATri.intersect(seg,AProper);
  }

  // Now we check if T2  intersect the plane containing T1
  Vector3d N1 = this->getNormal();
  TCoord D1 = -N1.dot(Vector3d({P1[0].X(),P1[0].Y(),P1[0].Z()}));
  
  TCoord dist_P20_to_Plane1 = N1.dot(Vector3d({P2[0].X(),P2[0].Y(),P2[0].Z()}))+D1;
  TCoord dist_P21_to_Plane1 = N1.dot(Vector3d({P2[1].X(),P2[1].Y(),P2[1].Z()}))+D1;
  TCoord dist_P22_to_Plane1 = N1.dot(Vector3d({P2[2].X(),P2[2].Y(),P2[2].Z()}))+D1;
  
  // Triangle T2 does not intersect the plane containing T1
  if ( dist_P20_to_Plane1>0.0 && dist_P21_to_Plane1>0.0 && dist_P22_to_Plane1>0.0) {
    return false;
  }
  if ( dist_P20_to_Plane1<0.0 && dist_P21_to_Plane1<0.0 && dist_P22_to_Plane1<0.0) {
    return false;
  }
  
  // check if T2 is in contact with the plane containing T1 but does not cross the plane
  // only one point in contact
  if ( dist_P20_to_Plane1==0.0 && ((dist_P21_to_Plane1>0.0 && dist_P22_to_Plane1>0.0) || (dist_P21_to_Plane1<0.0 && dist_P22_to_Plane1<0.0))) {
    if(!this->isIn(P2[0])) {
      return false;
    } else {
      return true;
    }
  }
  if ( dist_P21_to_Plane1==0.0 && ((dist_P20_to_Plane1>0.0 && dist_P22_to_Plane1>0.0) || (dist_P20_to_Plane1<0.0 && dist_P22_to_Plane1<0.0))) {
    if(!this->isIn(P2[1])) {
      return false;
    } else {
      return true;
    }
  }
  if ( dist_P22_to_Plane1==0.0 && ((dist_P21_to_Plane1>0.0 && dist_P20_to_Plane1>0.0) || (dist_P21_to_Plane1<0.0 && dist_P20_to_Plane1<0.0))) {
    if(!this->isIn(P2[2])) {
      return false;
    } else {
      return true;
    }
  }
  // segment in contact
  if ( dist_P20_to_Plane1==0.0 && dist_P21_to_Plane1==0.0) {
    Segment seg(P2[0],P2[1]);
    return this->intersect(seg,AProper);
  }
  if ( dist_P21_to_Plane1==0.0 && dist_P22_to_Plane1==0.0) {
    Segment seg(P2[1],P2[2]);
    return this->intersect(seg,AProper);
  }
  if ( dist_P20_to_Plane1==0.0 && dist_P22_to_Plane1==0.0) {
    Segment seg(P2[0],P2[2]);
    return this->intersect(seg,AProper);
  }
  
  // Intersection line between the two planes that intersect for sure
    Vector3d D = N1.cross(N2);
  
  /* we project the vertices of T1 on the intersection line. In fact to
   * reduce the number of computations, we project it onto the line
   * parallel to the intersection line and going through the origin O.
   * As we want parameter and not projected point, it doesn't matter.
   */
  Point C0,C1,C2;
  TCoord d0,d1,d2;
  //C0 and C2 will be on the same side and C1 on the other side
  
  if ((dist_P10_to_Plane2 >= 0.0 && dist_P11_to_Plane2 >= 0.0) || (dist_P10_to_Plane2 <= 0.0 && dist_P11_to_Plane2 <= 0.0)) {
    C0 = P1[0]; C2 = P1[1]; C1 = P1[2];
    d0 = dist_P10_to_Plane2;
    d1 = dist_P12_to_Plane2;
    d2 = dist_P11_to_Plane2;
  } else if((dist_P10_to_Plane2 >= 0.0 && dist_P12_to_Plane2 >= 0.0) || (dist_P10_to_Plane2 <= 0.0 && dist_P12_to_Plane2 <= 0.0)) {
    C0 = P1[0]; C2 = P1[2]; C1 = P1[1];
    d0 = dist_P10_to_Plane2;
    d1 = dist_P11_to_Plane2;
    d2 = dist_P12_to_Plane2;
    
  } else {
    C0 = P1[1]; C2 = P1[2]; C1 = P1[0];
    d0 = dist_P11_to_Plane2;
    d1 = dist_P10_to_Plane2;
    d2 = dist_P12_to_Plane2;
  }
  
  TCoord proj_C0_L = D.dot({C0.X(),C0.Y(),C0.Z()});
  TCoord proj_C1_L = D.dot({C1.X(),C1.Y(),C1.Z()});
  TCoord proj_C2_L = D.dot({C2.X(),C2.Y(),C2.Z()});
  
  TCoord a = proj_C0_L + (proj_C1_L-proj_C0_L)*d0/(d0-d1);
  TCoord b = proj_C2_L + (proj_C1_L-proj_C2_L)*d2/(d2-d1);
  
  /* As previously, we project the vertices of T2 on the intersection line. */
  
  if ((dist_P20_to_Plane1 >= 0.0 && dist_P21_to_Plane1 >= 0.0) || (dist_P20_to_Plane1 <= 0.0 && dist_P21_to_Plane1 <= 0.0)) {
    C0 = P2[0]; C2 = P2[1]; C1 = P2[2];
    d0 = dist_P20_to_Plane1;
    d1 = dist_P22_to_Plane1;
    d2 = dist_P21_to_Plane1;
  }
  else if((dist_P20_to_Plane1 >= 0.0 && dist_P22_to_Plane1 >= 0.0) || (dist_P20_to_Plane1 <= 0.0 && dist_P22_to_Plane1 <= 0.0)) {
    C0 = P2[0]; C2 = P2[2]; C1 = P2[1];
    d0 = dist_P20_to_Plane1;
    d1 = dist_P21_to_Plane1;
    d2 = dist_P22_to_Plane1;
  } else {
    C0 = P2[1]; C2 = P2[2]; C1 = P2[0];
    d0 = dist_P21_to_Plane1;
    d1 = dist_P20_to_Plane1;
    d2 = dist_P22_to_Plane1;
  }
  
  proj_C0_L = D.dot({C0.X(),C0.Y(),C0.Z()});
  proj_C1_L = D.dot({C1.X(),C1.Y(),C1.Z()});
  proj_C2_L = D.dot({C2.X(),C2.Y(),C2.Z()});
  
  TCoord c = proj_C0_L + (proj_C1_L-proj_C0_L)*d0/(d0-d1);
  TCoord d = proj_C2_L + (proj_C1_L-proj_C2_L)*d2/(d2-d1);
  
  TCoord i1,i2,i3,i4;
  if(a<b) {
    i1=a;i2=b;
  } else {
    i1=b;i2=a;
  }
  
  if(c<d) {
    i3=c;i4=d;
  } else {
    i3=d;i4=c;
  }

  if(i3>i2 || i4<i1) {
    return false;
  } else if (AProper && (i3==i2 || i4==i1 || i1==i2 || i3==i4)) {
    //return GEOM_UNDEF;
    return true;
  } else {
    return true;
  }
}
/*----------------------------------------------------------------------------*/
bool
Triangle::intersect2D(const Triangle& ATri, const bool AProper) const
{
	Segment seg1(this->getPoint(0), this->getPoint(1));
	bool r1 = ATri.intersect2D(seg1, AProper);
	if (r1) {
		return r1;
	}
	Segment seg2(this->getPoint(1), this->getPoint(2));
	bool r2 = ATri.intersect2D(seg2, AProper);
	if (r2) {
		return r2;
	}
	Segment seg3(this->getPoint(2), this->getPoint(0));
	bool r3 = ATri.intersect2D(seg3, AProper);
	if (r3) {
		return r3;
	}
	if (!AProper) {
		return r3;
	}

	return false;
}
/*----------------------------------------------------------------------------*/
bool
Triangle::intersect(const Segment& AS, const bool AProper) const {
        // Do the segment line and the plane containing the triangle intersect
        // each other?
        Plane pl(*this);
        if (!pl.intersect(AS, AProper))
            return false;

	// we have no to check if the segment line lies onto the plane
	if (!AS.getPoint(0).areCoplanar(this->getPoint(0), this->getPoint(1), this->getPoint(2))
	         || !AS.getPoint(1).areCoplanar(this->getPoint(0), this->getPoint(1), this->getPoint(2))) {
		// no coplanar
		Point p;
        double w0=0, w1=0;

        pl.intersect(AS, p, w0, w1);

        return (AProper) ? this->isStrictlyIn(p) : this->isIn(p);
	} else {
		//coplanar

        Vector3d normal = this->getNormal();
		int maxIndex = normal.getMaxAbsComponentIndex();
		if (maxIndex == 2) {
			// we can project on plane Oxy
			Point p1(AS.getPoint(0).X(), AS.getPoint(0).Y());
			Point p2(AS.getPoint(1).X(), AS.getPoint(1).Y());
			Point a(this->getPoint(0).X(), this->getPoint(0).Y());
			Point b(this->getPoint(1).X(), this->getPoint(1).Y());
			Point c(this->getPoint(2).X(), this->getPoint(2).Y());
			Segment l2d(p1, p2);
			Triangle t2d(a, b, c);
			return t2d.intersect2D(l2d, AProper);
		} else if (maxIndex == 1) {
			// we can project on plane Oxz
			Point p1(AS.getPoint(0).X(), AS.getPoint(0).Z());
			Point p2(AS.getPoint(1).X(), AS.getPoint(1).Z());
			Point a(this->getPoint(0).X(), this->getPoint(0).Z());
			Point b(this->getPoint(1).X(), this->getPoint(1).Z());
			Point c(this->getPoint(2).X(), this->getPoint(2).Z());
			Segment l2d(p1, p2);
			Triangle t2d(a, b, c);
			return t2d.intersect2D(l2d, AProper);
		} else {
			// we can project on plane Oyz
			Point p1(AS.getPoint(0).Y(), AS.getPoint(0).Z());
			Point p2(AS.getPoint(1).Y(), AS.getPoint(1).Z());
			Point a(this->getPoint(0).Y(), this->getPoint(0).Z());
			Point b(this->getPoint(1).Y(), this->getPoint(1).Z());
			Point c(this->getPoint(2).Y(), this->getPoint(2).Z());
			Segment l2d(p1, p2);
			Triangle t2d(a, b, c);
			return t2d.intersect2D(l2d, AProper);
		}
	}
}
/*----------------------------------------------------------------------------*/
bool
Triangle::intersect2D(const Segment& ASeg, const bool AProper) const
{
	Segment s1(this->getPoint(0), this->getPoint(1));
	Segment s2(this->getPoint(1), this->getPoint(2));
	Segment s3(this->getPoint(2), this->getPoint(0));

	bool r1 = s1.intersect2D(ASeg, AProper);
	if (r1) {
		return r1;
	}
	bool r2 = s2.intersect2D(ASeg, AProper);
	if (r2) {
		return r2;
	}
	bool r3 = s3.intersect2D(ASeg, AProper);
	if (r3) {
		return r3;
	}

	return false;
}
/*----------------------------------------------------------------------------*/
bool
Triangle::intersect(const Ray& ARay, const bool AProper) const
{
	/* Check the case where the directional vector is null */
	if (ARay.getDir().isZero()) {
		// we have now to check if the point lies onto the plane
		if (!ARay.getPoint().areCoplanar(this->getPoint(0), this->getPoint(1), this->getPoint(2))) {
			return false;
		} else {
		// coplanar
			return (AProper) ? this->isStrictlyIn(ARay.getPoint()) : this->isIn(ARay.getPoint());
		}	
	}

	/* Check if the segment starting from the ray point, along the ray direction and of
  	 * length twice the sum of the distance between the ray point and the triangle points
  	 * intersects the triangle
  	 */
	TCoord maxDist = (this->getPoint(0)-ARay.getPoint()).norm();
	maxDist = maxDist + (this->getPoint(1)-ARay.getPoint()).norm();
	maxDist = maxDist + (this->getPoint(2)-ARay.getPoint()).norm();

	Segment seg(ARay.getPoint(), ARay.getPoint() + (maxDist * 2) * ARay.getDirUnit());
	return this->intersect(seg, AProper);
}
/*----------------------------------------------------------------------------*/
bool Triangle::isIn(const Point& AP) const
{
	if(!this->getPlaneIncluding().isIn(AP)) {
		return false;
	}
	
	TCoord coordX, coordY, coordZ;
	std::vector<TCoord> coordNew;
	std::vector<gmds::math::Point> AT(3);
	AT[0] = this->getPoint(0);
	AT[1] = this->getPoint(1);
	AT[2] = this->getPoint(2);
	Point::computeBarycentric(this->getPoint(0),this->getPoint(1),this->getPoint(2),AP,coordX,coordY,coordZ);
	
	if (coordX<0.0 || coordY<0.0 || coordZ<0.0) {
        return false;
     } 
     else {
        return true;
     }
        
	/* Point::computeBarycentric(AT,AP,coordNew);
	 * if (coordX<0.0 || coordY<0.0 || coordZ<0.0) {
	  	if (coordNew[0]>=0.0 && coordNew[1]>=0.0 && coordNew[2]>=0.0) 
		  throw GMDSException("Triangle::isIn computeBarycentric 2 functions 2 different results");
		return false;
	} 
	else {
	  	if (coordNew[0]<0.0 && coordNew[1]<0.0 && coordNew[2]<0.0) 
		  	throw GMDSException("Triangle::isIn computeBarycentric 2 functions 2 different results");
		return true;
	}*/

}
/*----------------------------------------------------------------------------*/
bool Triangle::isIn2ndMethod(const Point& AP) const
{
	if(!this->getPlaneIncluding().isIn(AP)) {
		return false;
	}
	
	TCoord coordX, coordY, coordZ;
	std::vector<TCoord> coordNew;
	std::vector<gmds::math::Point> AT(3);
	AT[0] = this->getPoint(0);
	AT[1] = this->getPoint(1);
	AT[2] = this->getPoint(2);
	
	Point::computeBarycentric2ndMethod(this->getPoint(0),this->getPoint(1),this->getPoint(2),AP,coordX,coordY,coordZ);
	
	if (coordX<0.0 || coordY<0.0 || coordZ<0.0) {
        return false;
     } 
     else {
        return true;
     }
     

}
/*----------------------------------------------------------------------------*/
bool
Triangle::isStrictlyIn(const Point& AP) const
{
        if(!this->getPlaneIncluding().isIn(AP)) {
                return false;
        }

        TCoord coordX, coordY, coordZ;
        Point::computeBarycentric(this->getPoint(0),this->getPoint(1),this->getPoint(2),AP,coordX,coordY,coordZ);
        if (coordX<=0.0 || coordY<=0.0 || coordZ<=0.0) {
        return false;
        } else {
        return true;
        }

}
/*----------------------------------------------------------------------------*/
Plane
Triangle::getPlaneIncluding() const
{
	return Plane(*this);
}
/*----------------------------------------------------------------------------*/
Point
Triangle::project(const Point& APoint) const
{
	Point X = this->getPlaneIncluding().project(APoint);
	
	TCoord coordX,coordY,coordZ;
	Point::computeBarycentric(this->getPoint(0),this->getPoint(1),this->getPoint(2),X,coordX,coordY,coordZ);

	if(coordX<0.0) {
		if(coordY<0.0) {
			return this->getPoint(2);
		} else if(coordZ<0.0) {
			return this->getPoint(1);
		} else {
			return Segment(this->getPoint(1),this->getPoint(2)).project(X);
		}
	} else if(coordY<0.0) {
		if(coordZ<0.0) {
			return this->getPoint(0);
		} else {
			return Segment(this->getPoint(0),this->getPoint(2)).project(X);
		}
	} else if(coordZ<0.0) {
		return Segment(this->getPoint(0),this->getPoint(1)).project(X);
	}

	// we are in the triangle
	return X;
}
/*----------------------------------------------------------------------------*/
TCoord
Triangle::distance2(const Point& APoint) const
{
        return project(APoint).distance2(APoint);
}
/*----------------------------------------------------------------------------*/
TCoord
Triangle::distance(const Point& APoint) const
{
        return sqrt(this->distance2(APoint));
}
/*---------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& AStr, const Triangle& ATri){
        AStr<<"Triangle"<<std::endl;
        for(int iPoint=0; iPoint<3; iPoint++) {
                AStr<<std::endl<<ATri.getPoint(iPoint);
        }
        return AStr;
}
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/

