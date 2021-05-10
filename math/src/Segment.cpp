
/*----------------------------------------------------------------------------*/
/*
 * Segment.cpp
 */
/*----------------------------------------------------------------------------*/
#include <gmds/math/Segment.h>
/*----------------------------------------------------------------------------*/
#include <gmds/utils/Exception.h>
#include <gmds/math/Plane.h>
#include <gmds/math/Constants.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
void Segment::_reset()
{
	isVunit_ = false;
}
/*----------------------------------------------------------------------------*/
Segment::Segment()
{
	_reset();
}
/*----------------------------------------------------------------------------*/
Segment::Segment(const Point &AP1, const Point &AP2)
{
	pnts_[0] = AP1;
	pnts_[1] = AP2;
	_reset();
}
/*----------------------------------------------------------------------------*/
void Segment::set(const Point &AP1, const Point &AP2)
{
	pnts_[0] = AP1;
	pnts_[1] = AP2;
	_reset();
}
/*----------------------------------------------------------------------------*/
Segment& Segment::operator=(const Segment& ASeg)
{
	if (ASeg == *this)
		return *this;

	set(ASeg.getPoint(0), ASeg.getPoint(1));
	return *this;
}
/*----------------------------------------------------------------------------*/
Point Segment::getPoint(const int AIndex) const
{
	return pnts_[AIndex];
}
/*----------------------------------------------------------------------------*/
	Vector3d Segment::getDir() const
{
	return Vector3d(pnts_[0], pnts_[1]);
}
/*----------------------------------------------------------------------------*/
	Vector3d& Segment::getUnitVector() const
{
	if (isVunit_) {
		return vunit_;
	}
	isVunit_ = true;
	vunit_ = Vector3d(pnts_[0], pnts_[1]);
	vunit_.normalize();
	return vunit_;
}
/*----------------------------------------------------------------------------*/
Point Segment::computeCenter() const
{
	return (0.5 * (pnts_[0]+ pnts_[1]));
}
/*----------------------------------------------------------------------------*/
TCoord Segment::computeLength() const
{
	return pnts_[0].distance(pnts_[1]);
}
/*----------------------------------------------------------------------------*/
void Segment::computeBoundingBox(TCoord AMinXYZ[3], TCoord AMaxXYZ[3]) const
{
    AMinXYZ[0] = std::min(pnts_[0].X(), pnts_[1].X());
    AMinXYZ[1] = std::min(pnts_[0].Y(), pnts_[1].Y());
    AMinXYZ[2] = std::min(pnts_[0].Z(), pnts_[1].Z());
    AMaxXYZ[0] = std::max(pnts_[0].X(), pnts_[1].X());
    AMaxXYZ[1] = std::max(pnts_[0].Y(), pnts_[1].Y());
    AMaxXYZ[2] = std::max(pnts_[0].Z(), pnts_[1].Z());
}
/*----------------------------------------------------------------------------*/
TCoord Segment::computeLength2() const
{
	return pnts_[0].distance2(pnts_[1]);
}
/*----------------------------------------------------------------------------*/
bool Segment::operator==(const Segment& ASeg) const
{
	if (&ASeg == this)
		return true;

	return (ASeg.pnts_[0]== pnts_[0]) && (ASeg.pnts_[1] ==  pnts_[1]);
}
/*----------------------------------------------------------------------------*/
bool Segment::operator!= (const Segment& ASeg) const
{
	if (&ASeg == this)
		return false;

	return ASeg.pnts_[0] != pnts_[0] || ASeg.pnts_[1] != pnts_[1];
}
/*----------------------------------------------------------------------------*/
bool Segment::isIn(const Point& AP, const bool AProper) const
{
	if(!AP.areColinear(pnts_[0], pnts_[1])) {
		return false;
	}
	
	TCoord a = Vector3d(pnts_[0], AP).norm2();
	TCoord b = Vector3d(pnts_[1], AP).norm2();
	TCoord c = Vector3d(pnts_[0], pnts_[1]).norm2();
	
	return (a <= c) && (b <= c);
}
/*----------------------------------------------------------------------------*/
bool Segment::isIn2ndMethod(const Point& AP, const bool AProper) const
{
  	//std::cout<<"AP.areColinear(pnts_[0], pnts_[1]) "<<AP.areColinear(pnts_[0], pnts_[1])<<std::endl;
	if(!AP.areColinear2ndMethod(pnts_[0], pnts_[1])) {
		return false;
	}
	
	TCoord a = Vector3d(pnts_[0], AP).norm2();
	TCoord b = Vector3d(pnts_[1], AP).norm2();
	TCoord c = Vector3d(pnts_[0], pnts_[1]).norm2();
	//std::cout<<"a= "<<a<<" , b= "<<b<<" , c= "<<c<<" , c+gmds::math::Constants::EPSILON="<<c+gmds::math::Constants::EPSILON<<std::endl;
	return (a <= c+gmds::math::Constants::EPSILON && b <= c+gmds::math::Constants::EPSILON);
}
/*----------------------------------------------------------------------------*/
bool Segment::intersect(const Plane& AP, const bool AProper) const {
	return AP.intersect(*this, AProper);
}
/*----------------------------------------------------------------------------*/
bool Segment::
intersect(const Plane& AP, Point &API, const bool AProper) const {
	double w0=0, w1=0;
	return AP.intersect(*this, API, w0, w1, AProper);
}
/*----------------------------------------------------------------------------*/
TCoord Segment::distance(const Point& AP) const
{
	return AP.distance(project(AP));
}
/*----------------------------------------------------------------------------*/
TCoord Segment::distanceInf(const Segment& ASegment) const
{

	TCoord dist1 = ASegment.distance(project(this->getPoint(0)));
	TCoord dist2 = ASegment.distance(project(this->getPoint(1)));
	TCoord dist3 = this->distance(project(ASegment.getPoint(0)));
	TCoord dist4 = this->distance(project(ASegment.getPoint(1)));

	TCoord dist = 0.;
	if(dist<dist1)
		dist = dist1;
	if(dist<dist2)
		dist = dist2;
	if(dist<dist3)
		dist = dist3;
	if(dist<dist4)
		dist = dist4;

	return dist;
}
/*----------------------------------------------------------------------------*/
Point Segment::project(const Point& AP) const {
		Vector3d v1(pnts_[0], pnts_[1]);
		Vector3d v2(pnts_[0], AP);
	TCoord a = v1.dot(v2);
	if (a <= 0.0) {
		return pnts_[0];
	}
	TCoord b = v1.dot(v1);
	if (a >= b) {
		return pnts_[1];
	}
	TCoord ratio = a / b;
	return pnts_[0] + ratio * v1;
}
/*----------------------------------------------------------------------------*/
bool 
Segment::intersect2D(const Segment& ASeg, const bool AProper) const
{
	Point pA = this->getPoint(0);
	Point pB = this->getPoint(1);

	Point pC = ASeg.getPoint(0);
	Point pD = ASeg.getPoint(1);

	if((pA.X() < pC.X() && pB.X() < pC.X() && pA.X() < pD.X() && pB.X() < pD.X())
        || (pA.Y() < pC.Y() && pB.Y() < pC.Y() && pA.Y() < pD.Y() && pB.Y() < pD.Y())
        || (pA.X() > pC.X() && pB.X() > pC.X() && pA.X() > pD.X() && pB.X() > pD.X())
        || (pA.Y() > pC.Y() && pB.Y() > pC.Y() && pA.Y() > pD.Y() && pB.Y() > pD.Y())) {
                return false;
        }


	if (ASeg.isIn(pA) || ASeg.isIn(pB) || this->isIn(pC) || this->isIn(pD)) {
		return true;
	}
	bool left_A_CD = pA.isStrictlyOnLeft2D(pC, pD);
	bool left_B_CD = pB.isStrictlyOnLeft2D(pC, pD);
	bool left_C_AB = pC.isStrictlyOnLeft2D(pA, pB);
	bool left_D_AB = pD.isStrictlyOnLeft2D(pA, pB);

	if ((left_A_CD && left_B_CD) || (!left_A_CD && !left_B_CD)) {
		return false;
	}
	if ((left_C_AB && left_D_AB) || (!left_C_AB && !left_D_AB)) {
		return false;
	}
	return true;        
}
    /*----------------------------------------------------------------------------*/
    bool Segment::intersect3D(const Segment& AS, Point& AP, double& AParamSeg,
                          double& AParamThis) const {
        // see details in http://www.lucidarme.me/?p=1872
        AParamSeg=-1;
        AParamThis=-1;
        Point  p1= AS.getPoint(0);
        Point  p2= AS.getPoint(1);
        
        Point p3 = pnts_[0];
        Point p4 = pnts_[1];

		Vector3d v13(p1,p3);
		Vector3d v12(p1,p2);
		Vector3d v34(p3,p4);
        
        double tol = 0.01;
        
        if(v12.isColinear(v34)){
            //std::cout<<"colinear"<<std::endl;
            return false; // No intersection
        }
        
        double coplanar = std::abs(v13.dot(v12.cross(v34)));
        
        if(coplanar<tol){//Coplanar
			Vector3d v12_c_v34 = v12.cross(v34);
            double t = (v13.cross(v34)).dot(v12_c_v34)/(v12_c_v34.dot(v12_c_v34));
            double s = (v13.cross(v12)).dot(v12_c_v34)/(v12_c_v34.dot(v12_c_v34));
            if(0<=t && t<=1 && 0<=s && s<=1){
                //std::cout<<"\t t,s = "<<t<<", "<<s<<std::endl;
                AParamSeg = t;
                AParamThis = s;
                AP = p1 +t*v12;
                return true;
            }
        }
//        else
//            std::cout<<"\t not coplanar: "<<coplanar<<std::endl;
        return false;
        
    }

/*----------------------------------------------------------------------------*/
bool Segment::intersect2D(const Segment& ASeg, Point& APnt) const
{
  Point pA = this->getPoint(0);
  Point pB = this->getPoint(1);

  Point pC = ASeg.getPoint(0);
  Point pD = ASeg.getPoint(1);

  //OK We try to intersect so;

	Vector3d vAB(pA,pB);
	Vector3d vCD(pC,pD);
  if(vAB.isColinear(vCD)){
    return false; // No intersection
  }
  
  double x1 = pA.X();
  double y1 = pA.Y();
  double x2 = pB.X();
  double y2 = pB.Y();
  double xs = pC.X();
  double ys = pC.Y();
  double xd = pD.X();
  double yd = pD.Y();

  double D = x1*(yd-ys) + x2*(ys-yd) + xd*(y2-y1) + xs* (y1-y2);
  double N = x1*(yd-ys) + xs*(y1-yd) + xd*(ys-y1);

  //param in [p1,p2]
  double s = N/D;

  N =  -(x1*(ys- y2) + x2*(y1-ys) + xs*(y2-y1));

  //param in [src_pnt,dir_pnt]
  double t = N/D;

  if ((0.0 <= s) && (s <= 1.0) && (0.0 <= t) && (t <= 1.0)){
    APnt = pA + s*(pB-pA);
    return true;
  }
  return false;

  
}

/*----------------------------------------------------------------------------*/
bool Segment::SecondMetIntersect2D(const Segment& ASeg, Point& APnt, double& AParam, double& tempEpsilon) const
{
  	Point pA = this->getPoint(0);
  	Point pB = this->getPoint(1);

  	Point pC = ASeg.getPoint(0);
  	Point pD = ASeg.getPoint(1);
  

  	//OK We try to intersect so;

	Vector3d vAB(pA,pB);
	Vector3d vCD(pC,pD);
  	if(vAB.isColinear(vCD)){
    		return false; // No intersection
  	}
  
  	double x1 = pA.X();
  	double y1 = pA.Y();
  	double x2 = pB.X();
  	double y2 = pB.Y();
  	double xs = pC.X();
  	double ys = pC.Y();
  	double xd = pD.X();
  	double yd = pD.Y();

  	double D = x1*(yd-ys) + x2*(ys-yd) + xd*(y2-y1) + xs* (y1-y2);
  	double N = x1*(yd-ys) + xs*(y1-yd) + xd*(ys-y1);

  	//param in [p1,p2]

  	double t = N/D; 

  	N =  -(x1*(ys- y2) + x2*(y1-ys) + xs*(y2-y1));
	
	//double t = N/D; 
  	//param in [src_pnt,dir_pnt]
  	AParam = N/D;
		
    /*	std::cout<<-tempEpsilon<<" <= "<< AParam<<" = "<<(-tempEpsilon <= AParam)<<std::endl;
	
	std::cout<<AParam<<" <= "<<1.0+tempEpsilon<<" = "<<(AParam <= 1.0+tempEpsilon)<<std::endl;
	std::cout.precision(17);
	std::cout << "AParam " << std::fixed << AParam << std::endl;
	std::cout << "1+tempEpsilon " << std::fixed << 1+tempEpsilon << std::endl;
	
	std::cout<<-tempEpsilon<<" <= "<<t<<" = "<<(-tempEpsilon <= t)<<std::endl;
	std::cout<<t<<" <= "<<(1.0+tempEpsilon)<<" = "<<(t <= 1.0+tempEpsilon)<<std::endl;
    	*/
	// if ((-0.00001 <= AParam) && (AParam <= 1.00001) && (-0.00001 <= t) && (t <= 1.00001)){
   	if ((-tempEpsilon <= AParam) && (AParam <= 1+tempEpsilon) && (-tempEpsilon <= t) && (t <= 1+tempEpsilon)){
    		APnt = pA + t*(pB-pA);
    		return true;
  	}
  	return false;

  
}
/*---------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& AStr, const Segment& AS){
        AStr<<"("<<AS.pnts_[0]<<", "<<AS.pnts_[1]<<")";
                return AStr;
}
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
