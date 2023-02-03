
/*----------------------------------------------------------------------------*/
/*
 * Ray.cpp
 *
 */
/*----------------------------------------------------------------------------*/
#include <gmds/math/Ray.h>
/*----------------------------------------------------------------------------*/
#include <gmds/math/Segment.h>
#include <gmds/math/Line.h>
#include <gmds/math/Plane.h>
#include <gmds/math/Triangle.h>
using namespace std;
/*----------------------------------------------------------------------------*/
namespace gmds{
  /*--------------------------------------------------------------------------*/
  namespace math{
    /*------------------------------------------------------------------------*/
    const Point& Ray::getPoint() const
    {
      return m_pnt;
    }
    /*------------------------------------------------------------------------*/
    const Vector3d& Ray::getDir() const {
      return m_dir;
    }
    /*------------------------------------------------------------------------*/
    const Vector3d& Ray::getDirUnit() const {
      if(m_isDirUnit)
	return m_dirUnit;

      m_isDirUnit = true;
      m_dirUnit = m_dir.getNormalize();
      return m_dirUnit;
    }
    /*------------------------------------------------------------------------*/
    Ray& Ray::operator=(const Ray& ARay) = default;
      /*----------------------------------------------------------------------------*/
      bool Ray::intersect2D(const Segment& AS, Point& AP, double& AParam) const {
          
          math::Point  p1= AS.getPoint(0);
          math::Point  p2= AS.getPoint(1);
	  
          math::Point src_pnt = m_pnt;
	  //Point dir_pnt m_pnt + m_dirUnit; // here was the problem...not addition
          math::Point dir_pnt;
	  dir_pnt.X() = m_pnt.X() + m_dirUnit.X(); 
	  dir_pnt.Y() = m_pnt.Y() + m_dirUnit.Y(); 
	  dir_pnt.Z() = m_pnt.Z() + m_dirUnit.Z(); 
	  
          Vector3d v12= p2-p1;
          if(v12.isColinear(m_dirUnit)){
              return false; // No intersection
          }
          
          double x1 = p1.X();
          double y1 = p1.Y();
          double x2 = p2.X();
          double y2 = p2.Y();
          double xs = src_pnt.X();
          double ys = src_pnt.Y();
          double xd = dir_pnt.X();  	  
          double yd = dir_pnt.Y();
          
          double D = x1*(yd-ys) + x2*(ys-yd) + xd*(y2-y1) + xs* (y1-y2);
          double N = x1*(yd-ys) + xs*(y1-yd) + xd*(ys-y1);
          
          //param in [p1,p2]
          double s = N/D;
          
          N =  -(x1*(ys- y2) + x2*(y1-ys) + xs*(y2-y1));
          
          //param in [src_pnt,dir_pnt]
          double t = N/D;
          
          if ((0.0 <= s) && (s <= 1.0) && (0.0 < t)){
              AP = p1 + s*(p2-p1);
              AParam = 1-s;
              return true;
          }
          return false;
          
      }
   
       bool Ray::SecondMetIntersect2D(const Segment& AS, Point& AP, double& AParam,  double& tempEpsilon) const {
         
          Point  p1 = AS.getPoint(0);
          Point  p2 = AS.getPoint(1);
	            
          Point q1 = m_pnt;
          Point q2 = m_pnt + m_dirUnit;	  
		
          Vector3d v12=p2-p1;
          if(v12.isColinear(m_dirUnit)){
		  //cout<<"v12.isColinear(m_dirUnit)"<<endl;
              return false; // No intersection
          }
                      
          math::Point r(v12.X(),v12.Y(),v12.Z());
	/*  
	  math::Point s = q2 - q1;
	  
	  double rx = r.X();
	  double ry = r.Y();
	  double sx = s.X();
	  double sy = s.Y();
	  double D = rx*sy - ry*sx;
	  
	  
	  
	  if(D==0)
	    return false;
	  
	  math::Point q1p1 = q1 - p1;
	  
	  double N = q1p1.X()*sy - q1p1.Y()*sx;
	  
	  std::cout<<"N "<<N<<std::endl;
	  std::cout<<"D "<<D<<std::endl;
	  
	  	 
	  if(std::fabs(N)<math::Constants::EPSILON)
	    return false;
	  
	  	//  if((D>1.00-math::Constants::EPSILON))
	    //return false; 
	  
	  AParam = N / D;
	  std::cout<<"AParam "<<AParam<<std::endl;
	  
	  //  if((AParam>=-0.000000001) && (AParam<=1.000000000-math::Constants::EPSILON)){
	  if((AParam>=-0.000000001) && (AParam<=1.000000001)){
	   	AP = p1 + AParam * r;
		std::cout<<"true"<<std::endl;
	   	return true;
	  }        
	  */
	  
	  ////////////////////////////////////////////
	  
	 math::Vector3d v1=q1-p1;
	  math::Vector3d v2=p2-p1;
	  math::Vector3d d=q2-q1;
	  math::Vector3d v3={-d[1], d[0],0};
	  double NN = v1.dot(v3);
	  double DD = v2.dot(v3);
	  //cout<<"v1 "<<v1<<endl;
	  //cout<<"v2 "<<v2<<endl;
	  //cout<<"v3 "<<v3<<endl;
	  //cout<<"d "<<d<<endl;
	  //cout<<"NN "<<NN<<endl;
	  //cout<<"DD "<<DD<<endl;
	  if(DD==0){	    
	    	return false;
	  }
	  
	  AParam = NN/DD;
	  
	  math::Vector3d crossProd = v2.cross(v1);
	  
	  double t1 = (crossProd.norm())/DD;
	  if(crossProd[2]<0)
	    	t1 = - t1;
	  //cout<<"v2.cross(v1) "<<v2.cross(v1)<<endl;
	 // std::cout<<"t1 "<<t1<<std::endl;
	  //std::cout<<"AParam "<<AParam<<std::endl;
	  //cout<<"p1 + AParam*r "<<p1 + AParam*r<<endl;
	  //cout<<"q1 + t1*d "<<q1 + t1*d<<endl;
	  
	  if(t1<0)
	    	return false;
	  
	  if((AParam>=0)&&(AParam<=1)){
	  	AP = p1 + AParam*r;
		return true;
	  }
          
          return false;
      }
     
      /*----------------------------------------------------------------------------*/
      bool Ray::intersect3D(const Segment& AS, Point& AP, double& AParamSeg,
                            double& AParamRay) const {
          // see details in http://www.lucidarme.me/?p=1872
          
          Point  p1= AS.getPoint(0);
          Point  p2= AS.getPoint(1);
          Point p3 = m_pnt;
          Point p4 = m_pnt + m_dirUnit;

          Vector3d v13=p3-p1;
          Vector3d v12=p2-p1;
          Vector3d v34=p4-p3;

          double tol = 0.001;
          
          if(v12.isColinear(m_dirUnit)){
          //    std::cout<<"\t colinear"<<std::endl;
              return false; // No intersection
          }

          double coplanar = std::abs(v13.dot(v12.cross(v34)));

          if(coplanar<tol){//Coplanar
              Vector3d v12_c_v34 = v12.cross(v34);
              double t = (v13.cross(v34)).dot(v12_c_v34)/(v12_c_v34.dot(v12_c_v34));
              double s = (v13.cross(v12)).dot(v12_c_v34)/(v12_c_v34.dot(v12_c_v34));
//              std::cout<<"Param seg = "<<t<<std::endl;
//              std::cout<<"Param ray = "<<s<<std::endl;
              if(0<=t && t<=1 && 0<s){
                  AParamSeg = t;
                  AP = p1 +t*v12;
                  AParamRay =s;
                  
                  return true;
              }
          }
//          else
//              std::cout<<"\t not coplanar: "<<std::abs(v13.dot(v12.cross(v34)))<<std::endl;
          return false;
          
      }
      /*----------------------------------------------------------------------------*/
      bool Ray::intersect2D(const Ray& AR, Point& AP) const {
          
          Point  p1= AR.getPoint();
          Point  p2= p1 + AR.getDirUnit();
          
          Point src_pnt = m_pnt;
          Point dir_pnt = m_pnt + m_dirUnit;

          Vector3d v12=p2-p1;
      if(v12.isColinear(m_dirUnit)){
	return false; // No intersection
      }

      double x1 = p1.X();
      double y1 = p1.Y();
      double x2 = p2.X();
      double y2 = p2.Y();
      double xs = src_pnt.X();
      double ys = src_pnt.Y();
      double xd = dir_pnt.X();
      double yd = dir_pnt.Y();

      double D = x1*(yd-ys) + x2*(ys-yd) + xd*(y2-y1) + xs* (y1-y2);
      double N = x1*(yd-ys) + xs*(y1-yd) + xd*(ys-y1);

      //param in [p1,p2]
      double s = N/D;

      N =  -(x1*(ys- y2) + x2*(y1-ys) + xs*(y2-y1));

      //param in [src_pnt,dir_pnt]
      double t = N/D;

      if ((0.0 <= s)  && (0.0 <= t)){
	AP = p1 + s*(p2-p1);
	return true;
      }
      return false;

    }
/*----------------------------------------------------------------------------*/
bool Ray::intersect3D(const Plane& APlane, Point& AP) const
{
	// check whether the line intersects the plane
	Line ln(m_pnt,m_dirUnit);
	double param(0.);
	if(!ln.intersect3D(APlane,AP,param)) {
		return false;
	}
	
	// check whether the intersection point is on the ray or on the other side
	if(param < 0.) {
		return false;
	}
	
	return true;
}
/*----------------------------------------------------------------------------*/
bool Ray::intersect3D(const Triangle& ATri, Point& AP) const
{
        // check whether the line intersects the triangle
        Line ln(m_pnt,m_dirUnit);
        double param(0.);
        if(!ln.intersect3D(ATri,AP,param)) {
                return false;
        }

        // check whether the intersection point is on the ray or on the other side
        if(param < 0.) {
                return false;
        }

        return true;
}
/*----------------------------------------------------------------------------*/
Point
Ray::project(const Point& AP) const
{
	gmds::math::Line line(m_pnt,m_dirUnit);
	gmds::math::Point projectedPoint = line.project(AP);

	if(AP.distance2(projectedPoint) < AP.distance2(m_pnt)) {
		return projectedPoint;
	} else {
		return m_pnt;
	}
}
/*---------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& AStr, const Ray& ARay){
        AStr<<"("<<ARay.m_pnt<<" | "<<ARay.m_dirUnit<<")";
                return AStr;
}
    /*----------------------------------------------------------------------------*/
  } // namespace math
  /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
