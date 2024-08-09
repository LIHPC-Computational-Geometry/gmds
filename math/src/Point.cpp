/*----------------------------------------------------------------------------*/
/*
 * Point.cpp
 *
 *  Created on: 6 févr. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <gmds/math/Point.h>
#include <gmds/math/Matrix.h>
#include <gmds/math/Numerics.h>
/*-----------------------------------------------------------------------------*/
#include <cmath>
#include <gmds/math/Vector.h>
/*-----------------------------------------------------------------------------*/
namespace gmds{
    /*-------------------------------------------------------------------------*/
    namespace math{
        /*---------------------------------------------------------------------*/
        Point::Point(TCoord AX, TCoord AY, TCoord AZ)
        {
            m_coord[0]=AX;
            m_coord[1]=AY;
            m_coord[2]=AZ;
        }
        /*---------------------------------------------------------------------*/
        Point::~Point() = default;
        /*---------------------------------------------------------------------*/
        Point operator+(const Point& AP1, const Point& AP2){
            return Point(AP1.m_coord[0] + AP2.m_coord[0],
                    AP1.m_coord[1] + AP2.m_coord[1],
                    AP1.m_coord[2] + AP2.m_coord[2]);
        }
        /*---------------------------------------------------------------------*/
        bool Point::operator==(const Point& AP) const
        {
            if (&AP == this)
                return true;

            return AP.m_coord[0]==m_coord[0] && AP.m_coord[1]==m_coord[1] && AP.m_coord[2]==m_coord[2];
        }
        /*---------------------------------------------------------------------*/
        bool Point::operator<(const Point& AP) const
        {
            if (&AP == this)
                return false;

            return m_coord[0]<AP.m_coord[0] ||
                   (m_coord[0]==AP.m_coord[0] && m_coord[1]<AP.m_coord[1])||
                   (m_coord[0]==AP.m_coord[0] && m_coord[1]==AP.m_coord[1] &&m_coord[2]<AP.m_coord[2]);
        }
        /*---------------------------------------------------------------------*/
        bool Point::operator<=(const Point& AP) const
        {
            if (&AP == this)
                return true;

            return m_coord[0]<=AP.m_coord[0] ||
                   (m_coord[0]==AP.m_coord[0] && m_coord[1]<=AP.m_coord[1])||
                   (m_coord[0]==AP.m_coord[0] && m_coord[1]==AP.m_coord[1] &&m_coord[2]<=AP.m_coord[2]);
        }
        /*---------------------------------------------------------------------*/
        bool Point::operator!= (const Point& AP) const
        {
            if (&AP == this)
                return false;
            
            return AP.m_coord[0]!=m_coord[0] || AP.m_coord[1]!=m_coord[1] || AP.m_coord[2]!=m_coord[2];
        }
        /*---------------------------------------------------------------------*/
        TCoord Point::distance(const Point& AP) const{
            return sqrt(distance2(AP));
        }
        /*---------------------------------------------------------------------*/
        TCoord Point::distance2(const Point& AP) const{
            return (AP-(*this)).norm2();
        }
        /*---------------------------------------------------------------------*/
        bool Point::areColinear(const Point& AP2, const Point& AP3) const
        {
            Vector3d v = (AP2-(*this)).cross(AP3-(*this));
	         return (isZero(v[0]) && isZero(v[1]) && isZero(v[2]) );
        }
         /*---------------------------------------------------------------------*/
        bool Point::areColinear2ndMethod(const Point& AP2, const Point& AP3) const
        {
	        Vector3d v = (AP2-(*this)).cross(AP3-(*this));
	        return (isZero2ndMethod(v[0]) && isZero2ndMethod(v[1]) && isZero2ndMethod(v[2]) );
        }
        /*---------------------------------------------------------------------*/
        bool Point::
        areCoplanar(const Point& AP2, const Point& AP3, const Point& AP4) const{
            return isZero((AP2-(*this)).dot((AP3-(*this)).cross((AP4-(*this)))));
        }
        /*---------------------------------------------------------------------*/
        bool
        Point::isStrictlyOnLeft2D(const Point& AP1, const Point& AP2) const
        {

            double r = (AP2.X() - AP1.X()) * (this->Y() - AP1.Y()) - (this->X() - AP1.X()) * (AP2.Y() - AP1.Y());

            return (r > 0.0);
        }
        /*---------------------------------------------------------------------*/
        void Point::computeBarycentric2D(
                                         const math::Point& AT1,
                                         const math::Point& AT2,
                                         const math::Point& AT3,
                                         const math::Point& AP,
                                         TCoord& AX, TCoord& AY, TCoord& AZ)
        {
            TCoord det;

            if (AT1.areColinear(AT2,AT3))
            {
                throw GMDSException("flat triangle in barycentric computation");
            }
            TCoord x1 = AT1.X();
            TCoord y1 = AT1.Y();
            TCoord x2 = AT2.X();
            TCoord y2 = AT2.Y();
            TCoord x3 = AT3.X();
            TCoord y3 = AT3.Y();
            TCoord x = AP.X();
            TCoord y = AP.Y();

            Matrix<2, 2,double> T;
            T.set(0, 0, x1 - x3);
            T.set(0, 1, x2 - x3);
            T.set(1, 0, y1 - y3);
            T.set(1, 1, y2 - y3);

            det = T.det();

            if (isZero(det)) {
                throw GMDSException("error in barycentric computation");
            }

            TCoord a1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3));
            TCoord a2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3));
            //TCoord a3 = det - a1 - a2;
		  AX = a1 / det;
            AY = a2 / det;
		 
            AZ = 1.0 - AX - AY;
        }
	 /*---------------------------------------------------------------------*/
        void Point::computeBarycentric2D2ndMethod(
                                         const math::Point& AT1,
                                         const math::Point& AT2,
                                         const math::Point& AT3,
                                         const math::Point& AP,
                                         TCoord& AX, TCoord& AY, TCoord& AZ)
        {
            TCoord det;

            if (AT1.areColinear2ndMethod(AT2,AT3))
            {
                throw GMDSException("flat triangle in barycentric computation");
            }
            TCoord x1 = AT1.X();
            TCoord y1 = AT1.Y();
            TCoord x2 = AT2.X();
            TCoord y2 = AT2.Y();
            TCoord x3 = AT3.X();
            TCoord y3 = AT3.Y();
            TCoord x = AP.X();
            TCoord y = AP.Y();

            Matrix<2, 2,double> T;
            T.set(0, 0, x1 - x3);
            T.set(0, 1, x2 - x3);
            T.set(1, 0, y1 - y3);
            T.set(1, 1, y2 - y3);

            det = T.det();

            if (isZero2ndMethod(det)) {
                throw GMDSException("error in barycentric computation");
            }

            TCoord a1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3));
            TCoord a2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3));
            //TCoord a3 = det - a1 - a2;
		  AX = a1 / det;
            AY = a2 / det;
		 
            AZ = 1.0 - AX - AY;
        }
        /*---------------------------------------------------------------------*/
        void Point::computeBarycentric(const math::Point& AT0,
                                       const math::Point& AT1,
                                       const math::Point& AT2,
                                       const math::Point& AT3,
                                       const math::Point& AP,
                                       TCoord& A0, TCoord& A1,
                                       TCoord& A2, TCoord& A3) {


            VectorND<4, double>   b={AP.X(),AP.Y(),AP.Z(),1.0};
            Matrix<4, 4, double>  A={
	            AT0.X(), AT1.X(), AT2.X(), AT3.X(),
	            AT0.Y(), AT1.Y(), AT2.Y(), AT3.Y(),
	            AT0.Z(), AT1.Z(), AT2.Z(), AT3.Z(),
	            1.0    , 1.0    , 1.0    , 1.0    };
            //We solve AX=b
            VectorND<4, double> x = A.solve(b);

            A0=x[0]; A1=x[1]; A2=x[2]; A3=x[3];
        }

        /*---------------------------------------------------------------------*/
        void Point::computeBarycentric(const std::vector<math::Point>& AT,
                                       const math::Point& AP,
                                       std::vector<TCoord>& ACoeff) {
            if(AT.size()==4){
		         Matrix<4, 4, double>  A= {
                    AT[0].X(), AT[1].X(), AT[2].X(), AT[3].X(),
                    AT[0].Y(), AT[1].Y(), AT[2].Y(), AT[3].Y(),
                    AT[0].Z(), AT[1].Z(), AT[2].Z(), AT[3].Z(),
                    1.0      , 1.0      , 1.0      , 1.0      };

		         VectorND<4, double>   b = {AP.X(),AP.Y(),AP.Z(),1.0};
                //We solve AX=b
                VectorND<4, double> x = A.solve(b);
                ACoeff.resize(4);
                ACoeff[0]=x[0];
                ACoeff[1]=x[1];
                ACoeff[2]=x[2];
                ACoeff[3]=x[3];
            }
            else if (AT.size()==3){
		         Matrix<3, 3, double> A= {
                    AT[0].X(), AT[1].X(), AT[2].X(),
                    AT[0].Y(), AT[1].Y(), AT[2].Y(),
                    AT[0].Z(), AT[1].Z(), AT[2].Z()};

                VectorND<3, double>   b= {AP.X(),AP.Y(),AP.Z()};

                //We solve AX=b
                VectorND<3, double> x = A.solve(b);
                ACoeff.resize(3);
                ACoeff[0]=x[0];
                ACoeff[1]=x[1];
                ACoeff[2]=x[2];


            }
            else
                throw GMDSException("Point::computeBarycentric only for 3 or 4 points");
        }
        /*---------------------------------------------------------------------*/
        Point Point::massCenter(const std::vector<math::Point>& AP) {
            double x=0, y=0, z=0;
            for(auto p:AP){
                x+=p.X();
                y+=p.Y();
                z+=p.Z();
            }
            return Point(x/AP.size(),y/AP.size(),z/AP.size());
        }
        /*---------------------------------------------------------------------*/
        void Point::computeBarycentric(
                                       const math::Point& AT1,
                                       const math::Point& AT2,
                                       const math::Point& AT3,
                                       const math::Point& AP,
                                       TCoord& AX, TCoord& AY, TCoord& AZ)
        {
            const math::Point& p0 = AT1;
            const math::Point& p1 = AT2;
            const math::Point& p2 = AT3;
			//std::cout<<"p0.areColinear(p1,p2) "<<p0.areColinear(p1,p2)<<std::endl;
            if(p0.areColinear(p1,p2))
            {
                AX=-1.0;
                AY=-1.0;
                AZ=-1.0;
               return;
             //   throw GMDSException("flat triangle in barycentric computation");
            }
            if(!AP.areCoplanar(p0,p1,p2))
            {
                //throw GMDSException("Coplanarity is mandatory to compute barycentric coordinates of a point into a 3D triangle");
            }

            Vector3d v1= p1-p0;
            Vector3d v3= p2-p0;

            Vector3d normal = (v1.cross(v3));
            normal.normalize();

            int maxIndex = normal.getMaxAbsComponentIndex();
		//std::cout<<"maxIndex "<<maxIndex<<std::endl;
            if(maxIndex==2) {
                // we can project on plane Oxy
                Point p(AP.X() ,AP.Y());
                Point a(AT1.X(),AT1.Y());
                Point b(AT2.X(),AT2.Y());
                Point c(AT3.X(),AT3.Y());
                computeBarycentric2D(a,b,c,p,AX,AY,AZ);
            } else if(maxIndex==1) {
                // we can project on plane Oxz
                Point p(AP.X() ,AP.Z());
                Point a(AT1.X(),AT1.Z());
                Point b(AT2.X(),AT2.Z());
                Point c(AT3.X(),AT3.Z());
                computeBarycentric2D(a,b,c,p,AX,AY,AZ);
            } else {
                // we can project on plane Oyz
                Point p(AP.Y() ,AP.Z());
                Point a(AT1.Y(),AT1.Z());
                Point b(AT2.Y(),AT2.Z());
                Point c(AT3.Y(),AT3.Z());
                computeBarycentric2D(a,b,c,p,AX,AY,AZ);
            }
        }
        
         /*---------------------------------------------------------------------*/
        void Point::computeBarycentric2ndMethod(
                                       const math::Point& AT1,
                                       const math::Point& AT2,
                                       const math::Point& AT3,
                                       const math::Point& AP,
                                       TCoord& AX, TCoord& AY, TCoord& AZ)
        {
            const math::Point& p0 = AT1;
            const math::Point& p1 = AT2;
            const math::Point& p2 = AT3;
			//std::cout<<"p0.areColinear(p1,p2) "<<p0.areColinear(p1,p2)<<std::endl;
            if(p0.areColinear2ndMethod(p1,p2))
            {
                AX=-1.0;
                AY=-1.0;
                AZ=-1.0;
               return;
             //   throw GMDSException("flat triangle in barycentric computation");
            }
            if(!AP.areCoplanar(p0,p1,p2))
            {
                //throw GMDSException("Coplanarity is mandatory to compute barycentric coordinates of a point into a 3D triangle");
            }

            Vector3d v1= p1-p0;
            Vector3d v3= p2-p0;

            Vector3d normal = (v1.cross(v3));
            normal.normalize();

            int maxIndex = normal.getMaxAbsComponentIndex();
		//std::cout<<"maxIndex "<<maxIndex<<std::endl;
            if(maxIndex==2) {
                // we can project on plane Oxy
                Point p(AP.X() ,AP.Y());
                Point a(AT1.X(),AT1.Y());
                Point b(AT2.X(),AT2.Y());
                Point c(AT3.X(),AT3.Y());
                computeBarycentric2D2ndMethod(a,b,c,p,AX,AY,AZ);
            } else if(maxIndex==1) {
                // we can project on plane Oxz
                Point p(AP.X() ,AP.Z());
                Point a(AT1.X(),AT1.Z());
                Point b(AT2.X(),AT2.Z());
                Point c(AT3.X(),AT3.Z());
                computeBarycentric2D2ndMethod(a,b,c,p,AX,AY,AZ);
            } else {
                // we can project on plane Oyz
                Point p(AP.Y() ,AP.Z());
                Point a(AT1.Y(),AT1.Z());
                Point b(AT2.Y(),AT2.Z());
                Point c(AT3.Y(),AT3.Z());
                computeBarycentric2D2ndMethod(a,b,c,p,AX,AY,AZ);
            }
        }
        
        /*---------------------------------------------------------------------*/
        Point operator*(const TCoord& AK, const Point& AP){
            return Point(AK*AP.m_coord[0], AK*AP.m_coord[1], AK*AP.m_coord[2]);
        }
        /*---------------------------------------------------------------------*/
        Point operator*(const Point& AP, const TCoord& AK){
            return Point(AK*AP.m_coord[0], AK*AP.m_coord[1], AK*AP.m_coord[2]);
        }
        /*---------------------------------------------------------------------*/
        std::ostream& operator<<(std::ostream& AStr, const Point& AP){
            AStr<<"("<<AP.m_coord[0]<<", "<<AP.m_coord[1]<<", "<<AP.m_coord[2]<<")";
            return AStr;
        }
        
        /*---------------------------------------------------------------------*/
    } // namespace math
    /*-------------------------------------------------------------------------*/
} // namespace gmds
bool operator<(const gmds::math::Point& AP1, const gmds::math::Point& AP2){
    return ((AP1.X()<AP2.X()) ||
            (AP1.X()==AP2.X() && AP1.Y()<AP2.Y()) ||
            (AP1.X()==AP2.X() && AP1.Y()==AP2.Y() && AP1.Z()<AP2.Z()) );
}
/*----------------------------------------------------------------------------*/
