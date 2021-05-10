/*----------------------------------------------------------------------------*/
/*
 * Line.h
 *
 *  Created on: sept. 3, 2015
 *      Author: franck ledoux
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_LINE_H_
#define GMDS_MATH_LINE_H_
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
class Segment;
class Plane;
class Triangle;
        /*----------------------------------------------------------------------------*/
        /** \class Line
         *  \brief Geometrical class implementing a line defining by 2 points
         */
        /*----------------------------------------------------------------------------*/
        class EXPORT_GMDS Line
        {
        public:
            
            /*------------------------------------------------------------------------*/
            /** \brief  constructor.
             *
             * \param AP1 first point of the line
             * \param AP2 second point of the line
             */
            Line(const Point& AP1, const Point& AP2):m_p1(AP1),m_p2(AP2) {;}
           
	    /*------------------------------------------------------------------------*/
            /** \brief  constructor.
             *
             * \param AP a point on the line
             * \param AVec direction vector of the line
             */
            Line(const Point& AP, const Vector3d& AVec):m_p1(AP),m_p2(AP+AVec) {;}
 
            /*------------------------------------------------------------------------*/
            /** \brief Copy constructor.
             *
             * \param AL another line
             */
            Line(const Line& AL):m_p1(AL.m_p1),m_p2(AL.m_p2) {;}
            
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator=
             */
            virtual Line& operator= (const Line&);
            
            /*------------------------------------------------------------------------*/
            /** \brief  Getter for the first point
             *
             * \return a point
             */
            const Point& getFirstPoint() const;
            
            /*------------------------------------------------------------------------*/
            /** \brief  Getter for the first point
             *
             * \return a point
             */
            const Point& getSecondPoint() const;
            
            /*------------------------------------------------------------------------*/
            /** \brief 2D Intersection with another line
             *
             * \param[IN]  AL the Line we want to get the intersection with
             * \param[OUT] AP the intersection point if it exists
             * \param[OUT] AParam the a parameter such that AP = a AS[0]+ (1-a) AS[1]
             *
             * \return true if (*this) intersects AL
             */
            bool intersect2D(const Line& AL, Point& AP, double& AParam) const;
           
	    /*------------------------------------------------------------------------*/
            /** \brief 3D Intersection with a plane
             *
             * \param[IN]  APlane the line we want to get the intersection with
             * \param[OUT] AP the intersection point if it exists
             * \param[OUT] AParam the a parameter such that AP = a AS[0]+ (1-a) AS[1]
             *
             * \return true if (*this) intersects APlane
             */
            bool intersect3D(const Plane& APlane, Point& AP, double& AParam) const;

	    /*------------------------------------------------------------------------*/
            /** \brief 3D Intersection with a triangle
             *
             * \param[IN]  ATri the triangle we want to get the intersection with
             * \param[OUT] AP the intersection point if it exists
             * \param[OUT] AParam the a parameter such that AP = a AS[0]+ (1-a) AS[1]
             *
             * \return true if (*this) intersects ATri
             */
            bool intersect3D(const Triangle& ATri, Point& AP, double& AParam) const;
 
            /*------------------------------------------------------------------------*/
            /** \brief 2D distance to a segment. The distance is the minimal orthogonal
             *         distance to @AS
             *
             * \param[IN]  AS  the segment we want to get the distance to
             * \param[OUT] AP1 the point of *this where the min. distance is computed
             * \param[OUT] AP2 the point of AS where the minimal distance is computed
             *
             * \return true if (*this) intersects AS
             */
            TCoord distance2D(const Segment& AS, Point& AP1, Point& AP2) const;

            /*------------------------------------------------------------------------*/
            /** \brief Compute the orthogonal projection of AP on (*this)
             *
             * \param[IN]  AP the point we want to project onto (*this)
             *
             * \return the projected point
             */
            Point project(const Point& AP) const;
            
        private:
            
            /* first point */
            Point m_p1;
            /* second point */
            Point m_p2;
            
        };
        /*----------------------------------------------------------------------------*/
    } // namespace math
    /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_RAY_H_ */
/*----------------------------------------------------------------------------*/
