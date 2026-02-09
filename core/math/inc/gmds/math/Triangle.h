/*----------------------------------------------------------------------------*/
/*
 * Triangle.h
 *
 *  Created on: 3 juil. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_TRIANGLE_H_
#define GMDS_MATH_TRIANGLE_H_
/*----------------------------------------------------------------------------*/
// gmds file headers
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
#include "GMDSMath_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
class Segment;
class Ray;
class Plane;
/*----------------------------------------------------------------------------*/
/** \class Triangle
 *  \brief Defines a 3D point
 */
/*----------------------------------------------------------------------------*/
	class GMDSMath_API Triangle
{

public:

	/*------------------------------------------------------------------------*/
	/** \brief  constructor
	 *
	 */
	Triangle();

	/*------------------------------------------------------------------------*/
	/** \brief  constructor
	 *
	 * \param AP1 a point of the triangle
	 * \param AP2 a point of the triangle (after AP1)
	 * \param AP3 a point of the triangle (after AP2, before AP1)
	 *
	 */
	Triangle(const Point& AP1, const Point& AP2, const Point& AP3);

	/*------------------------------------------------------------------------*/
	/** \brief  constructor
	 *
	 * \param AT a triangle
	 *
	 */
	Triangle(const Triangle& AT);

	/*------------------------------------------------------------------------*/
	/** \brief  destructor
	 */
	virtual ~Triangle();

        /*------------------------------------------------------------------------*/
        /** \brief  Getter for the triangle point
         *
         * \param AIndex an integer
         *
         * \return a point
         */
        const Point& getPoint(const TInt& AIndex) const;
        /*------------------------------------------------------------------------*/
        /** \brief  Setter for the triangle point
         *
         * \param AIndex an integer
         * \param APnt the new location
         *
         */
        void setPoint(const TInt& AIndex, Point& APnt);
    /*------------------------------------------------------------------------*/
    /** \brief  Computes the area of the triangle
     */
    double area() const;

    /*------------------------------------------------------------------------*/
    /** \brief  Computes the signed area of the triangle (works for triangles in the plane)
     */
    double signedArea() const;

    /*------------------------------------------------------------------------*/
    /** \brief  Computes the angle (in rad) of the triangle, as seen by its first vertex
     */
    double angle() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Return the number of points
         *
         * \return the number of points
         */
        int getNbPoints() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Computes if the triangle is not degenerated (i.e. flat)
	 *
	 * \return 	false if the triangle points are aligned
	 * 			else true
	 */
	bool isGood() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Computes the normal of the triangle
         *
         * \return the normal of the triangle
         */
        Vector3d getNormal() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Computes the center of the triangle
         *
         * \return the center of the triangle
         */
        Point getCenter() const;

	     /*------------------------------------------------------------------------*/
	     /** \brief  Computes the center of the circumcircle triangle
         *
         * \return the center of the circumcircle of the triangle
	      */
	     Point getCircumcenter() const;

	/*------------------------------------------------------------------------*/ 
	/** \brief  Compute the scaled jacobian of the quadrilateral
         *
         * \return the scaled jacobian
         */
        double computeScaledJacobian2D() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the scaled jacobian of the quadrilateral,
	 *          normalized between [-1., 1.]
         *
         * \return the scaled jacobian
         */
        double computeNormalizedScaledJacobian2D() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the mean edge length of the triangle.
         *
         * \return the mean edge length
         */
        double computeMeanEdgeLength() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Computes the bounding box of the triangle.
         *
         * \param AMinXYZ the lower front left coordinates
         * \param AMaxXYZ the upper back right coordinates
         */
        void computeBoundingBox(TCoord AMinXYZ[3], TCoord AMaxXYZ[3]) const;

        /*------------------------------------------------------------------------*/
        /** \brief  predicate indicating if a triangle and a triangle intersect each
         *          other.
         * \param ATri a triangle
         */
        bool intersect(const Triangle& ATri, const bool AProper = false) const;

        /*------------------------------------------------------------------------*/
        /** \brief  predicate indicating if a triangle and a triangle intersect each
         *          other, in 2D.
         * \param ATri a triangle
         */
        bool intersect2D(const Triangle& ATri, const bool AProper = false) const;

        /*------------------------------------------------------------------------*/
        /** \brief  predicate indicating if a triangle and a segment intersect each
         *          other.
         * \param ASeg a segment
         */
        bool intersect(const Segment& ASeg, const bool AProper = false) const;

        /*------------------------------------------------------------------------*/
        /** \brief  predicate indicating if a triangle and a segment intersect each
         *          other, in 2D.
         * \param ASeg a segment
         */
        bool intersect2D(const Segment& ASeg, const bool AProper = false) const;

        /*------------------------------------------------------------------------*/
        /** \brief  predicate indicating if a triangle and a ray intersect each
         *          other.
         * \param ARay a ray
         */
        bool intersect(const Ray& ARay, const bool AProper = false) const;

	
        /*------------------------------------------------------------------------*/
        /** \brief  predicate indicating if a point is in a triangle
         *  
         * \param AP a point
         */
	bool isIn(const Point& AP) const;
	
	  /*------------------------------------------------------------------------*/
        /** \brief  predicate indicating if a point is in a triangle (account for numerical innacuracies with Epsilon)
         *  
         * \param AP a point
         */
	bool isIn2ndMethod(const Point& AP) const;

        /*------------------------------------------------------------------------*/
        /** \brief  predicate indicating if a point is strictly in a triangle
         *
         * \param AP a point
         */
	bool isStrictlyIn(const Point& AP) const;

	/*------------------------------------------------------------------------*/
        /** \brief  Return a plane built from the triangle
         *
         * \return a plane
         */
	Plane getPlaneIncluding() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Return the square distance between a point and a triangle
         *
         * \return a scalar
         */
        Point project(const Point& APoint) const;


        /*------------------------------------------------------------------------*/
        /** \brief  Return the square distance between a point and a triangle
         *
         * \return a scalar
         */
        TCoord distance2(const Point& APoint) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Return the distance between a point and a triangle
         *
         * \return a scalar
         */
        TCoord distance(const Point& APoint) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Overloaded operator<< for output
         */
        friend std::ostream& operator<<(std::ostream&, const Triangle&);


protected:
	Point m_pnts[3];

};
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_TRIANGLE_H_ */
/*----------------------------------------------------------------------------*/
