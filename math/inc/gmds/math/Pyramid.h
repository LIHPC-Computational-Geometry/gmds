/*----------------------------------------------------------------------------*/
/*
 * Pyramid.h
 *
 *  Created on: 27 nov 2014
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_PYRAMID_H_
#define GMDS_MATH_PYRAMID_H_
/*----------------------------------------------------------------------------*/
#include <iostream>
/*----------------------------------------------------------------------------*/
// gmds file headers
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Point.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
class Triangle;
/*----------------------------------------------------------------------------*/
/** \class Pyramid
 *  \brief Defines a pyramid
 */
/*----------------------------------------------------------------------------*/
	class EXPORT_GMDS Pyramid
{

public:

	/*------------------------------------------------------------------------*/
	/** \brief  constructor
	 *
	 */
	Pyramid();

	/*------------------------------------------------------------------------*/
	/** \brief  constructor
  	 *
  	 * \param AP0 a point of the tet
  	 * \param AP1 a point of the tet
  	 * \param AP2 a point of the tet
  	 * \param AP3 a point of the tet
  	 *				 4
  	 *
         *
         *                          3 ----------- 2
         *                         /             /
         *                        /             /
         *                       0 ----------- 1
         *                                
	 */
	Pyramid(const Point& AP0, const Point& AP1, const Point& AP2, const Point& AP3, const Point& AP4);

        /*------------------------------------------------------------------------*/
        /** \brief  constructor
         *
         * \param APoints an array of points
         *
         */
        Pyramid(Point APoints[5]);

	/*------------------------------------------------------------------------*/
        /** \brief  constructor
         *
         * \param APoints a vector of points
         *
         */
        Pyramid(const std::vector<Point>& APoints);

	/*------------------------------------------------------------------------*/
	/** \brief  constructor
	 *
	 * \param ATet a tet
	 *
	 */
	Pyramid(const Pyramid& APyr);

	/*------------------------------------------------------------------------*/
	/** \brief  destructor
	 */
	virtual ~Pyramid();

        /*------------------------------------------------------------------------*/
        /** \brief  operator=
         */
        void operator=(const Pyramid& APyr);

	/*------------------------------------------------------------------------*/
	/** \brief  Getter for the Pyramid point
	 *
	 * \param AIndex an integer
	 *
	 * \return a point
	 */
	const Point& getPoint(const TInt& AIndex) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Compute the center of the pyramid.
         *
         * \return a point
         */
        const Point getCenter() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Compute the signed volume of the pyramid. It is computed as 
	 *          the sum of 4 tetrahedra
	 *
	 * \return the signed volume of the pyramid
	 */
	TCoord getVolume() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Compute the scaled jacobian of the pyramid.
         *
         * \return the scaled jacobian
         */
        double computeScaledJacobian() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the scaled jacobian of the pyramid
	 *          normalized between [-1., 1.]
         *
         * \return the scaled jacobian
         */
        double computeNormalizedScaledJacobian() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the mean ratio of the pyramid. 
	 *          It is computed on the four base vertices only
         *
         * \return the mean ratio
         */
        double computeMeanRatio() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the mean edge length of the pyranid.
         *
         * \return the mean edge length
         */
        double computeMeanEdgeLength() const;

        /*------------------------------------------------------------------------*/
        /** \brief  predicate indicating if a triangle and a pyramid intersect each
 	 *          other.
         * \param AT a triangle
         */
        bool intersect(const Triangle& AT, const bool AProper = false) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Overloaded operator<< for output
         */
        friend std::ostream& operator<<(std::ostream&, const Pyramid&);


protected:
	Point m_pnts[5];

};
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_PYRAMID_H_ */
/*----------------------------------------------------------------------------*/
