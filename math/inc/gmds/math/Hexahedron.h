/*----------------------------------------------------------------------------*/
/*
 * Hexahedron.h
 *
 *  Created on: 16 oct 2014
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_HEXAHEDRON_H_
#define GMDS_MATH_HEXAHEDRON_H_
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <set>
/*----------------------------------------------------------------------------*/
// gmds file headers
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Point.h>
#include <gmds/math/Matrix.h>
#include "GMDSMath_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
class Triangle;
/*----------------------------------------------------------------------------*/
/** \class Hexahedron
 *  \brief Defines a hexahedron
 */
/*----------------------------------------------------------------------------*/
	class GMDSMath_API Hexahedron
{
        typedef enum {
            BOTTOM,
            TOP,
            FRONT,
            BACK,
            LEFT,
            RIGHT
        } EHexahedronFace;

public:

	/*------------------------------------------------------------------------*/
	/** \brief  constructor
	 *
	 */
	Hexahedron();

	/*------------------------------------------------------------------------*/
	/** \brief  constructor
  	 *
  	 * \param AP1 a point of the hex
  	 * \param AP2 a point of the hex
  	 * \param AP3 a point of the hex
  	 * \param AP4 a point of the hex
  	 * \param AP5 a point of the hex
  	 * \param AP6 a point of the hex
  	 * \param AP7 a point of the hex
  	 * \param AP8 a point of the hex
  	 *			    7--------------6
  	 * 			   /|		  /|
  	 * 			  / |	  	 / |
  	 * 			 4--------------5  |
  	 * 			 |  |		|  |
  	 * 			 |  3-----------|--2   
  	 * 			 | /  		| /   
  	 * 			 |/		|/
  	 * 			 0--------------1
  	 */
	Hexahedron(const Point& AP1, const Point& AP2, const Point& AP3, const Point& AP4,
		const Point& AP5, const Point& AP6, const Point& AP7, const Point& AP8);

        /*------------------------------------------------------------------------*/
        /** \brief  constructor
         *
         * \param APoints an array of points
         *
         */
	Hexahedron(Point APoints[8]);

	/*------------------------------------------------------------------------*/
        /** \brief  constructor
         *
         * \param APoints a vector of points
         *
         */
        Hexahedron(const std::vector<Point>& APoints);

	/*------------------------------------------------------------------------*/
	/** \brief  copy constructor
	 *
	 * \param AHex a hexahedron
	 *
	 */
	Hexahedron(const Hexahedron& AHex);

	/*------------------------------------------------------------------------*/
	/** \brief  destructor
	 */
	virtual ~Hexahedron();

        /*------------------------------------------------------------------------*/
        /** \brief  operator=
         */
        void operator=(const Hexahedron& AHex);


	/*------------------------------------------------------------------------*/
	/** \brief  Getter for the hexahedron point
	 *
	 * \param AIndex an integer
	 *
	 * \return a point
	 */
	const Point& getPoint(const TInt& AIndex) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Compute the center of the hexahedron
         *
         * \return a point
         */
        const Point getCenter() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Compute the signed volume of the hexahedron. It is computed as 
	 *          the sum of 6 pyramids.
	 *
	 * \return the signed volume of the hexahedron
	 */
	TCoord getVolume() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the scaled jacobian of the hexahedron.
         *
         * \return the scaled jacobian
         */
        double computeScaledJacobian() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Compute the scaled jacobian of the hexahedron at one corner.
         *
         * \return the scaled jacobian
         */
        double computeScaledJacobianAt(int AVert) const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the scaled jacobian of the hexahedron
	 *          normalized between [-1., 1.]
         *
         * \return the scaled jacobian
         */
        double computeNormalizedScaledJacobian() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the mean ratio of the hexahedron.
         *
         * \return the mean ratio
         */
        double computeMeanRatio() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the mean edge length of the hexahedron.
         *
         * \return the mean edge length
         */
        double computeMeanEdgeLength() const;

        /*------------------------------------------------------------------------*/
		/** \brief  Computes the bounding box of the hexahedron.
         *
         * \param AMinXYZ the lower front left coordinates
         * \param AMaxXYZ the upper back right coordinates
         */
		void computeBoundingBox(TCoord AMinXYZ[3], TCoord AMaxXYZ[3]) const;

        /*------------------------------------------------------------------------*/
        /** \brief  predicate indicating if a triangle and a hex intersect each
 	 	 *          other.
         * \param AT a triangle
         */
        bool intersect(const Triangle& AT, const bool AProper = false) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Overloaded operator<< for output
         */
        friend std::ostream& operator<<(std::ostream&, const Hexahedron&);

protected:
        /*------------------------------------------------------------------------*/
        /** \brief  Return the jacobian matrix at a vertex.
         *
         * \return the jacobian matrix
         */
    math::Matrix<3,3,double> jacobian(const int iVertex) const;
    
    
protected:
	Point m_pnts[8];

};
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_HEXAHEDRON_H_ */
/*----------------------------------------------------------------------------*/
