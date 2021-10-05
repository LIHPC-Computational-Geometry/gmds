/*----------------------------------------------------------------------------*/
/*
 * Quadrilateral.h
 *
 *  Created on: 26 oct. 2015
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_QUADRILATERAL_H_
#define GMDS_MATH_QUADRILATERAL_H_
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
/** \class Quadrilateral
 *  \brief Defines a quadrilateral
 */
/*----------------------------------------------------------------------------*/
	class GMDSMath_API Quadrilateral
{

public:

	/*------------------------------------------------------------------------*/
	/** \brief  constructor
	 *
	 */
	Quadrilateral();

	/*------------------------------------------------------------------------*/
	/** \brief  constructor
	 *
	 * \param AP1 a point of the quad
	 * \param AP2 a point of the quad
	 * \param AP3 a point of the quad
	 * \param AP4 a point of the quad
	 *
	 */
	Quadrilateral(const Point& AP1, const Point& AP2, const Point& AP3, const Point& AP4);

	/*------------------------------------------------------------------------*/
	/** \brief  constructor
	 *
	 * \param AQuad a quad
	 *
	 */
	Quadrilateral(const Quadrilateral& AQuad);

	/*------------------------------------------------------------------------*/
	/** \brief  destructor
	 */
	virtual ~Quadrilateral();

	/*------------------------------------------------------------------------*/
	/** \brief  Getter for the quadrilateral point
	 *
	 * \param AIndex an integer
	 *
	 * \return a point
	 */
	const Point& getPoint(const TInt& AIndex) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Return the number of points
         *
         * \return the number of points
         */
        int getNbPoints() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Computes the center of the quadrilateral
         *
         * \return the center of the quadrilateral
         */
        Point getCenter() const;

		/*------------------------------------------------------------------------*/
		/** \brief  Computes the normal of the quadrilateral
         *
         * \return the normal of the quadrilateral
         */
		Vector3d getNormal() const;

		/*------------------------------------------------------------------------*/
		/** \brief  Computes the area of the quadrilateral
         */
		double area() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the scaled jacobian of the quadrilateral
         *
         * \return the scaled jacobian
         */
        double computeScaledJacobian2D() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Compute the scaled jacobian of the quadrilateral at one corner.
         *
         * \return the scaled jacobian
         */
        double computeScaledJacobian2DAt(int AVert) const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the scaled jacobian of the quadrilateral,
	 *          normalized between [-1., 1.]
         *
         * \return the scaled jacobian
         */
        double computeNormalizedScaledJacobian2D() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the mean edge length of the quadrilateral.
         *
         * \return the mean edge length
         */
        double computeMeanEdgeLength() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Overloaded operator<< for output
         */
        friend std::ostream& operator<<(std::ostream&, const Quadrilateral&);


protected:

	/*------------------------------------------------------------------------*/
        /** \brief  Return the jacobian matrix at a vertex.
         *
         * \return the jacobian matrix
         */
        math::Matrix<2,2,double> jacobian2D(const int iVertex) const;

	Point m_pnts[4];

};
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_QUADRILATERAL_H_ */
/*----------------------------------------------------------------------------*/
