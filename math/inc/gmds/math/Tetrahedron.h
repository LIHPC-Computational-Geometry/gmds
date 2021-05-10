/*----------------------------------------------------------------------------*/
/*
 * Tetrahedron.h
 *
 *  Created on: 27 nov 2014
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_TETRAHEDRON_H_
#define GMDS_MATH_TETRAHEDRON_H_
/*----------------------------------------------------------------------------*/
#include <iostream>
/*----------------------------------------------------------------------------*/
// gmds file headers
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Point.h>
#include <gmds/math/Matrix.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
class Triangle;
/*----------------------------------------------------------------------------*/
/** \class Tetrahedron
 *  \brief Defines a tetrahedron
 */
/*----------------------------------------------------------------------------*/
	class EXPORT_GMDS Tetrahedron
{

public:

	/*------------------------------------------------------------------------*/
	/** \brief  constructor
	 *
	 */
	Tetrahedron();

	/*------------------------------------------------------------------------*/
	/** \brief  constructor
  	 *
  	 * \param AP0 a point of the tet
  	 * \param AP1 a point of the tet
  	 * \param AP2 a point of the tet
  	 * \param AP3 a point of the tet
  	 *
  	 *	               3__
  	 *                / \ \____
  	 *	        	  _/   \     2
  	 * 		        _/      \   /
  	 *		       /         \ /
  	 *		     0 -----------1
	 */
	Tetrahedron(const Point& AP0, const Point& AP1, const Point& AP2, const Point& AP3);

	/*------------------------------------------------------------------------*/
        /** \brief  constructor
         *
         * \param APoints an array of points
         *
         */
        Tetrahedron(Point APoints[4]);

	/*------------------------------------------------------------------------*/
        /** \brief  constructor
         *
         * \param APoints a vector of points
         *
         */
        Tetrahedron(const std::vector<Point>& APoints);

	/*------------------------------------------------------------------------*/
	/** \brief  constructor
	 *
	 * \param ATet a tet
	 *
	 */
	Tetrahedron(const Tetrahedron& ATet);

	/*------------------------------------------------------------------------*/
	/** \brief  destructor
	 */
	virtual ~Tetrahedron();

        /*------------------------------------------------------------------------*/
        /** \brief  operator=
         */
        void operator=(const Tetrahedron& ATet);


	/*------------------------------------------------------------------------*/
	/** \brief  Getter for the Tetrahedron point
	 *
	 * \param AIndex an integer
	 *
	 * \return a point
	 */
	const Point& getPoint(const TInt& AIndex) const;

    void setPoint(const TInt& AIndex, const Point& AP){m_pnts[AIndex]=AP;};
    /*------------------------------------------------------------------------*/
    /** \brief  Compute the center of the tetrahedron
     *
     * \return a point
     */
    const Point getCenter() const;

    /*------------------------------------------------------------------------*/
    /** \brief  Compute the signed volume of the tetrahedron
     *
     * \return the signed volume of the tetrahedron
     */
    TCoord getVolume() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Compute the scaled jacobian of the tetrahedron.
         *
         * \return the scaled jacobian
         */
        double computeScaledJacobian() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the scaled jacobian of the tetrahedron
	 *          normalized between [-1., 1.]
         *
         * \return the scaled jacobian
         */
        double computeNormalizedScaledJacobian() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the mean ratio of the tetrahedron.
         *
         * \return the mean ratio
         */
        double computeMeanRatio() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the mean edge length of the tetrahedron.
         *
         * \return the mean edge length
         */
        double computeMeanEdgeLength() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Computes the bounding box of the tetrahedron.
         *
         * \param AMinXYZ the lower front left coordinates
         * \param AMaxXYZ the upper back right coordinates
         */
        void computeBoundingBox(TCoord AMinXYZ[3], TCoord AMaxXYZ[3]) const;

        /*------------------------------------------------------------------------*/
        /** \brief  predicate indicating if a triangle and a tet intersect each
 	 *          other.
         * \param AT a triangle
         */
        bool intersect(const Triangle& AT, const bool AProper = false) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Overloaded operator<< for output
         */
        friend std::ostream& operator<<(std::ostream&, const Tetrahedron&);


protected:

        /*------------------------------------------------------------------------*/
        /** \brief  Return the jacobian matrix.
         *
         * \return the jacobian matrix
         */
        math::Matrix<3,3,double> jacobian() const;


	Point m_pnts[4];

};
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_TETRAHEDRON_H_ */
/*----------------------------------------------------------------------------*/
