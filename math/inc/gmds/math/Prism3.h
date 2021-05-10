/*----------------------------------------------------------------------------*/
/*
 * Prism3.h
 *
 *  Created on: 27 nov 2014
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_PRISM3_H_
#define GMDS_MATH_PRISM3_H_
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
/** \class Prism3
 *  \brief Defines a prism3
 */
/*----------------------------------------------------------------------------*/
	class EXPORT_GMDS Prism3
{

public:

	/*------------------------------------------------------------------------*/
	/** \brief  constructor
	 *
	 */
	Prism3();
                                            
	/*------------------------------------------------------------------------*/
	/** \brief  constructor
  	 *
  	 * \param AP0 a point of the prism
  	 * \param AP1 a point of the prism
  	 * \param AP2 a point of the prism
  	 * \param AP3 a point of the prism
         * \param AP4 a point of the prism
         * \param AP5 a point of the prism
         *
         * 
	 *		         ____ 5
	 *		     ___/    / 
	 *		   _/       / |
	 *		 3 __      /  |
	 *		     \__ 4    |                            
	 *		 |            |                       
	 *		 |       |    |                     
	 *		 |       |    |                                 
	 *		 |       | 
	 *		 |       |    2                     
	 *		 |       |   /
	 *		         |  /                     
	 *		 0 __      / 
	 *		     \__ 1        
         *
         *                                
	 */
	Prism3(const Point& AP0, const Point& AP1, const Point& AP2,
		const Point& AP3, const Point& AP4, const Point& AP5);

        /*------------------------------------------------------------------------*/
        /** \brief  constructor
         *
         * \param APoints an array of points
         *
         */
        Prism3(Point APoints[6]);

	/*------------------------------------------------------------------------*/
        /** \brief  constructor
         *
         * \param APoints a vector of points
         *
         */
        Prism3(const std::vector<Point>& APoints);

	/*------------------------------------------------------------------------*/
	/** \brief  constructor
	 *
	 * \param APrsm a prsm
	 *
	 */
	Prism3(const Prism3& APrsm);

	/*------------------------------------------------------------------------*/
	/** \brief  destructor
	 */
	virtual ~Prism3();

        /*------------------------------------------------------------------------*/
        /** \brief  operator=
         */
        void operator=(const Prism3& APrsm);

	/*------------------------------------------------------------------------*/
	/** \brief  Getter for the Prism3 point
	 *
	 * \param AIndex an integer
	 *
	 * \return a point
	 */
	const Point& getPoint(const TInt& AIndex) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Compute the center of the prism3
         *
         * \return a point
         */
        const Point getCenter() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Compute the signed volume of the prism3. It is computed as 
	 *          the sum of 2 tetrahedra and 3 pyramids.
	 *
	 * \return the signed volume of the prism3
	 */
	TCoord getVolume() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Compute the scaled jacobian of the prism3.
         *
         * \return the scaled jacobian
         */
        double computeScaledJacobian() const;
	/*------------------------------------------------------------------------*/
        /** \brief  Compute the scaled jacobian of the prism3
	 *          normalized between [-1., 1.]
         *
         * \return the scaled jacobian
         */
        double computeNormalizedScaledJacobian() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the mean ratio of the prism3.
         *
         * \return the mean ratio
         */
        double computeMeanRatio() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the mean edge length of the prism3.
         *
         * \return the mean edge length
         */
        double computeMeanEdgeLength() const;

        /*------------------------------------------------------------------------*/
        /** \brief  predicate indicating if a triangle and a prism3 intersect each
 	 *          other.
         * \param AT a triangle
         */
        bool intersect(const Triangle& AT, const bool AProper = false) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Overloaded operator<< for output
         */
        friend std::ostream& operator<<(std::ostream&, const Prism3&);


protected:

        /*------------------------------------------------------------------------*/
        /** \brief  Return the jacobian matrix at a vertex.
         *
         * \return the jacobian matrix
         */
        math::Matrix<3,3,double> jacobian(const int iVertex) const;


	Point m_pnts[6];

};
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_PRISM3_H_ */
/*----------------------------------------------------------------------------*/
