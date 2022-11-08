//
// Created by rochec on 28/10/2022.
//
/*----------------------------------------------------------------------------*/
#ifndef GMDS_BEZIERSURFACE_H
#define GMDS_BEZIERSURFACE_H
/*----------------------------------------------------------------------------*/
#include <cmath>
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/utils/Array.h>
#include <gmds/math/Point.h>
#include "GMDSMath_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*--------------------------------------------------------------------------*/
namespace math {
/*------------------------------------------------------------------------*/
/** \class BezierSurface
     *  \brief Defines a Bezier surface in 3D.
 */
class GMDSMath_API BezierSurface {


 public:

	/*------------------------------------------------------------------------*/
	/** \brief Constructor of a bezier surface from an 2D array of ordered
	    * 	set of control points.
       *
       * \param[in] APts the set of control points to define the curve
	 */
	BezierSurface(const Array2D<math::Point> & APts);

	/*------------------------------------------------------------------------*/
	/** \brief Returns the point located on the parametric surface with:
       *
       * \param[in] Au the first parameter in [0,1]
       * \param[in] Av the second parameter in [0,1]
       *
       * \return The point located at (Au,Av) on the parametric surface.
	 */
	math::Point operator()(const double& Au, const double& Av) const;

	/*------------------------------------------------------------------------*/
	/** \brief Returns a set of point that discretize the surface in
	    * 	ANb_m x ANb_n segments in the parametric space.
       *
       * \param[in] ANb_m the number of segments in the first direction
       * \param[in] ANb_n the number of segments in the second direction
       *
       * \return The set of point discretizing the surface (*this).
	 */
	Array2D<math::Point> getDiscretization(const int ANb_m, const int ANb_n) const;

 private:
	/*------------------------------------------------------------------------*/
	/** \brief Returns the binomial coefficient C_Ak^An. This is the number of
	    *		part of Ak elements in a set of An elements.
       *
       * \param[in] An the number of elements in the set
       * \param[in] Ak
       *
       * \return The binomial coefficient C_Ak^An.
	 */
	double BinomialCoefficient(const int An, const int Ak) const;

	/*------------------------------------------------------------------------*/
	/** \brief Returns the value of the Bernstein's Ai-th polynomial of degree An
	    * 	for Au in [0,1].
       *
       * \param[in] An degree of the Bernstein's polynomial
       * \param[in] Ai Bernstein's i-th polynomial of degree n
       * \param[in] Au parameter in [0,1]
       *
       * \return Bernstein's Ai-th polynomial of degree An evaluated at Au
	 */
	double BernsteinPolynomial(const int An, const int Ai, const double Au) const;

 private:
	/** ordered 2D array of control points of the curve */
	Array2D<math::Point> m_control_points;
	/** polynomial degree of the surface in the first direction */
	int m_degree_m;
	/** polynomial degree of the surface in the second direction */
	int m_degree_n;

};
/*----------------------------------------------------------------------------*/
}
  /*----------------------------------------------------------------------------*/
}

#endif     // GMDS_BEZIERSURFACE_H
