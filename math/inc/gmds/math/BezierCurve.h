/*----------------------------------------------------------------------------*/
/*
 * BezierCurve.h
 *
 *  Created on: 07/02/2015
 *      Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_BEZIERCURVE_H_
#define GMDS_MATH_BEZIERCURVE_H_
/*----------------------------------------------------------------------------*/
#include <cmath>
#include <vector>
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
  /*--------------------------------------------------------------------------*/
  namespace math {
    /*------------------------------------------------------------------------*/
    /** \class BezierCurve
     *  \brief Defines a Bezier curve in 3D. Underlying computations are based
     *         on the simple de Casteljau algorithm
     */
    class EXPORT_GMDS BezierCurve {
      
    
    public:
      
      /*------------------------------------------------------------------------*/
      /** \brief Constructor of a quadratic bezier curve from 3 control points.
       * 
       * \param AP1 first control point
       * \param AP2 second control point
       * \param AP3 third control point
       */
      BezierCurve(const Point& AP1, const Point& AP2, const Point& AP3);

      /*------------------------------------------------------------------------*/
      /** \brief Constructor of a cubic bezier curve from 4 control points.
       * 
       * \param AP1 first control point
       * \param AP2 second control point
       * \param AP3 third control point
       * \param AP4 fourth control point
       */
      BezierCurve(const Point& AP1, const Point& AP2, const Point& AP3,
		  const Point& AP4);


      /*------------------------------------------------------------------------*/
      /** \brief Constructor of a cubic bezier curve from 2 end points and two
       *         tangent vectors
       * 
       * \param AP1 first end point
       * \param AV1 first derivative vector at AP1
       * \param AP2 second end point
       * \param AV2 second derivative vector at AP1
       */
      BezierCurve(const Point& AP1, const Vector3d& AV1,
		  const Point& AP2, const Vector3d& AV2);


      /*------------------------------------------------------------------------*/
      /** \brief Constructor of a bezier curve from an ordered set of control 
       *         points.
       * 
       * \param APts the set of control points to define the curve
       */
      BezierCurve(const std::vector<Point>& APts);

      /*------------------------------------------------------------------------*/
      /** \brief Returns the point located on the parametric curve wit
       * 
       * \param AT the parameter in [0..1]
       *
       * \return The point located at (*this)(AT)
       */
      Point operator()(const double& AT) const;

      /*------------------------------------------------------------------------*/
      /** \brief Returns a set of point that discretize the curve in ANb segment
       *         in the parametric space
       * 
       * \param ANb the number of segment we want to have
       *
       * \return The set of point discretizing (*this);
       */
      std::vector<Point> getDiscretization(const int ANb) const;
      
    private:
      /** ordered list of control points of the curve */
      std::vector<Point> m_control_points;
      
    };
    /*----------------------------------------------------------------------------*/
  }
  /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_BEZIERCURVE_H_ */
/*----------------------------------------------------------------------------*/
