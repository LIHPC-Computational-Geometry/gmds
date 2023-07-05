/*----------------------------------------------------------------------------*/
/*
 * Numerics.h
 *
 *  Created on:March 01, 2015
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_NUMERICS_H_
#define GMDS_MATH_NUMERICS_H_
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Constants.h>
#include <gmds/math/Point.h>
#include <gmds/math/Matrix.h>
#include "GMDSMath_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
  /*----------------------------------------------------------------------------*/
  namespace math{
      /*----------------------------------------------------------------------------*/
      bool GMDSMath_API isZero(TCoord a, TCoord AEpsilon = Constants::EPSILON);
	 
	  bool isZero2ndMethod(TCoord a);

      bool GMDSMath_API areEquals(TCoord a, TCoord b, TCoord AEpsilon = Constants::EPSILON);

    /*------------------------------------------------------------------------*/
    /** \brief Compute AVal modulo AMod
     */
    TCoord GMDSMath_API modulo(TCoord AVal, TCoord AMod);
    /*------------------------------------------------------------------------*/
    /** \brief Compute AVal modulo 2xPI
     */
    TCoord GMDSMath_API modulo2PI(TCoord AVal);
    /*------------------------------------------------------------------------*/
    /** \brief Returns the min value between AV1, AV2 and AV3
     */
    TCoord GMDSMath_API min3(TCoord AV1, TCoord AV2, TCoord AV3);
    /*------------------------------------------------------------------------*/
    /** \brief Returns the max value between AV1, AV2 and AV3
     */
    TCoord max3(TCoord AV1, TCoord AV2, TCoord AV3);
    /*------------------------------------------------------------------------*/
    /** \brief Returns the min value between AV1 and AV2
     */
    TCoord min2(TCoord AV1, TCoord AV2);
    /*------------------------------------------------------------------------*/
    /** \brief Returns the max value between AV1 and AV2
     */
    TCoord max2(TCoord AV1, TCoord AV2);
    /*------------------------------------------------------------------------*/
    /** \brief Returns the min value between AV1, AV2, AV3 and AV4
     */
    TCoord min4(TCoord AV1, TCoord AV2, TCoord AV3, TCoord AV4);
    /*------------------------------------------------------------------------*/
    /** \brief Returns the max value among AV1, AV2, AV3 and AV4
     */
    TCoord max4(TCoord AV1, TCoord AV2, TCoord AV3, TCoord AV4);
    /*------------------------------------------------------------------------*/
    /** \brief Returns the min value among eight values
     */
    TCoord min8(TCoord AV1, TCoord AV2, TCoord AV3, TCoord AV4, TCoord AV5, TCoord AV6, TCoord AV7, TCoord AV8);
    /*------------------------------------------------------------------------*/
    /** \brief Returns the max value among eight values
     */
    TCoord max8(TCoord AV1, TCoord AV2, TCoord AV3, TCoord AV4, TCoord AV5, TCoord AV6, TCoord AV7, TCoord AV8);
    /*------------------------------------------------------------------------*/
    /** \brief Returns true if AV1==AV2 [AEpsilon]
     */
    bool GMDSMath_API near(TCoord AV1, TCoord AV2, TCoord AEpsilon = Constants::EPSILON);

      /*------------------------------------------------------------------------*/
      /** \brief Returns the stiffness matrix for the triangle element defined by
       *         AP1, AP2 and AP3.
       *
       *    WARNING ONLY 2D!!!!! (x,y)
       *
       * \param AP1 a first point
       * \param AP2 a second point
       * \param AP3 a third point
       *
       * \return A 3x3 the stiffness matrix of (AP1,AP2,AP3)
       */
	  GMDSMath_API Matrix<3,3,double> stiffnessMatrix2D(const Point& AP1,
                                           const Point& AP2,
                                           const Point& AP3);
      
      /*------------------------------------------------------------------------*/
      /** \brief Compute an average plane going through the set of points \p AP
       *
       * \param[in]  AP           the set of points the plane must fit to
       * \param[out] APlanePnt    one point of the plane
       * \param[out] APlaneNormal normal vector to the plane
       *
       */
      void computeLeastSquarePlane(const std::vector<Point>& AP,
                                   Point& APlanePnt,
                                   math::Vector3d& APlaneNormal);
      
      
      /*------------------------------------------------------------------------*/
      /** \brief Solve a 2nd degree polynom ax2+bx+c=0
       *
       * \param[in]  AA coeff of X2
       * \param[in]  AB coeff of X
       * \param[in]  AC constant coeff
       * \param[out] AX solutions
       *
       * \return the number of real solutions (0, 1 or 2)
       */
	  GMDSMath_API int solve2ndDegreePolynomial( const double& AA,
                                   const double& AB,
                                   const double& AC,
                                   std::vector<double>& AX);
      
      /*------------------------------------------------------------------------*/
      /** \brief Considering tetrahedron T defined by points \p AP0, \p AP1,
       *         \p AP2 and \p AP3, this function returs the cotangent weight
       *         for edge[\p AP0, \p AP1] in T.
       *
       * \param[in]  AP0 a first tet point
       * \param[in]  AP1 a second tet point
       * \param[in]  AP2 a third tet point
       * \param[in]  AP3 a fourth tet point
       *
       * \return the cotangent weight of edge [\p AP0, \p AP1] in this tet.
       */
      
      double cotangentWeight(const Point& AP1,
                             const Point& AP2,
                             const Point& AP3,
                             const Point& AP4);
      
      /*------------------------------------------------------------------------*/
      /** \brief Considering tetrahedron T defined by points \p AP0, \p AP1,
       *         \p AP2 and \p AP3, this function returns the dihedral angle
       *         for edge[\p AP0, \p AP1] in T.
       *
       * \param[in]  AP0 a first tet point
       * \param[in]  AP1 a second tet point
       * \param[in]  AP2 a third tet point
       * \param[in]  AP3 a fourth tet point
       *
       * \return the dihedral angle of edge [\p AP0, \p AP1] in this tet.
       */
      
      double dihedralAngle(const Point& AP1,
                             const Point& AP2,
                             const Point& AP3,
                             const Point& AP4);
      
      

      /*------------------------------------------------------------------------*/
      /** \brief  Returns whether two axis-aligned bounding boxes intersect.
       *
       * \param AMinXYZ_0 the lower front left coordinates of the first aabbox
       * \param AMaxXYZ_0 the upper back right coordinates of the first aabbox
       * \param AMinXYZ_1 the upper back right coordinates of the second aabbox
       * \param AMaxXYZ_1 the upper back right coordinates of the second aabbox
       *
       * \return whether the two boxes intersect
       */
      bool intersectBoundingBox(const double AMinXYZ_0[3], const double AMaxXYZ_0[3], const double AMinXYZ_1[3], const double AMaxXYZ_1[3]);
      /*----------------------------------------------------------------------------*/
  } // namespace math
  /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_CONSTANTS_H_ */
/*----------------------------------------------------------------------------*/
