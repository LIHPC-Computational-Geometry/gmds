/*----------------------------------------------------------------------------*/
/*
 * QualityMeasure.h
 *
 *  Created on: 12 april 2017
 *      Author: ledoux f.
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_QUALITY_MEASURE_H_
#define GMDS_MATH_QUALITY_MEASURE_H_
/*----------------------------------------------------------------------------*/
// gmds file headers
/*----------------------------------------------------------------------------*/
#include <gmds/math/Point.h>
#include <gmds/math/Triangle.h>
#include <gmds/utils/CommonTypes.h>
#include "GMDSMath_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace math {
/*----------------------------------------------------------------------------*/
/** \class QualityMeasure
 *  \brief Defines a QualityMeasure
 **/
class GMDSMath_API QualityMeasure
{
public:
        /*------------------------------------------------------------------------*/
        /** \brief  Compute a shape metric for a tetrahedron T=(\p AP0, \p AP1,
         *          \p AP2, \p AP3) as 3*r/R with r the radius of the inscribed
         *          sphere of T and R the radius of the circumscribed sphere of T
         *
         * \param[in] AP0 first tet point
         * \param[in] AP1 second tet point
         * \param[in] AP2 third tet point
         * \param[in] AP3 fourth tet point
         *
         * \return  a value whose max 1 means a regular tetrahedron and a negative
         *          value means an inverted element
         */
        static double sphereRatio(const Point& AP0, const Point& AP1, const Point& AP2, const Point& AP3);

    /*------------------------------------------------------------------------*/
    /** \brief  Returns the  smallest and largest angle inside a 3D triangle.
     *          Angles are returned in degrees.
     *
     * \param[in]  AP0         first triangle point
     * \param[in]  AP1         second triangle point
     * \param[in]  AP2         third triangle point
     * \param[out] ASmallAngle smallest angle
     * \param[out] ALargeAngle largest angle
     */
    static void extremAngles(const Point& AP0, const Point& AP1, const Point& AP2,
                             double& ASmallAngle, double& ALargeAngle);
    /*------------------------------------------------------------------------*/
    /** \brief  Returns the  smallest and largest angle inside a 3D triangle.
     *          Angles are returned in degrees.
     *
     * \param[in]  AT          a triangle
     * \param[out] ASmallAngle smallest angle
     * \param[out] ALargeAngle largest angle
     */
    static void extremAngles(const Triangle& AT,
                             double& ASmallAngle, double& ALargeAngle);
};
}
/*----------------------------------------------------------------------------*/
}  // namespace math
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_QUALITY_MEASURE_H_ */
/*----------------------------------------------------------------------------*/
