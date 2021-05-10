/*----------------------------------------------------------------------------*/
/*
 * BezierTriangle.h
 *
 * Created on: 04/27/2017
 * Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDSSUITE_BEZIERTRIANGLE_H
#define GMDSSUITE_BEZIERTRIANGLE_H
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
        class EXPORT_GMDS BezierTriangle {


        public:

            /*------------------------------------------------------------------------*/
            /** \brief Constructor of a cubic bezier triangle from end points and
             *         their normals
             *
             * \param[in] AP1 first end point
             * \param[in] AP2 second end point
             * \param[in] AP3 third end point
             * \param[in] AN1 first normal at \p AP1
             * \param[in] AN2 first normal at \p AP2
             * \param[in] AN3 first normal at \p AP3
             */
            BezierTriangle(const Point& AP1, const Point& AP2, const Point& AP3,
                           const Vector3d& AN1, const Vector3d& AN2, const Vector3d& AN3);


            /*------------------------------------------------------------------------*/
            /** \brief Returns the point located on the parametric triangle (u,v) in
             *         (u,v,1-u-v)
             *
             * \param[in] AU the parameter u in [0..1]
             * \param[in] AV the parameter v in [0..1]
             *
             * \return The point located at (*this)(\p AU, \p AV, 1-\p AU -\p AV)
             */
            Point operator()(const double AU, const double AV) const;
            /*------------------------------------------------------------------------*/
            /** \brief Returns the normal located on the parametric triangle (u,v) in
             *         (u,v,1-u-v)
             *
             * \param[in] AU the parameter u in [0..1]
             * \param[in] AV the parameter v in [0..1]
             *
             * \return The normal to the triangle at (*this)(\p AU, \p AV, 1-\p AU -\p AV)
             */
            Vector3d normal(const double AU, const double AV) const;

            /*------------------------------------------------------------------------*/
            /** \brief Returns the point, normal, u and v derivatives of the parametric
             *         triangle (u,v) in (u,v,1-u-v)
             *
             * \param[in] AU the parameter u in [0..1]
             * \param[in] AV the parameter v in [0..1]
             *
             * \param[out] AP The point of the triangle at (*this)(\p AU, \p AV, 1-\p AU -\p AV)
             * \param[out] AN The normal to the triangle at (*this)(\p AU, \p AV, 1-\p AU -\p AV)
             * \param[out] ADU The u-derivative to the triangle at (*this)(\p AU, \p AV, 1-\p AU -\p AV)
             * \param[out] ADV The v-derivative to the triangle at (*this)(\p AU, \p AV, 1-\p AU -\p AV)
             */
            void geomInfo(const double AU, const double AV, Point& AP,
                          Vector3d& AN, Vector3d& ADU, Vector3d& ADV) const;

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
            Point m_B[10];
        };
        /*----------------------------------------------------------------------------*/
    }
    /*----------------------------------------------------------------------------*/
}
#endif //GMDSSUITE_BEZIERTRIANGLE_H
