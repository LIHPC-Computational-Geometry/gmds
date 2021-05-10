/*----------------------------------------------------------------------------*/
/*
 * FE.h
 *
 *  Created on: July 31, 2017
 *      Author: ledoux franck
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_FE_H
#define GMDS_FE_H
/*----------------------------------------------------------------------------*/
// gmds file headers
/*----------------------------------------------------------------------------*/
#include <gmds/math/Matrix.h>
#include <gmds/math/Vector.h>
#include <gmds/math/Point.h>
/*----------------------------------------------------------------------------*/
#include<cmath>
#include<string.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
    /*------------------------------------------------------------------------*/
    namespace math {

        /*--------------------------------------------------------------------*/
        /** \class FE
         *  \brief Describe the services provided by a FE Element built from
         *         a mesh cell using usual reference elements to perform
         *         computations. Points live in the physical
         *         elements while nodes live in the reference element.
         */
        class SegmentP1{
        public:
            static Matrix<2,2,double> stiffnessMatrix(const Point& AP1,
                                                      const Point& AP2) ;
        };
        class TriangleP1{
        public:
            static Matrix<3,3,double> stiffnessMatrix(const Point& AP1,
                                                      const Point& AP2,
                                                      const Point& AP3) ;
        private:
            static Vector2d grad_ref[3];
            static Matrix<2,2,double> B(const Point& AP1,
                                        const Point& AP2,
                                        const Point& AP3);
        };
        class TetrahedronP1{
        public:
            static Matrix<4,4,double> stiffnessMatrix(const Point& AP1,
                                                      const Point& AP2,
                                                      const Point& AP3,
                                                      const Point& AP4);
        private:
            static Vector3d grad_ref[4];

        };
        /*--------------------------------------------------------------------*/
    }//namespace math
/*----------------------------------------------------------------------------*/
}//namespace gmds
/*----------------------------------------------------------------------------*/
#endif //GMDS_FE_H
/*----------------------------------------------------------------------------*/
