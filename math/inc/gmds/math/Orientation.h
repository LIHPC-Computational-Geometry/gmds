/*----------------------------------------------------------------------------*/
#ifndef GMDS_ORIENTATION_H
#define GMDS_ORIENTATION_H
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Constants.h>
#include <gmds/math/Point.h>
#include <gmds/math/Matrix.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
    namespace math {
        /**@brief Provide orientation methods. Some of them are based on the
         *        exact predicate designed by B. Levy, INRIA.
         *
         */
        class Orientation {
        public:
            enum Sign {
                NEGATIVE = -1,  /*!< Negative value */
                ZERO = 0,       /*!< Null  value */
                POSITIVE = 1    /*!< Positive value */
            };
            /**@brief Need to be called once before using Orientation methods
             */
            static void initialize();

            /**@brief Need to be called at the end
             */
            static void finalize();

            /**@brief check the volume orientation
             * @param AP0 First tetrahedral point
             * @param AP1 Second tetrahedral point
             * @param AP2 Third tetrahedral point
             * @param AP3 Fourth tetrahedral point
             * @return POSITIVE if the tet volume is positive, i.e. v03.dot(v01.cross(v02))>0
             * @return NEGATIVE if the tet volume is negative, i.e. v03.dot(v01.cross(v02))<0
-             * @return ZERO if the all points are coplanar
             */
            static Sign orient3d(const Point &AP0,
                                 const Point &AP1,
                                 const Point &AP2,
                                 const Point &AP3);


        };

/*----------------------------------------------------------------------------*/
    }
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif //GMDS_ORIENTATION_H
