/*----------------------------------------------------------------------------*/
#ifndef GMDS_ORIENTATION_TESTSUITE_H
#define GMDS_ORIENTATION_TESTSUITE_H
/*----------------------------------------------------------------------------*/
#include "gtest/gtest.h"
#include <gmds/math/Orientation.h>
#include <gmds/math/Tetrahedron.h>
/*----------------------------------------------------------------------------*/
TEST(OrientationClass, orient3d_basic)
{
    gmds::math::Orientation::initialize();

    gmds::math::Point p0(0,0,0);
    gmds::math::Point p1(1,0,0);
    gmds::math::Point p2(0,1,0);

    gmds::math::Point guess(0,0,1);
    gmds::math::Point guess2(0.25,0.5,0);
    gmds::math::Tetrahedron t(guess,p0,p1,p2);
    std::cout<<"volume: "<<t.getVolume()<<std::endl;
    ASSERT_EQ(gmds::math::Orientation::NEGATIVE,
              gmds::math::Orientation::orient3d(guess,p0,p1,p2));
    ASSERT_EQ(gmds::math::Orientation::POSITIVE,
              gmds::math::Orientation::orient3d(guess,p0,p2,p1));
    ASSERT_EQ(gmds::math::Orientation::ZERO,
              gmds::math::Orientation::orient3d(p0,p0,p1,p2));
    ASSERT_EQ(gmds::math::Orientation::ZERO,
              gmds::math::Orientation::orient3d(p1,p0,p1,p2));
    ASSERT_EQ(gmds::math::Orientation::ZERO,
              gmds::math::Orientation::orient3d(p2,p0,p1,p2));
    ASSERT_EQ(gmds::math::Orientation::ZERO,
              gmds::math::Orientation::orient3d(guess2,p0,p1,p2));

}
/*----------------------------------------------------------------------------*/
#endif //GMDS_ORIENTATION_TESTSUITE_H
/*----------------------------------------------------------------------------*/
