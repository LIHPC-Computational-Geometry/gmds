/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/math/Chart.h>
#include <gmds/math/AxisAngleRotation.h>

using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(ChartClass, testBasic)
{
    math::Vector3d v1(0,1,1);
    math::Vector3d v2(1,0,1);
    math::Chart c(v1,v2);

    ASSERT_EQ(c.get(0), v1);
    ASSERT_EQ(c.get(1), v2);
    ASSERT_EQ(c.get(2), v1.cross(v2));

    ASSERT_ANY_THROW(c.get(4););
    ASSERT_ANY_THROW(c.get(-1););

    math::Chart default_ch;
}
/*----------------------------------------------------------------------------*/
TEST(ChartClass, testRotationMatrix)
{
    math::Chart ref; // 1,0,0 and 0,1,0 and 0,0,1

    math::Matrix<3,3,double> M = ref.computeRotationTo(ref);
    math::Matrix<3,3,double> I = math::Matrix<3,3,double>::identity();

    for(auto i=0;i<3;i++)
        for(auto j=0;j<3;j++)
            ASSERT_NEAR(std::abs(M(i,j)),std::abs(I(i,j)),1e-12);

}
/*----------------------------------------------------------------------------*/
