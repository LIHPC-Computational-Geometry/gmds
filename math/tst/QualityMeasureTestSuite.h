/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/math/QualityMeasure.h>

using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(QualityMeasureClass, testSphereRation)
{
    math::Point p1(1,0,0);
    math::Point p2(0,1,0);
    math::Point p3(0,0,1);
    math::Point p4(0,-1,0);
    double q = math::QualityMeasure::sphereRatio(p1,p2,p3,p4);
    ASSERT_NEAR(0.833779,q, 0.0001);
}
