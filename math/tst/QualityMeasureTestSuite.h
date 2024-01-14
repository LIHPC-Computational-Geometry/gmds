/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/math/QualityMeasure.h>
/*----------------------------------------------------------------------------*/
using namespace gmds::math;

/*----------------------------------------------------------------------------*/
TEST(QualityMeasureClass, testSphereRatio)
{
    // Test case 1
    {
        Point AP0(0, 0, 0);
        Point AP1(1, 0, 0);
        Point AP2(0, 1, 0);
        Point AP3(0, 0, 1);

        double ratio = QualityMeasure::sphereRatio(AP0, AP1, AP2, AP3);

        ASSERT_NEAR(ratio, 0.942809, pow(10,-6));
    }

    // Test case 2
    {
        Point AP0(0, 0, 4);
        Point AP1(0, 0, 0);
        Point AP2(1, 0, 4);
        Point AP3(1, 0, 0);

        double ratio = QualityMeasure::sphereRatio(AP0, AP1, AP2, AP3);

        ASSERT_NEAR(ratio, 1.0, pow(10,-6));
    }
}

/*----------------------------------------------------------------------------*/
TEST(QualityMeasureClass, testExtremAngles)
{
    // Test case 1
    {
        Triangle AT(Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0));

        double smallAngle, largeAngle;
        QualityMeasure::extremAngles(AT, smallAngle, largeAngle);

        ASSERT_NEAR(smallAngle, 0.25, pow(10,-6));
        ASSERT_NEAR(largeAngle, 0.75, pow(10,-6));
    }

    // Test case 2
    {
        Point AP0(0, 0, 0);
        Point AP1(1, 0, 0);
        Point AP2(0, 1, 0);

        double smallAngle, largeAngle;
        QualityMeasure::extremAngles(AP0, AP1, AP2, smallAngle, largeAngle);

        ASSERT_NEAR(smallAngle, 0.25, pow(10,-6));
        ASSERT_NEAR(largeAngle, 0.75, pow(10,-6));
    }
}
/*----------------------------------------------------------------------------*/
