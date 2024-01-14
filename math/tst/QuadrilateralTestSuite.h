/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/math/QualityMeasure.h>
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
#include <gmds/math/Triangle.h>
/*----------------------------------------------------------------------------*/
using namespace gmds::math;
/*----------------------------------------------------------------------------*/
TEST(QualityMeasureClass, testSphereRatio)
{
    // Define points for a tetrahedron
    Point p0(0, 0, 0), p1(1, 0, 0), p2(0.5, 1, 0), p3(0.5, 0.5, 1);

    double ratio = QualityMeasure::sphereRatio(p0, p1, p2, p3);

    // Assuming some expected value based on known input
    ASSERT_NEAR(ratio, 0.5, pow(10, -6));
}

TEST(QualityMeasureClass, testExtremAngles)
{
    // Define points for a triangle
    Point p0(0, 0, 0), p1(1, 0, 0), p2(0.5, 1, 0);
    Triangle triangle(p0, p1, p2);

    double smallAngle, largeAngle;
    QualityMeasure::extremAngles(triangle, smallAngle, largeAngle);

    // Assuming some expected values based on known input
    ASSERT_NEAR(smallAngle, 0.0, pow(10, -6));
    ASSERT_NEAR(largeAngle, 1.0, pow(10, -6));
}

TEST(QualityMeasureClass, testExtremAnglesWithPoints)
{
    // Define points for a triangle
    Point p0(0, 0, 0), p1(1, 0, 0), p2(0.5, 1, 0);

    double smallAngle, largeAngle;
    QualityMeasure::extremAngles(p0, p1, p2, smallAngle, largeAngle);

    // Assuming some expected values based on known input
    ASSERT_NEAR(smallAngle, 0.0, pow(10, -6));
    ASSERT_NEAR(largeAngle, 1.0, pow(10, -6));
}

/*----------------------------------------------------------------------------*/
