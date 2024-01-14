/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/math/Prism3.h>
#include <gmds/math/Point.h>
#include <gmds/math/Tetrahedron.h>
#include <gmds/math/Pyramid.h>
#include <gmds/math/Triangle.h>
/*----------------------------------------------------------------------------*/
TEST(Prism3, testPrism3Constructor)
{
    // Test default constructor
    gmds::math::Prism3 prism1;
    ASSERT_EQ(gmds::math::Point(0, 0, 0), prism1.getPoint(0));

    // Test constructor with specified points
    gmds::math::Point points[6] = {
        gmds::math::Point(1, 2, 3),
        gmds::math::Point(4, 5, 6),
        gmds::math::Point(7, 8, 9),
        gmds::math::Point(10, 11, 12),
        gmds::math::Point(13, 14, 15),
        gmds::math::Point(16, 17, 18)
    };

    gmds::math::Prism3 prism2(points);
    ASSERT_EQ(points[0], prism2.getPoint(0));
    ASSERT_EQ(points[1], prism2.getPoint(1));
    ASSERT_EQ(points[2], prism2.getPoint(2));
    ASSERT_EQ(points[3], prism2.getPoint(3));
    ASSERT_EQ(points[4], prism2.getPoint(4));
    ASSERT_EQ(points[5], prism2.getPoint(5));
}

TEST(Prism3, testPrism3Center)
{
    gmds::math::Point points[6] = {
        gmds::math::Point(1, 2, 3),
        gmds::math::Point(4, 5, 6),
        gmds::math::Point(7, 8, 9),
        gmds::math::Point(10, 11, 12),
        gmds::math::Point(13, 14, 15),
        gmds::math::Point(16, 17, 18)
    };

    gmds::math::Prism3 prism(points);

    // Test center calculation
    gmds::math::Point center = prism.getCenter();
    ASSERT_EQ(gmds::math::Point(8.5, 9.5, 10.5), center);
}

TEST(Prism3, testPrism3Volume)
{
    gmds::math::Point points[6] = {
        gmds::math::Point(0, 0, 0),
        gmds::math::Point(1, 0, 0),
        gmds::math::Point(0, 1, 0),
        gmds::math::Point(0, 0, 1),
        gmds::math::Point(1, 0, 1),
        gmds::math::Point(0, 1, 1)
    };

    gmds::math::Prism3 prism(points);

    // Test volume calculation
    double volume = prism.getVolume();
    ASSERT_DOUBLE_EQ(0.5, volume);
}

TEST(Prism3, testPrism3ScaledJacobian)
{
    gmds::math::Point points[6] = {
        gmds::math::Point(0, 0, 0),
        gmds::math::Point(1, 0, 0),
        gmds::math::Point(0, 1, 0),
        gmds::math::Point(0, 0, 1),
        gmds::math::Point(1, 0, 1),
        gmds::math::Point(0, 1, 1)
    };

    gmds::math::Prism3 prism(points);

    // Test scaled Jacobian calculation
    double scaledJacobian = prism.computeScaledJacobian();
    ASSERT_DOUBLE_EQ(1.0, scaledJacobian);
}

TEST(Prism3, testPrism3NormalizedScaledJacobian)
{
    gmds::math::Point points[6] = {
        gmds::math::Point(0, 0, 0),
        gmds::math::Point(1, 0, 0),
        gmds::math::Point(0, 1, 0),
        gmds::math::Point(0, 0, 1),
        gmds::math::Point(1, 0, 1),
        gmds::math::Point(0, 1, 1)
    };

    gmds::math::Prism3 prism(points);

    // Test normalized scaled Jacobian calculation
    double normalizedScaledJacobian = prism.computeNormalizedScaledJacobian();
    ASSERT_DOUBLE_EQ(sqrt(3) / 2.0, normalizedScaledJacobian);
}