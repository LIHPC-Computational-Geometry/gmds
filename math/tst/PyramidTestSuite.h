/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/math/Pyramid.h>
#include <gmds/math/Point.h>
#include <gmds/math/Tetrahedron.h>
#include <gmds/math/Triangle.h>
/*----------------------------------------------------------------------------*/
TEST(Pyramid, testPyramidConstructor)
{
    // Test default constructor
    gmds::math::Pyramid pyramid1;
    ASSERT_EQ(gmds::math::Point(0, 0, 0), pyramid1.getPoint(0));

    // Test constructor with specified points
    gmds::math::Point points[5] = {
        gmds::math::Point(1, 2, 3),
        gmds::math::Point(4, 5, 6),
        gmds::math::Point(7, 8, 9),
        gmds::math::Point(10, 11, 12),
        gmds::math::Point(13, 14, 15)
    };

    gmds::math::Pyramid pyramid2(points);
    ASSERT_EQ(points[0], pyramid2.getPoint(0));
    ASSERT_EQ(points[1], pyramid2.getPoint(1));
    ASSERT_EQ(points[2], pyramid2.getPoint(2));
    ASSERT_EQ(points[3], pyramid2.getPoint(3));
    ASSERT_EQ(points[4], pyramid2.getPoint(4));
}

TEST(Pyramid, testPyramidCenter)
{
    gmds::math::Point points[5] = {
        gmds::math::Point(1, 2, 3),
        gmds::math::Point(4, 5, 6),
        gmds::math::Point(7, 8, 9),
        gmds::math::Point(10, 11, 12),
        gmds::math::Point(13, 14, 15)
    };

    gmds::math::Pyramid pyramid(points);

    // Test center calculation
    gmds::math::Point center = pyramid.getCenter();
    ASSERT_EQ(gmds::math::Point(7, 8, 9), center);
}

TEST(Pyramid, testPyramidVolume)
{
    gmds::math::Point points[5] = {
        gmds::math::Point(0, 0, 0),
        gmds::math::Point(1, 0, 0),
        gmds::math::Point(1, 1, 0),
        gmds::math::Point(0, 1, 0),
        gmds::math::Point(0, 0, 1)
    };

    gmds::math::Pyramid pyramid(points);

    // Test volume calculation
    double volume = pyramid.getVolume();
    ASSERT_DOUBLE_EQ(1. / 6., volume);
}

TEST(Pyramid, testPyramidScaledJacobian)
{
    gmds::math::Point points[5] = {
        gmds::math::Point(0, 0, 0),
        gmds::math::Point(1, 0, 0),
        gmds::math::Point(1, 1, 0),
        gmds::math::Point(0, 1, 0),
        gmds::math::Point(0, 0, 1)
    };

    gmds::math::Pyramid pyramid(points);

    // Test scaled Jacobian calculation
    double scaledJacobian = pyramid.computeScaledJacobian();
    ASSERT_DOUBLE_EQ(1. / 3., scaledJacobian);
}

TEST(Pyramid, testPyramidNormalizedScaledJacobian)
{
    gmds::math::Point points[5] = {
        gmds::math::Point(0, 0, 0),
        gmds::math::Point(1, 0, 0),
        gmds::math::Point(1, 1, 0),
        gmds::math::Point(0, 1, 0),
        gmds::math::Point(0, 0, 1)
    };

    gmds::math::Pyramid pyramid(points);

    // Test normalized scaled Jacobian calculation
    double normalizedScaledJacobian = pyramid.computeNormalizedScaledJacobian();
    ASSERT_DOUBLE_EQ(sqrt(2. / 3.), normalizedScaledJacobian);
}