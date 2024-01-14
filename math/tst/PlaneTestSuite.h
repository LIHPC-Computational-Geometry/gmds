/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/math/Plane.h>
#include <gmds/math/Point.h>
#include <gmds/math/Segment.h>
/*----------------------------------------------------------------------------*/
TEST(Plane, testPlaneConstructors)
{
    // Test default constructor
    gmds::math::Plane plane1;
    ASSERT_EQ(gmds::math::Point(), plane1.getPoint());
    ASSERT_EQ(gmds::math::Vector3d(), plane1.getNormal());

    // Test constructor with a point and normal vector
    gmds::math::Point point2(1.0, 2.0, 3.0);
    gmds::math::Vector3d normal2(4.0, 5.0, 6.0);
    gmds::math::Plane plane2(point2, normal2);
    ASSERT_EQ(point2, plane2.getPoint());
    ASSERT_EQ(normal2.normalized(), plane2.getNormal());

    // Test constructor with three points
    gmds::math::Point point3(1.0, 2.0, 3.0);
    gmds::math::Point point4(4.0, 5.0, 6.0);
    gmds::math::Point point5(7.0, 8.0, 9.0);
    gmds::math::Plane plane3(point3, point4, point5);
    gmds::math::Vector3d expectedNormal3 = (point4 - point3).cross(point5 - point3).normalized();
    ASSERT_EQ(point3, plane3.getPoint());
    ASSERT_EQ(expectedNormal3, plane3.getNormal());

    // Test constructor with a Triangle
    gmds::math::Triangle triangle(point3, point4, point5);
    gmds::math::Plane plane4(triangle);
    ASSERT_EQ(point3, plane4.getPoint());
    ASSERT_EQ(expectedNormal3, plane4.getNormal());
}

TEST(Plane, testEqualityOperators)
{
    gmds::math::Point point1(1.0, 2.0, 3.0);
    gmds::math::Vector3d normal1(4.0, 5.0, 6.0);
    gmds::math::Plane plane1(point1, normal1);

    gmds::math::Point point2(1.0, 2.0, 3.0);
    gmds::math::Vector3d normal2(4.0, 5.0, 6.0);
    gmds::math::Plane plane2(point2, normal2);

    gmds::math::Point point3(7.0, 8.0, 9.0);
    gmds::math::Vector3d normal3(10.0, 11.0, 12.0);
    gmds::math::Plane plane3(point3, normal3);

    ASSERT_EQ(plane1, plane2);
    ASSERT_NE(plane1, plane3);
}

TEST(Plane, testIsIn)
{
    gmds::math::Point point(1.0, 2.0, 3.0);
    gmds::math::Vector3d normal(4.0, 5.0, 6.0);
    gmds::math::Plane plane(point, normal);

    gmds::math::Point pointInPlane(5.0, 4.0, 3.0);
    gmds::math::Point pointNotInPlane(7.0, 8.0, 9.0);

    ASSERT_TRUE(plane.isIn(pointInPlane));
    ASSERT_FALSE(plane.isIn(pointNotInPlane));
}

TEST(Plane, testIsStrictlyOnLeft)
{
    gmds::math::Point point(1.0, 2.0, 3.0);
    gmds::math::Vector3d normal(4.0, 5.0, 6.0);
    gmds::math::Plane plane(point, normal);

    gmds::math::Point pointLeft(2.0, 3.0, 4.0);
    gmds::math::Point pointRight(0.0, 1.0, 2.0);

    ASSERT_TRUE(plane.isStrictlyOnLeft(pointLeft));
    ASSERT_FALSE(plane.isStrictlyOnLeft(pointRight));
}

TEST(Plane, testDistance)
{
    gmds::math::Point point(1.0, 2.0, 3.0);
    gmds::math::Vector3d normal(4.0, 5.0, 6.0);
    gmds::math::Plane plane(point, normal);

    gmds::math::Point pointOnPlane(5.0, 4.0, 3.0);
    gmds::math::Point pointNotOnPlane(7.0, 8.0, 9.0);

    ASSERT_DOUBLE_EQ(0.0, plane.distance(pointOnPlane));
    ASSERT_DOUBLE_EQ(plane.distance(pointNotOnPlane), plane.distance(-pointNotOnPlane));
}

TEST(Plane, testProject)
{
    gmds::math::Point point(1.0, 2.0, 3.0);
    gmds::math::Vector3d normal(4.0, 5.0, 6.0);
    gmds::math::Plane plane(point, normal);

    gmds::math::Point pointOnPlane(5.0, 4.0, 3.0);
    gmds::math::Point pointNotOnPlane(7.0, 8.0, 9.0);

    ASSERT_EQ(pointOnPlane, plane.project(pointOnPlane));
    ASSERT_EQ(point, plane.project(pointNotOnPlane));
}

TEST(Plane, testIntersect)
{
    // Define a plane (1,0,0) passing through (0,0,0)
    gmds::math::Point planePoint(0.0, 0.0, 0.0);
    gmds::math::Vector3d planeNormal(1.0, 0.0, 0.0);
    gmds::math::Plane plane(planePoint, planeNormal);

    // Define a segment parallel to the plane
    gmds::math::Point segPoint1(0.0, 1.0, 2.0);
    gmds::math::Point segPoint2(0.0, 3.0, 4.0);
    gmds::math::Segment segmentParallel(segPoint1, segPoint2);

    // Define a segment intersecting the plane
    gmds::math::Point segPoint3(2.0, 1.0, 0.0);
    gmds::math::Point segPoint4(4.0, 3.0, 0.0);
    gmds::math::Segment segmentIntersecting(segPoint3, segPoint4);

    // Define a segment not intersecting the plane
    gmds::math::Point segPoint5(0.0, 1.0, 4.0);
    gmds::math::Point segPoint6(0.0, 3.0, 6.0);
    gmds::math::Segment segmentNotIntersecting(segPoint5, segPoint6);

    // Test intersection
    ASSERT_TRUE(plane.intersect(segmentParallel, true));
    ASSERT_TRUE(plane.intersect(segmentIntersecting, true));
    ASSERT_FALSE(plane.intersect(segmentNotIntersecting, true));
}

TEST(Plane, testIntersectDetails)
{
    // Define a plane (1,0,0) passing through (0,0,0)
    gmds::math::Point planePoint(0.0, 0.0, 0.0);
    gmds::math::Vector3d planeNormal(1.0, 0.0, 0.0);
    gmds::math::Plane plane(planePoint, planeNormal);

    // Define a segment parallel to the plane
    gmds::math::Point segPoint1(0.0, 1.0, 2.0);
    gmds::math::Point segPoint2(0.0, 3.0, 4.0);
    gmds::math::Segment segmentParallel(segPoint1, segPoint2);

    // Define a segment intersecting the plane
    gmds::math::Point segPoint3(2.0, 1.0, 0.0);
    gmds::math::Point segPoint4(4.0, 3.0, 0.0);
    gmds::math::Segment segmentIntersecting(segPoint3, segPoint4);

    // Define a segment not intersecting the plane
    gmds::math::Point segPoint5(0.0, 1.0, 4.0);
    gmds::math::Point segPoint6(0.0, 3.0, 6.0);
    gmds::math::Segment segmentNotIntersecting(segPoint5, segPoint6);

    // Variables to store intersection details
    gmds::math::Point intersectionPoint;
    double w0, w1;

    // Test intersection details
    ASSERT_EQ(gmds::math::Plane::NO_INTERSECTION, plane.intersect(segmentParallel, intersectionPoint, w0, w1, true));
    ASSERT_EQ(gmds::math::Plane::SEGMENT_MIDDLE, plane.intersect(segmentIntersecting, intersectionPoint, w0, w1, true));
    ASSERT_EQ(gmds::math::Plane::NO_INTERSECTION, plane.intersect(segmentNotIntersecting, intersectionPoint, w0, w1, true));
}
/*----------------------------------------------------------------------------*/