/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/math/Ray.h>
#include <gmds/math/Segment.h>
#include <gmds/math/Plane.h>
#include <gmds/math/Triangle.h>
/*----------------------------------------------------------------------------*/
using namespace gmds::math;

/*----------------------------------------------------------------------------*/
TEST(RayClass, testGetters)
{
    // Test case 1
    {
        Point pnt(1, 2, 3);
        Vector3d dir(0.5, 0.5, 0.5);

        Ray ray(pnt, dir);

        ASSERT_EQ(ray.getPoint(), pnt);
        ASSERT_EQ(ray.getDir(), dir);

        // Additional checks for dirUnit
        ASSERT_EQ(ray.getDirUnit(), dir.getNormalize());
        ASSERT_EQ(ray.getDirUnit(), ray.getDirUnit()); // Check if the same instance is returned
    }

}

/*----------------------------------------------------------------------------*/
TEST(RayClass, testIntersect2D)
{
    // Test case 1
    {
        Ray ray(Point(0, 0, 0), Vector3d(1, 1, 0));
        Segment segment(Point(1, 0, 0), Point(0, 1, 0));

        Point intersectionPoint;
        double param;

        bool result = ray.intersect2D(segment, intersectionPoint, param);

        ASSERT_TRUE(result);
        ASSERT_NEAR(intersectionPoint.X(), 0.5, pow(10, -6));
        ASSERT_NEAR(intersectionPoint.Y(), 0.5, pow(10, -6));
        ASSERT_NEAR(param, 0.5, pow(10, -6));
    }
}

/*----------------------------------------------------------------------------*/
TEST(RayClass, testSecondMetIntersect2D)
{
    // Test case 1
    {
        Ray ray(Point(0, 0, 0), Vector3d(1, 1, 0));
        Segment segment(Point(1, 0, 0), Point(0, 1, 0));

        Point intersectionPoint;
        double paramRay, paramSegment;

        bool result = ray.SecondMetIntersect2D(segment, intersectionPoint, paramRay, paramSegment);

        ASSERT_TRUE(result);
        ASSERT_NEAR(intersectionPoint.X(), 0.5, pow(10, -6));
        ASSERT_NEAR(intersectionPoint.Y(), 0.5, pow(10, -6));
        ASSERT_NEAR(paramRay, 0.5, pow(10, -6));
        ASSERT_NEAR(paramSegment, 0.5, pow(10, -6));
    }
}

/*----------------------------------------------------------------------------*/
TEST(RayClass, testIntersect3D)
{
    // Test case 1
    {
        Ray ray(Point(0, 0, 0), Vector3d(1, 1, 1));
        Segment segment(Point(0, 0, 0), Point(1, 1, 1));

        Point intersectionPoint;
        double paramRay, paramSegment;

        bool result = ray.intersect3D(segment, intersectionPoint, paramRay, paramSegment);

        ASSERT_TRUE(result);
        ASSERT_NEAR(intersectionPoint.X(), 0.5, pow(10, -6));
        ASSERT_NEAR(intersectionPoint.Y(), 0.5, pow(10, -6));
        ASSERT_NEAR(intersectionPoint.Z(), 0.5, pow(10, -6));
        ASSERT_NEAR(paramRay, 0.5, pow(10, -6));
        ASSERT_NEAR(paramSegment, 0.5, pow(10, -6));
    }
}

/*----------------------------------------------------------------------------*/
TEST(RayClass, testIntersect2DWithRay)
{
    // Test case 1
    {
        Ray ray1(Point(0, 0, 0), Vector3d(1, 1, 0));
        Ray ray2(Point(1, 0, 0), Vector3d(0, 1, 0));

        Point intersectionPoint;

        bool result = ray1.intersect2D(ray2, intersectionPoint);

        ASSERT_TRUE(result);
        ASSERT_NEAR(intersectionPoint.X(), 0.5, pow(10, -6));
        ASSERT_NEAR(intersectionPoint.Y(), 0.5, pow(10, -6));
        ASSERT_NEAR(intersectionPoint.Z(), 0.0, pow(10, -6));
    }
}

/*----------------------------------------------------------------------------*/
TEST(RayClass, testIntersect3DWithPlane)
{
    // Test case 1
    {
        Ray ray(Point(0, 0, 0), Vector3d(1, 1, 1));
        Plane plane(Point(0, 0, 1), Vector3d(0, 0, 1));

        Point intersectionPoint;

        bool result = ray.intersect3D(plane, intersectionPoint);

        ASSERT_TRUE(result);
        ASSERT_NEAR(intersectionPoint.X(), 0.5, pow(10, -6));
        ASSERT_NEAR(intersectionPoint.Y(), 0.5, pow(10, -6));
        ASSERT_NEAR(intersectionPoint.Z(), 0.5, pow(10, -6));
    }
}

/*----------------------------------------------------------------------------*/
TEST(RayClass, testIntersect3DWithTriangle)
{
    // Test case 1
    {
        Ray ray(Point(0, 0, 0), Vector3d(1, 1, 1));
        Triangle triangle(Point(0, 0, 1), Point(1, 0, 0), Point(0, 1, 0));

        Point intersectionPoint;

        bool result = ray.intersect3D(triangle, intersectionPoint);

        ASSERT_TRUE(result);
        ASSERT_NEAR(intersectionPoint.X(), 0.5, pow(10, -6));
        ASSERT_NEAR(intersectionPoint.Y(), 0.5, pow(10, -6));
        ASSERT_NEAR(intersectionPoint.Z(), 0.5, pow(10, -6));
    }
}

/*----------------------------------------------------------------------------*/
TEST(RayClass, testProject)
{
    // Test case 1
    {
        Ray ray(Point(0, 0, 0), Vector3d(1, 1, 1));
        Point pointToProject(2, 2, 2);

        Point projectedPoint = ray.project(pointToProject);

        ASSERT_NEAR(projectedPoint.X(), 1.0, pow(10, -6));
        ASSERT_NEAR(projectedPoint.Y(), 1.0, pow(10, -6));
        ASSERT_NEAR(projectedPoint.Z(), 1.0, pow(10, -6));
    }
}
/*----------------------------------------------------------------------------*/
