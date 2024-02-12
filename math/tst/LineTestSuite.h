#include <gtest/gtest.h>
#include "YourGeometryLibrary.h" 

// Test case for Line-Triangle intersection
TEST(LineTriangleIntersection, IntersectionCheck) {
    // Create a Line and a Triangle for testing
    Line line(Point(0, 0, 0), Point(1, 1, 1));
    Triangle triangle(Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0));

    // Check intersection
    bool result = line.intersect(triangle, someParameterValue);
    ASSERT_TRUE(result);
}

// Test case for Line projection
TEST(LineProjection, ProjectPoint) {
    // Create a Line and a Point for testing
    Line line(Point(0, 0, 0), Point(1, 1, 1));
    Point pointToProject(2, 2, 2);

    // Project the point onto the line
    Point projectedPoint = line.project(pointToProject);

    // Check if the projected point is correct
    ASSERT_EQ(projectedPoint, Point(1, 1, 1));
}

// Test case for Line-Segment 2D distance computation
TEST(LineSegmentDistance, Distance2DCheck) {
    // Create a Line and a Segment for testing
    Line line(Point(0, 0, 0), Point(1, 1, 1));
    Segment segment(Point(0, 0, 0), Point(2, 2, 2));

    // Variables to store the closest points
    Point closestPoint1, closestPoint2;

    // Compute the 2D distance between the Line and the Segment
    TCoord distance = line.distance2D(segment, closestPoint1, closestPoint2);

    // Check if the computed distance is correct
    ASSERT_EQ(distance, sqrt(3));

    // Check if the closest points are correct
    ASSERT_EQ(closestPoint1, Point(1, 1, 1));
    ASSERT_EQ(closestPoint2, Point(1, 1, 1));
}
