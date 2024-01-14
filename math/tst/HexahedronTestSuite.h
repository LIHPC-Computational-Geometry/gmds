#include <gtest/gtest.h>
#include <gmds/math/Hexahedron.h>
#include <gmds/math/Point.h>
#include <gmds/math/Triangle.h>

using namespace gmds::math;

TEST(Hexahedron, Constructors)
{
    // Hexahedron created using different constructors

    Hexahedron hex1; // Default constructor
    Hexahedron hex2(Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0),
                    Point(0, 0, 1), Point(1, 0, 1), Point(1, 1, 1), Point(0, 1, 1)); // Constructor with individual points
    Hexahedron hex3({Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0),
                     Point(0, 0, 1), Point(1, 0, 1), Point(1, 1, 1), Point(0, 1, 1)}); // Constructor with an array of points
    Hexahedron hex4(std::vector<Point>{Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0),
                                       Point(0, 0, 1), Point(1, 0, 1), Point(1, 1, 1), Point(0, 1, 1)}); // Constructor with a vector of points
    Hexahedron hex5(hex1); // Copy constructor

    // Verify that the Hexahedron objects are correctly initialized

    ASSERT_TRUE(hex1.isValid());
    ASSERT_TRUE(hex2.isValid());
    ASSERT_TRUE(hex3.isValid());
    ASSERT_TRUE(hex4.isValid());
    ASSERT_TRUE(hex5.isValid());
}

TEST(Hexahedron, Getters)
{
    // Given a Hexahedron

    Hexahedron hex(Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0),
                   Point(0, 0, 1), Point(1, 0, 1), Point(1, 1, 1), Point(0, 1, 1));

    // When retrieving the center, volume, and other properties

    Point center = hex.getCenter();
    double volume = hex.getVolume();
    bool isValid = hex.isValid();
    double meanEdgeLength = hex.computeMeanEdgeLength();

    // Verify the correctness of the results

    ASSERT_EQ(center, Point(0.5, 0.5, 0.5));
    ASSERT_GT(volume, 0);
    ASSERT_TRUE(isValid);
    ASSERT_GT(meanEdgeLength, 0);
}

TEST(Hexahedron, Intersect)
{
    // Given a Hexahedron and a Triangle

    Hexahedron hex(Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0),
                   Point(0, 0, 1), Point(1, 0, 1), Point(1, 1, 1), Point(0, 1, 1));

    Triangle tri(Point(0.5, 0.5, -0.5), Point(1.5, 0.5, -0.5), Point(1.5, 1.5, -0.5));

    // When checking for intersection

    bool isIntersecting = hex.intersect(tri);

    // Then, verify the correctness of the result

    ASSERT_TRUE(isIntersecting);
}

TEST(Hexahedron, Jacobian)
{
    // Given a Hexahedron

    Hexahedron hex(Point(0, 0, 0), Point(1, 0, 0), Point(1, 1, 0), Point(0, 1, 0),
                   Point(0, 0, 1), Point(1, 0, 1), Point(1, 1, 1), Point(0, 1, 1));

    // When calculating the Jacobian matrix for each vertex

    Matrix<3, 3, double> jacobian0 = hex.jacobian(0);
    Matrix<3, 3, double> jacobian1 = hex.jacobian(1);
    Matrix<3, 3, double> jacobian2 = hex.jacobian(2);
    Matrix<3, 3, double> jacobian3 = hex.jacobian(3);
    Matrix<3, 3, double> jacobian4 = hex.jacobian(4);
    Matrix<3, 3, double> jacobian5 = hex.jacobian(5);
    Matrix<3, 3, double> jacobian6 = hex.jacobian(6);
    Matrix<3, 3, double> jacobian7 = hex.jacobian(7);

}