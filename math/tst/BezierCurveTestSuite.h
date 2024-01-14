#include <gtest/gtest.h>
#include <gmds/math/BezierCurve.h>

// Test the constructor with three points
TEST(BezierCurveTest, ConstructorWithThreePoints) {
    gmds::math::Point p1(0, 0, 0);
    gmds::math::Point p2(1, 1, 1);
    gmds::math::Point p3(2, 2, 2);
    
    gmds::math::BezierCurve curve(p1, p2, p3);
    
    // Check if control points are correctly initialized
    ASSERT_EQ(curve.getControlPoints().size(), 3);
    ASSERT_EQ(curve.getControlPoints()[0], p1);
    ASSERT_EQ(curve.getControlPoints()[1], p2);
    ASSERT_EQ(curve.getControlPoints()[2], p3);
}

// Test the constructor with four points
TEST(BezierCurveTest, ConstructorWithFourPoints) {
    gmds::math::Point p1(0, 0, 0);
    gmds::math::Point p2(1, 1, 1);
    gmds::math::Point p3(2, 2, 2);
    gmds::math::Point p4(3, 3, 3);
    
    gmds::math::BezierCurve curve(p1, p2, p3, p4);
    
    // Check if control points are correctly initialized
    ASSERT_EQ(curve.getControlPoints().size(), 4);
    ASSERT_EQ(curve.getControlPoints()[0], p1);
    ASSERT_EQ(curve.getControlPoints()[1], p2);
    ASSERT_EQ(curve.getControlPoints()[2], p3);
    ASSERT_EQ(curve.getControlPoints()[3], p4);
}

// Test the operator () for different t parameters
TEST(BezierCurveTest, OperatorFunction) {
    gmds::math::Point p1(0, 0, 0);
    gmds::math::Point p2(1, 1, 1);
    gmds::math::Point p3(2, 2, 2);
    
    gmds::math::BezierCurve curve(p1, p2, p3);
    
    // Check if the return value is correct for t = 0, t = 0.5, t = 1
    ASSERT_EQ(curve(0), p1);
    ASSERT_NEAR(curve(0.5), gmds::math::Point(1, 1, 1), 1e-6);
    ASSERT_EQ(curve(1), p3);
}

// Test the getDiscretization function for a given number of samples
TEST(BezierCurveTest, GetDiscretization) {
    gmds::math::Point p1(0, 0, 0);
    gmds::math::Point p2(1, 1, 1);
    gmds::math::Point p3(2, 2, 2);
    
    gmds::math::BezierCurve curve(p1, p2, p3);
    
    // Get the discretization with a certain number of samples and check its size
    auto discretization = curve.getDiscretization(10);
    ASSERT_EQ(discretization.size(), 11);  // Including 0
}
