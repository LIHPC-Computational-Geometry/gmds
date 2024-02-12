#include <gtest/gtest.h>
#include <gmds/math/Numerics.h>

using namespace gmds::math;

// Test for the modulo function
TEST(NumericsTest, ModuloTest) {
    // Test with positive values
    ASSERT_DOUBLE_EQ(modulo(10.0, 3.0), 1.0);
    // Test with negative values
    ASSERT_DOUBLE_EQ(modulo(-10.0, 3.0), 2.0);
}

// Test for the modulo2PI function
TEST(NumericsTest, Modulo2PITest) {
    ASSERT_DOUBLE_EQ(modulo2PI(7.0), 1.0);
    ASSERT_DOUBLE_EQ(modulo2PI(-7.0), 5.0);
}

// Tests for the min and max functions
TEST(NumericsTest, MinMaxTest) {
    // Test with two values
    ASSERT_DOUBLE_EQ(min2(3.0, 5.0), 3.0);
    ASSERT_DOUBLE_EQ(max2(3.0, 5.0), 5.0);

    // Test with three values
    ASSERT_DOUBLE_EQ(min3(3.0, 5.0, 2.0), 2.0);
    ASSERT_DOUBLE_EQ(max3(3.0, 5.0, 2.0), 5.0);

    // Test with four values
    ASSERT_DOUBLE_EQ(min4(3.0, 5.0, 2.0, 8.0), 2.0);
    ASSERT_DOUBLE_EQ(max4(3.0, 5.0, 2.0, 8.0), 8.0);

    // Test with eight values
    ASSERT_DOUBLE_EQ(min8(3.0, 5.0, 2.0, 8.0, 1.0, 6.0, 7.0, 4.0), 1.0);
    ASSERT_DOUBLE_EQ(max8(3.0, 5.0, 2.0, 8.0, 1.0, 6.0, 7.0, 4.0), 8.0);
}

// Test for the isZero function
TEST(NumericsTest, IsZeroTest) {
    ASSERT_TRUE(isZero(0.0, 1e-6));
    ASSERT_TRUE(isZero(0.0, Constants::EPSILON));
    ASSERT_FALSE(isZero(0.00001, 1e-6));
}

// Test for the areEquals function
TEST(NumericsTest, AreEqualsTest) {
    ASSERT_TRUE(areEquals(1.0, 1.0, 1e-6));
    ASSERT_TRUE(areEquals(1.0, 1.00001, 1e-6));
    ASSERT_FALSE(areEquals(1.0, 2.0, 1e-6));
}

// Test for the dihedralAngle function
TEST(NumericsTest, DihedralAngleTest) {
    Point A(0, 0, 0);
    Point B(1, 0, 0);
    Point C(1, 1, 0);
    Point D(0, 1, 1);

    ASSERT_DOUBLE_EQ(dihedralAngle(A, B, C, D), Constants::PI / 4.0);
}

// Test for the cotangentWeight function
TEST(NumericsTest, CotangentWeightTest) {
    Point A(0, 0, 0);
    Point B(1, 0, 0);
    Point C(1, 1, 0);
    Point D(0, 1, 1);

    ASSERT_DOUBLE_EQ(cotangentWeight(A, B, C, D), 0.125);
}

// Test for the intersectBoundingBox function
TEST(NumericsTest, IntersectBoundingBoxTest) {
    double minXYZ0[3] = {0.0, 0.0, 0.0};
    double maxXYZ0[3] = {1.0, 1.0, 1.0};

    double minXYZ1[3] = {0.5, 0.5, 0.5};
    double maxXYZ1[3] = {1.5, 1.5, 1.5};

    double minXYZ2[3] = {2.0, 2.0, 2.0};
    double maxXYZ2[3] = {3.0, 3.0, 3.0};

    ASSERT_TRUE(intersectBoundingBox(minXYZ0, maxXYZ0, minXYZ1, maxXYZ1));
    ASSERT_FALSE(intersectBoundingBox(minXYZ0, maxXYZ0, minXYZ2, maxXYZ2));
}
