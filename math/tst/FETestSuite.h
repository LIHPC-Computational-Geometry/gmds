#include <gtest/gtest.h>
#include <gmds/math/Segment.h>
#include <gmds/math/FE.h>

using namespace gmds::math;


TEST(TriangleP1, BMatrixCalculation)
{
    // Given three points representing a triangle
    Point P1({0.0, 0.0});
    Point P2({1.0, 0.0});
    Point P3({0.0, 1.0});

    // When calculating the B matrix
    auto B_matrix = TriangleP1::B(P1, P2, P3);

    // Then the B matrix should be correctly calculated
    ASSERT_EQ(B_matrix.rows(), 2);
    ASSERT_EQ(B_matrix.columns(), 2);
}

TEST(TriangleP1, StiffnessMatrixCalculation)
{
    // Given three points representing a triangle
    Point P1({0.0, 0.0});
    Point P2({1.0, 0.0});
    Point P3({0.0, 1.0});

    // When calculating the stiffness matrix
    auto stiffness_matrix = TriangleP1::stiffnessMatrix(P1, P2, P3);

    // Then the stiffness matrix should be correctly calculated
    ASSERT_EQ(stiffness_matrix.rows(), 3);
    ASSERT_EQ(stiffness_matrix.columns(), 3);
}

TEST(TetrahedronP1, StiffnessMatrixCalculation)
{
    // Given four points representing a tetrahedron
    Point P1({0.0, 0.0});
    Point P2({1.0, 0.0});
    Point P3({0.0, 1.0});
    Point P4({0.0, 0.0, 1.0});

    // When calculating the stiffness matrix
    auto stiffness_matrix = TetrahedronP1::stiffnessMatrix(P1, P2, P3, P4);

    // Then the stiffness matrix should be correctly calculated
    ASSERT_EQ(stiffness_matrix.rows(), 4);
    ASSERT_EQ(stiffness_matrix.columns(), 4);
}