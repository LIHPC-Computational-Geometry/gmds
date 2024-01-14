#include <gtest/gtest.h>
#include <gmds/math/AxisAngleRotation.h>

using namespace gmds::math;

// Test default constructor
TEST(AxisAngleRotationClass, DefaultConstructor)
{
    AxisAngleRotation axisAngleRotation;
    ASSERT_NO_THROW(axisAngleRotation.quaternion());
    ASSERT_NO_THROW(axisAngleRotation.toRotationMatrix());
    ASSERT_NO_THROW(axisAngleRotation.toChart());
    ASSERT_NO_THROW(axisAngleRotation.toRotationAxis(0));
    ASSERT_NO_THROW(axisAngleRotation * Vector3d({1.0, 0.0, 0.0}));
}

// Test constructor with axis and angle
TEST(AxisAngleRotationClass, ConstructorWithAxisAndAngle)
{
    Vector3d axis({1.0, 0.0, 0.0});
    double angle = M_PI / 2.0;
    AxisAngleRotation axisAngleRotation(axis, angle);

    ASSERT_NO_THROW(axisAngleRotation.quaternion());
    ASSERT_NO_THROW(axisAngleRotation.toRotationMatrix());
    ASSERT_NO_THROW(axisAngleRotation.toChart());
    ASSERT_NO_THROW(axisAngleRotation.toRotationAxis(0));
    ASSERT_NO_THROW(axisAngleRotation * Vector3d({1.0, 0.0, 0.0}));
}

// Test constructor with axis
TEST(AxisAngleRotationClass, ConstructorWithAxis)
{
    Vector3d axis({1.0, 0.0, 0.0});
    AxisAngleRotation axisAngleRotation(axis);

    ASSERT_NO_THROW(axisAngleRotation.quaternion());
    ASSERT_NO_THROW(axisAngleRotation.toRotationMatrix());
    ASSERT_NO_THROW(axisAngleRotation.toChart());
    ASSERT_NO_THROW(axisAngleRotation.toRotationAxis(0));
    ASSERT_NO_THROW(axisAngleRotation * Vector3d({1.0, 0.0, 0.0}));
}

// Test constructor with quaternion
TEST(AxisAngleRotationClass, ConstructorWithQuaternion)
{
    Quaternion quaternion(1.0, 0.0, 0.0, 0.0);
    AxisAngleRotation axisAngleRotation(quaternion);

    ASSERT_NO_THROW(axisAngleRotation.toRotationMatrix());
    ASSERT_NO_THROW(axisAngleRotation.toChart());
    ASSERT_NO_THROW(axisAngleRotation.toRotationAxis(0));
    ASSERT_NO_THROW(axisAngleRotation * Vector3d({1.0, 0.0, 0.0}));
}

// Test constructor with chart
TEST(AxisAngleRotationClass, ConstructorWithChart)
{
    Chart chart(Vector3d({1.0, 0.0, 0.0}), Vector3d({0.0, 1.0, 0.0}), Vector3d({0.0, 0.0, 1.0}));
    AxisAngleRotation axisAngleRotation(chart);

    ASSERT_NO_THROW(axisAngleRotation.quaternion());
    ASSERT_NO_THROW(axisAngleRotation.toRotationMatrix());
    ASSERT_NO_THROW(axisAngleRotation.toRotationAxis(0));
    ASSERT_NO_THROW(axisAngleRotation * Vector3d({1.0, 0.0, 0.0}));
}

// Test constructor between two vectors
TEST(AxisAngleRotationClass, ConstructorBetweenTwoVectors)
{
    Vector3d fromVector({1.0, 0.0, 0.0});
    Vector3d toVector({0.0, 1.0, 0.0});
    AxisAngleRotation axisAngleRotation(fromVector, toVector);

    ASSERT_NO_THROW(axisAngleRotation.quaternion());
    ASSERT_NO_THROW(axisAngleRotation.toRotationMatrix());
    ASSERT_NO_THROW(axisAngleRotation.toChart());
    ASSERT_NO_THROW(axisAngleRotation.toRotationAxis(0));
    ASSERT_NO_THROW(axisAngleRotation * Vector3d({1.0, 0.0, 0.0}));
}

// Test alignZ method
TEST(AxisAngleRotationClass, AlignZMethod)
{
    Vector3d axis({1.0, 0.0, 0.0});
    AxisAngleRotation axisAngleRotation = AxisAngleRotation::alignZ(axis);

    ASSERT_NO_THROW(axisAngleRotation.quaternion());
    ASSERT_NO_THROW(axisAngleRotation.toRotationMatrix());
    ASSERT_NO_THROW(axisAngleRotation.toChart());
    ASSERT_NO_THROW(axisAngleRotation.toRotationAxis(0));
    ASSERT_NO_THROW(axisAngleRotation * Vector3d({1.0, 0.0, 0.0}));
}

// Test alignYZ method
TEST(AxisAngleRotationClass, AlignYZMethod)
{
    Vector3d yAxis({0.0, 1.0, 0.0});
    Vector3d zAxis({0.0, 0.0, 1.0});
    AxisAngleRotation axisAngleRotation = AxisAngleRotation::alignYZ(yAxis, zAxis);

    ASSERT_NO_THROW(axisAngleRotation.quaternion());
    ASSERT_NO_THROW(axisAngleRotation.toRotationMatrix());
    ASSERT_NO_THROW(axisAngleRotation.toChart());
    ASSERT_NO_THROW(axisAngleRotation.toRotationAxis(0));
    ASSERT_NO_THROW(axisAngleRotation * Vector3d({1.0, 0.0, 0.0}));
}

// Test toRotationAxis method
TEST(AxisAngleRotationClass, ToRotationAxisMethod)
{
    Vector3d axis({1.0, 0.0, 0.0});
    AxisAngleRotation axisAngleRotation(axis);

    Vector3d rotationAxis = axisAngleRotation.toRotationAxis(0);

    ASSERT_EQ(rotationAxis, Vector3d({1.0, 0.0, 0.0}));
}

// Test multiplication operator
TEST(AxisAngleRotationClass, MultiplicationOperator)
{
    Vector3d axis1({1.0, 0.0, 0.0});
    AxisAngleRotation axisAngleRotation1(axis1);

    Vector3d axis2({0.0, 1.0, 0.0});
    AxisAngleRotation axisAngleRotation2(axis2);

    AxisAngleRotation result = axisAngleRotation1 * axisAngleRotation2;

    ASSERT_NO_THROW(result.quaternion());
    ASSERT_NO_THROW(result.toRotationMatrix());
    ASSERT_NO_THROW(result.toChart());
    ASSERT_NO_THROW(result.toRotationAxis(0));
    ASSERT_NO_THROW(result * Vector3d({1.0, 0.0, 0.0}));
}
