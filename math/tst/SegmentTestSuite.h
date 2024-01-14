/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/math/Segment.h>
#include <gmds/math/Point.h>
/*----------------------------------------------------------------------------*/
using namespace gmds::math;
/*----------------------------------------------------------------------------*/

TEST(SegmentClass, Constructors)
{
    // Test default constructor
    {
        Segment seg;
        ASSERT_EQ(seg.getPoint(0), Point());
        ASSERT_EQ(seg.getPoint(1), Point());
    }

    // Test constructor with points
    {
        Point p1(1, 2, 3);
        Point p2(4, 5, 6);
        Segment seg(p1, p2);
        ASSERT_EQ(seg.getPoint(0), p1);
        ASSERT_EQ(seg.getPoint(1), p2);
    }
}

TEST(SegmentClass, SetMethod)
{
    Point p1(1, 2, 3);
    Point p2(4, 5, 6);
    Segment seg;
    seg.set(p1, p2);
    ASSERT_EQ(seg.getPoint(0), p1);
    ASSERT_EQ(seg.getPoint(1), p2);
}

TEST(SegmentClass, Getters)
{
    Point p1(1, 2, 3);
    Point p2(4, 5, 6);
    Segment seg(p1, p2);

    // Test getPoint
    ASSERT_EQ(seg.getPoint(0), p1);
    ASSERT_EQ(seg.getPoint(1), p2);

    // Test getDir
    Vector3d dir = p2 - p1;
    ASSERT_EQ(seg.getDir(), dir);

    // Test getUnitVector
    Vector3d unitVector = dir.normalized();
    ASSERT_EQ(seg.getUnitVector(), unitVector);

    // Test computeCenter
    Point center = (p1 + p2) * 0.5;
    ASSERT_EQ(seg.computeCenter(), center);

    // Test computeLength
    TCoord length = p1.distance(p2);
    ASSERT_EQ(seg.computeLength(), length);
}

TEST(SegmentClass, AssignmentOperator)
{
    Point p1(1, 2, 3);
    Point p2(4, 5, 6);
    Segment seg1(p1, p2);
    Segment seg2;
    seg2 = seg1;
    ASSERT_EQ(seg2.getPoint(0), p1);
    ASSERT_EQ(seg2.getPoint(1), p2);
}

TEST(SegmentClass, ComparisonOperators)
{
    Point p1(1, 2, 3);
    Point p2(4, 5, 6);
    Point p3(7, 8, 9);
    Point p4(10, 11, 12);

    Segment seg1(p1, p2);
    Segment seg2(p1, p2);
    Segment seg3(p3, p4);

    // Test operator==
    ASSERT_TRUE(seg1 == seg2);
    ASSERT_FALSE(seg1 == seg3);

    // Test operator!=
    ASSERT_FALSE(seg1 != seg2);
    ASSERT_TRUE(seg1 != seg3);
}

TEST(SegmentClass, IsIn)
{
    Point p1(1, 2, 3);
    Point p2(4, 5, 6);
    Point p3(7, 8, 9);

    Segment seg(p1, p2);

    // Test isIn
    ASSERT_TRUE(seg.isIn(p1));
    ASSERT_TRUE(seg.isIn(p2));
    ASSERT_TRUE(seg.isIn(seg.computeCenter()));

    ASSERT_FALSE(seg.isIn(p3));
}

/*----------------------------------------------------------------------------*/
