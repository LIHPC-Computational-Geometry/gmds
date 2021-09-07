#ifndef GMDS_POINTTESTSUITE_H
#define GMDS_POINTTESTSUITE_H

#include "gtest/gtest.h"
#include <gmds/math/Point.h>

#include <string>


TEST(PointClass, Setter)
{
    gmds::math::Point p(1,2,3);
    gmds::TCoord z = p.Z();
    ASSERT_NEAR(p.X(), 1, 1e-6);
    ASSERT_NEAR(p.Y(), 2, 1e-6);
    ASSERT_NEAR(z, 3, 1e-6);
}
TEST(PointClass, Setter2)
{
    gmds::math::Point p(1,2,3);
    p.X()=2*p.X();
    p.Y()=3*p.Y();
    p.Z()=0*p.Z();

    ASSERT_NEAR(p.X(), 2, 1e-6);
    ASSERT_NEAR(p.Y(), 6, 1e-6);
    ASSERT_NEAR(p.Z(), 0, 1e-6);
}
TEST(PointClass, SetXYZ)
{
    gmds::math::Point p(1,2,3);
    p.setXYZ(4,4,4);


    ASSERT_NEAR(p.X(), 4, 1e-6);
    ASSERT_NEAR(p.Y(), 4, 1e-6);
    ASSERT_NEAR(p.Z(), 4, 1e-6);
}
TEST(PointClass, Add)
{
    gmds::math::Point p1(1,2,3);
    gmds::math::Point p2(2,1,4);

    ASSERT_NEAR((p1+p2).X(), 3, 1e-6);
    ASSERT_NEAR((p1+p2).Y(), 3, 1e-6);
    ASSERT_NEAR((p1+p2).Z(), 7, 1e-6);
}

TEST(PointClass, Minus)
{
    gmds::math::Point p1(1,2,3);
    gmds::math::Point p2(2,1,3);

    ASSERT_NEAR((p1-p2).X(), -1, 1e-6);
    ASSERT_NEAR((p1-p2).Y(), 1, 1e-6);
    ASSERT_NEAR((p1-p2).Z(), 0, 1e-6);
}
TEST(PointClass, Multi)
{
    gmds::math::Point p(1,2,3);
    gmds::math::Point p2 = 2*p;
    gmds::math::Point p3 = p*3;

    ASSERT_NEAR(p2.X(),2, 1e-6);
    ASSERT_NEAR(p2.Y(),4, 1e-6);
    ASSERT_NEAR(p2.Z(),6, 1e-6);

    ASSERT_NEAR(p3.X(),3, 1e-6);
    ASSERT_NEAR(p3.Y(),6, 1e-6);
    ASSERT_NEAR(p3.Z(),9, 1e-6);
}

TEST(PointClass, EqualAndDiff)
{
    gmds::math::Point p1(2,1,3);
    gmds::math::Point p2(2,1,3);
    gmds::math::Point p3(2,1,4);

    ASSERT_TRUE(p1==p2);
    ASSERT_TRUE(p1==p1);
    ASSERT_FALSE(p1!=p2);
    ASSERT_FALSE(p1!=p1);
    ASSERT_FALSE(p1<p1);
    ASSERT_TRUE(p1<p3);
    ASSERT_TRUE(p1<=p1);

    ASSERT_EQ(p1.distance(p1),0);
    ASSERT_EQ(p1.distance(p2),0);
    ASSERT_EQ(p1.distance(p3),1);

    ASSERT_EQ(p1.distance2(p1),0);
    ASSERT_EQ(p1.distance2(p2),0);
    ASSERT_EQ(p1.distance2(p3),1);
}

TEST(PointClass, STREAM_OUT)
{
    gmds::math::Point p(2,1,0);
    std::stringstream s1;
    s1<<p;
    std::string s2="(2, 1, 0)";

    ASSERT_STREQ(s1.str().c_str(),s2.c_str());

}

TEST(PointClass, Comparison)
{
    gmds::math::Point p1(2,1,0);
    gmds::math::Point p2(2,2,0);

    ASSERT_TRUE(p1<p2);

}

TEST(PointClass, Addition)
{
    gmds::math::Point p1(2,1,-5);
    gmds::math::Point p2(2,10,6);

    gmds::math::Point p3 = p1+p2;
    ASSERT_EQ(p3.X(), p1.X()+p2.X());
    ASSERT_EQ(p3.Y(), p1.Y()+p2.Y());
    ASSERT_EQ(p3.Z(), p1.Z()+p2.Z());

}


#endif //GMDS_POINTTESTSUITE_H
