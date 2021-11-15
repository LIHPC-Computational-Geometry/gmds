/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/math/DiscretizationScheme1D.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(Discretization1DTest, test_uniform) {
    math::Point ori(0,0,0);
    math::Point des(10,0,0);
    math::DiscretizationScheme1DUniform d(ori ,des,11);
    for(auto i=0;i<11;i++)
        EXPECT_EQ(d(i).X(),i);
}
/*----------------------------------------------------------------------------*/
TEST(Discretization1DTest, test_geo) {
    math::Point ori(0,0,0);
    math::Point des(10,0,0);
    math::DiscretizationScheme1DGeometric d(.5,ori ,des,11);
    EXPECT_FLOAT_EQ(d(0).X(),0);
    EXPECT_FLOAT_EQ(d(1).X(),5);
    EXPECT_FLOAT_EQ(d(2).X(),7.5);
    EXPECT_FLOAT_EQ(d(3).X(),8.75);
    EXPECT_FLOAT_EQ(d(4).X(),9.375);
    EXPECT_FLOAT_EQ(d(5).X(),9.6875);
    EXPECT_FLOAT_EQ(d(6).X(),9.84375);
    EXPECT_FLOAT_EQ(d(7).X(),9.921875);
    EXPECT_FLOAT_EQ(d(8).X(),9.96094);
    EXPECT_FLOAT_EQ(d(9).X(),9.98047);
    EXPECT_FLOAT_EQ(d(10).X(),9.9902344);
    EXPECT_FLOAT_EQ(d(11).X(),10);
}

/*----------------------------------------------------------------------------*/
TEST(Discretization1DTest, test_geo2) {
    math::Point ori(10,0);
    math::Point des(0,0,0);
    math::DiscretizationScheme1DGeometric d(.5,ori ,des,11);
    EXPECT_FLOAT_EQ(d(0).X(),10);
    EXPECT_FLOAT_EQ(d(1).X(),5);
    EXPECT_FLOAT_EQ(d(2).X(),2.5);
    EXPECT_FLOAT_EQ(d(3).X(),1.25);
    EXPECT_FLOAT_EQ(d(4).X(),0.625);
    EXPECT_FLOAT_EQ(d(5).X(),0.3125);
    EXPECT_FLOAT_EQ(d(6).X(),0.15625);
    EXPECT_FLOAT_EQ(d(7).X(),0.078125);
    EXPECT_FLOAT_EQ(d(8).X(),0.0390625);
    EXPECT_FLOAT_EQ(d(9).X(),0.01953125);
    EXPECT_FLOAT_EQ(d(10).X(),0.009765625);
    EXPECT_FLOAT_EQ(d(11).X(),0.0);
}
/*----------------------------------------------------------------------------*/
TEST(Discretization1DTest, test_geo3) {
    math::Point ori(0,0,0);
    math::Point des(10,0,0);
    math::DiscretizationScheme1DGeometric d(.5,ori ,des,11);
    d.setInverse(true);
    EXPECT_FLOAT_EQ(d(0).X(),0);
    EXPECT_FLOAT_EQ(d(10).X(),5);
    EXPECT_FLOAT_EQ(d(9).X(),2.5);
    EXPECT_FLOAT_EQ(d(8).X(),1.25);
    EXPECT_FLOAT_EQ(d(7).X(),0.625);
    EXPECT_FLOAT_EQ(d(6).X(),0.3125);
    EXPECT_FLOAT_EQ(d(5).X(),0.15625);
    EXPECT_FLOAT_EQ(d(4).X(),0.078125);
    EXPECT_FLOAT_EQ(d(3).X(),0.0390625);
    EXPECT_FLOAT_EQ(d(2).X(),0.01953125);
    EXPECT_FLOAT_EQ(d(1).X(),0.009765625);
    EXPECT_FLOAT_EQ(d(11).X(),10);
}
