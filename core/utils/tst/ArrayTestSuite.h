/*----------------------------------------------------------------------------*/
#ifndef GMDS_ARRAY_TESTSUITE_H
#define GMDS_ARRAY_TESTSUITE_H
/*----------------------------------------------------------------------------*/
#include "gtest/gtest.h"
#include <unit_test_config.h>
#include <gmds/utils/Array.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(ArrayTestSuite, test_int){
    Array2D<int> array(10,4);
    ASSERT_EQ(array.nbLines(),10);
    ASSERT_EQ(array.nbColumns(),4);
    ASSERT_EQ(array(1,1),0);
    array(1,1)=1;
    ASSERT_EQ(array(1,1),1);
}
/*----------------------------------------------------------------------------*/
TEST(ArrayTestSuite, test_double3D){
    Array3D<double> array(3,3,3);
    ASSERT_EQ(array.nbElements(1),3);
    ASSERT_EQ(array(1,1,1),0);
    array(1,1,1)=1;
    ASSERT_EQ(array(1,1,1),1);
}
/*----------------------------------------------------------------------------*/
TEST(ArrayTestSuite, test_TryArray){
    TriArray<double> a(3);
    ASSERT_EQ(a.nbElements(),6);
    ASSERT_EQ(a.size(),3);
    a(0,0)=1;
    a(0,1)=2;
    a(0,2)=3;

    a(1,0)=4;
    a(1,1)=5;
    a(2,0)=6;

    ASSERT_EQ(a(2,0,0),1);
    ASSERT_EQ(a(0,2,0),3);
    ASSERT_EQ(a(0,0,2),6);

    ASSERT_EQ(a(1,1,0),2);
    ASSERT_EQ(a(1,0,1),4);
    ASSERT_EQ(a(0,1,1),5);

}
/*----------------------------------------------------------------------------*/
#endif //GMDS_ARRAY_TESTSUITE_H
/*----------------------------------------------------------------------------*/
