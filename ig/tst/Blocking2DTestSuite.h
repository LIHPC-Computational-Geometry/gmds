/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Blocking2D.h>
#include <iostream>

/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(Blocking2DTestSuite, test_blocking2D_1)
{
    Blocking2D m;
    Node n1 = m.newBlockCorner(0,0);
    Node n2 = m.newBlockCorner(1,0);
    Node n3 = m.newBlockCorner(1,1);
    Node n4=  m.newBlockCorner(0,1);

    Blocking2D::Block b0 = m.newBlock(n1,n2,n3,n4);

    Node n5 = m.newBlockCorner(2,0,0);
    Node n6 = m.newBlockCorner(2,1.5,0);
    Blocking2D::Block b1 = m.newBlock(n2,n5,n6,n3);

    ASSERT_EQ(b0.id(), m.block(0).id());
    b0.seNbDiscretizationI(11);
    b0.seNbDiscretizationJ(11);
    b1.seNbDiscretizationI(11);
    b1.seNbDiscretizationI(11);

    ASSERT_EQ(n1.id(),b0.origin());
    ASSERT_FLOAT_EQ(b0.getUnitVectorI().X(),1.0);
    ASSERT_FLOAT_EQ(b0.getUnitVectorI().Y(),0.0);
    ASSERT_FLOAT_EQ(b0.getUnitVectorI().Z(),0.0);
    ASSERT_FLOAT_EQ(b0.getUnitVectorJ().X(),0.0);
    ASSERT_FLOAT_EQ(b0.getUnitVectorJ().Y(),1.0);
    ASSERT_FLOAT_EQ(b0.getUnitVectorJ().Z(),0.0);

    ASSERT_EQ(b0.getNbDiscretizationI(),11);
    ASSERT_EQ(b0.getNbDiscretizationJ(),11);

    math::Point p00 = b0(0,0).point();
    ASSERT_FLOAT_EQ(p00.X(),0);
    ASSERT_FLOAT_EQ(p00.Y(),0);

    math::Point p04 = b0(0,4).point();
    ASSERT_FLOAT_EQ(p04.X(),0);
    ASSERT_FLOAT_EQ(p04.Y(),0.4);

    math::Point p31 = b0(3,1).point();
    ASSERT_FLOAT_EQ(p31.X(),0.3);
    ASSERT_FLOAT_EQ(p31.Y(),0.1);

    math::Point p76 = b0(7,6).point();
    ASSERT_FLOAT_EQ(p76.X(),0.7);
    ASSERT_FLOAT_EQ(p76.Y(),0.6);

}