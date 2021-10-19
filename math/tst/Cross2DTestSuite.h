/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include <gmds/math/Cross2D.h>
#include <gmds/math/Numerics.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(Cross2DTest, initVec1) {
    math::Vector3d v1(1, 0, 0);
    math::Vector3d v2(0, 1, 0);
    math::Cross2D c(v1, v2);

    EXPECT_EQ(0, c.referenceAngle());

    EXPECT_FALSE(c.hasStoredComponentVectors());
    c.computeComponentVectors();
    std::vector<math::Vector3d> cvc = c.componentVectors();
    EXPECT_TRUE(c.hasStoredComponentVectors());
    EXPECT_EQ(cvc.size(), 4);
    EXPECT_NEAR(cvc[0].X(),v1.X(),1e-9);
    EXPECT_NEAR(cvc[0].Y(),v1.Y(),1e-9);
    EXPECT_NEAR(cvc[0].Z(),v1.Z(),1e-9);

    EXPECT_NEAR(cvc[1].X(),v2.X(),1e-9);
    EXPECT_NEAR(cvc[1].Y(),v2.Y(),1e-9);
    EXPECT_NEAR(cvc[1].Z(),v2.Z(),1e-9);

    EXPECT_NEAR(cvc[2].X(),-v1.X(),1e-9);
    EXPECT_NEAR(cvc[2].Y(),-v1.Y(),1e-9);
    EXPECT_NEAR(cvc[2].Z(),-v1.Z(),1e-9);

    EXPECT_NEAR(cvc[3].X(),-v2.X(),1e-9);
    EXPECT_NEAR(cvc[3].Y(),-v2.Y(),1e-9);
    EXPECT_NEAR(cvc[3].Z(),-v2.Z(),1e-9);

}

/*----------------------------------------------------------------------------*/
TEST(Cross2DTest, initVec2) {
    math::Vector3d v1(0, 1, 0);
    math::Vector3d v2(-1, 0, 0);
    math::Cross2D c(v1, v2);

    EXPECT_EQ(0, c.referenceAngle());
}

/*----------------------------------------------------------------------------*/
TEST(Cross2DTest, initVec3) {
    math::Vector3d v1(-1, 0, 0);
    math::Vector3d v2(0, -1, 0);
    math::Cross2D c(v1, v2);

    EXPECT_EQ(0, c.referenceAngle());
}

/*----------------------------------------------------------------------------*/
TEST(Cross2DTest, initVec4) {
    math::Vector3d v1(1, 0, 0);
    math::Vector3d v2(0, -1, 0);
    math::Cross2D c(v1, v2);

    EXPECT_EQ(0, c.referenceAngle());
}

/*----------------------------------------------------------------------------*/
TEST(Cross2DTest, initVecFail) {
    math::Vector3d v1(1, 0, 0);
    math::Vector3d v2(0.1, 1, 0);
    EXPECT_THROW({
                     math::Cross2D c(v1, v2);

                 }, GMDSException);

}

/*----------------------------------------------------------------------------*/
TEST(Cross2DTest, initVectorsVsRef) {
    math::Vector3d v_ref(1, 0, 0);
    math::Cross2D c(v_ref);

    EXPECT_EQ(0, c.referenceAngle());

    math::Vector3d ref = c.referenceVector();

    EXPECT_EQ(1, ref.X());
    EXPECT_EQ(0, ref.Y());
    EXPECT_EQ(0, ref.Z());

    std::vector<math::Vector3d> c_vecs = c.componentVectors();

    EXPECT_EQ(4, c_vecs.size());

    for (unsigned int i = 0; i < c_vecs.size(); i++) {
        math::Vector3d current_vec = c_vecs[i];
        TCoord current_angle = current_vec.angle(v_ref);
        EXPECT_NEAR(0.0, math::modulo(current_angle, math::Constants::PIDIV2), 1e-12);
    }
}

/*----------------------------------------------------------------------------*/
TEST(Cross2DTest, Cross2DIndex1) {
    math::Cross2D c1(0);
    math::Cross2D c2(0);
    math::Cross2D c3(0);

    int index = math::Cross2D::index(c1, c2, c3);
    EXPECT_EQ(0, index);
}

/*----------------------------------------------------------------------------*/
TEST(Cross2DTest, Cross2DIndex2) {
    math::Cross2D c1(15);
    math::Cross2D c2(15);
    math::Cross2D c3(15);

    int index = math::Cross2D::index(c1, c2, c3);
    EXPECT_EQ(0, index);
}

/*----------------------------------------------------------------------------*/
TEST(Cross2DTest, Cross2DIndex3) {
    math::Cross2D c1(math::Constants::PI);
    math::Cross2D c2(math::Constants::PI * 1.1);
    math::Cross2D c3(math::Constants::PI);

    int index = math::Cross2D::index(c1, c2, c3);
    EXPECT_EQ(0, index);
}

/*----------------------------------------------------------------------------*/
TEST(Cross2DTest, Cross2DIndex4) {
    math::Cross2D c1(1);
    math::Cross2D c2(-2.1);
    math::Cross2D c3(-1.3);

    int index = math::Cross2D::index(c1, c2, c3);
    EXPECT_EQ(0, index);
}

/*----------------------------------------------------------------------------*/
TEST(Cross2DTest, Cross2DIndex5) {
    math::Cross2D c1(math::Constants::PIDIV2);
    math::Cross2D c2(0);
    math::Cross2D c3(math::Constants::PIDIV2 * 3);

    int index = math::Cross2D::index(c3, c2, c1);
    EXPECT_EQ(1, index);
}

/*----------------------------------------------------------------------------*/
TEST(Cross2DTest, Cross2DIndex6) {
    math::Cross2D c1(math::Constants::PIDIV2);
    math::Cross2D c2(0);
    math::Cross2D c3(math::Constants::PIDIV2 * 3);

    int index = math::Cross2D::index(c1, c2, c3);
    EXPECT_EQ(-1, index);
}
