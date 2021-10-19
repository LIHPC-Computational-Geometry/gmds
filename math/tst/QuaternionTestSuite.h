/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/math/Quaternion.h>
#include <gmds/math/Cross2D.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(QuaternionTest,init) {
  math::Quaternion q;
  EXPECT_DOUBLE_EQ(1, q.X());
  EXPECT_DOUBLE_EQ(0, q.I());
  EXPECT_DOUBLE_EQ(0, q.J());
  EXPECT_DOUBLE_EQ(0, q.K());
}
/*----------------------------------------------------------------------------*/
TEST(QuaternionTest,opposite) {
  math::Quaternion q1(1,2,1,1);
  math::Quaternion q2 = q1.opposite();

  EXPECT_DOUBLE_EQ(-q1.X(), q2.X());
  EXPECT_DOUBLE_EQ(-q1.I(), q2.I());
  EXPECT_DOUBLE_EQ(-q1.J(), q2.J());
  EXPECT_DOUBLE_EQ(-q1.K(), q2.K());
}

/*----------------------------------------------------------------------------*/
TEST(QuaternionTest,conjugate) {
  math::Quaternion q1(1,2,1,1);
  math::Quaternion q2 = q1.conjugate();

  EXPECT_DOUBLE_EQ(q1.X(), q2.X());
  EXPECT_DOUBLE_EQ(-q1.I(), q2.I());
  EXPECT_DOUBLE_EQ(-q1.J(), q2.J());
  EXPECT_DOUBLE_EQ(-q1.K(), q2.K());
}

/*----------------------------------------------------------------------------*/
TEST(QuaternionTest,chart) {
  math::Quaternion q1(1,2,1,1);
  math::Chart c1(q1);
  
  math::Quaternion q2(c1);
  
  EXPECT_DOUBLE_EQ(q1.X(), q2.X());
  EXPECT_DOUBLE_EQ(q1.I(), q2.I());
  EXPECT_DOUBLE_EQ(q1.J(), q2.J());
  EXPECT_DOUBLE_EQ(q1.K(), q2.K());
}
/*----------------------------------------------------------------------------*/
TEST(QuaternionTest, test_local_face) {

    math::Vector3d x1(0.99018058993906, -0.13979404255477, 0.00015803205361307);
    math::Vector3d y1(0.13979401007176, 0.99018061974691, 0.00022989646223284);
    math::Vector3d z1(0.000188619, 0.000205547, -1);

    math::Vector3d x2(-0.71667006172124, -0.69741232511488, -0.00026722711938157);
    math::Vector3d y2(0.69741237821558, -0.71667010968904, -1.7222600781856e-05);
    math::Vector3d z2(0.000179502, 0.00019871, -1);

    math::Vector3d x3(-0.45598714581052, 0.88998635453604, 0.00010768114638167);
    math::Vector3d y3(-0.86100802396785, -0.44113997782073, -0.000284806450452);
    math::Vector3d z3(0.000212904, 0.000230073, -1);

    math::Chart c1(x1,y1,z1);
    math::Chart c2(x2,y2,z2);
    math::Chart c3(x3,y3,z3);

    math::Chart::Mapping m12(c1,c2);
    math::Chart::Mapping m23(c2,c3);
    math::Chart::Mapping m31(c3,c1);

    math::Quaternion q1(c1);
    math::Quaternion q2(c2);
    math::Quaternion q3(c3);

    ASSERT_FALSE((m31*m23*m12).isIdentity());


    ASSERT_EQ(math::Quaternion::testSingularity(q1,q2,q3,q3),0);

    math::Cross2D cr1(x1,y1);
    math::Cross2D cr2(x2,y2);
    math::Cross2D cr3(x3,y3);

    int index = math::Cross2D::index(cr1,cr2,cr3);
    ASSERT_TRUE(index!=0);

    math::Vector3d vi = cr1.referenceVector();
    math::Vector3d vj = cr2.referenceVector();
    math::Vector3d vk = cr3.referenceVector();


    double wij = vi.orientedAngle(vj,c1.Z());
    double wjk = vj.orientedAngle(vk,c2.Z());
    double wki = vk.orientedAngle(vi,c3.Z());

    int index2 = std::round((wij+wjk+wki)/math::Constants::PI2);
    ASSERT_TRUE(index2!=0.0);

}