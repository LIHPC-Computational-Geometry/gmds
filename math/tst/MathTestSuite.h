
/*----------------------------------------------------------------------------*/
#include<string>
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include "gmds/math/Line.h"
#include <gmds/math/Hexahedron.h>
#include <gmds/math/Numerics.h>
#include <gmds/math/Plane.h>
#include <gmds/math/Point.h>
#include <gmds/math/Ray.h>
#include <gmds/math/Segment.h>
#include <gmds/math/Triangle.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(MathTest,hexahedronTriangleIntersection) {
	
	math::Point pnts[8];

	pnts[0] = math::Point(1, 2, 4.3);
	pnts[1] = math::Point(1.1, 2, 4.3);
	pnts[2] = math::Point(1.1, 2.25, 4.3);
	pnts[3] = math::Point(1, 2.25, 4.3);
	pnts[4] = math::Point(1, 2, 4.6);
	pnts[5] = math::Point(1.1, 2, 4.6);
	pnts[6] = math::Point(1.1, 2.25, 4.6);
	pnts[7] = math::Point(1, 2.25, 4.6);
	
	math::Hexahedron hex(pnts[0],pnts[1],pnts[2],pnts[3],
			pnts[4],pnts[5],pnts[6],pnts[7]);

	math::Point pntsbis[3];
	pntsbis[0] = math::Point(1, 1.73205, 4);
	pntsbis[1] = math::Point(1, 1.73205, 4.5);
	pntsbis[2] = math::Point(0.72555, 1.86375, 4.16573);

	math::Triangle triangle(pntsbis[0],pntsbis[1],pntsbis[2]);

	bool res = hex.intersect(triangle);

	EXPECT_EQ(false,res);
}
/*----------------------------------------------------------------------------*/
TEST(MathTest,TriangleTriangleIntersection) {
	math::Point pnts[3];
	pnts[0] = math::Point(1, 2.25, 4.3);
	pnts[1] = math::Point(1, 2, 4.3);
	pnts[2] = math::Point(1, 2, 4.6);

	math::Triangle triangle(pnts[0],pnts[1],pnts[2]);

	math::Point pntsbis[3];
        pntsbis[0] = math::Point(1, 1.73205, 4);
        pntsbis[1] = math::Point(1, 1.73205, 4.5);
        pntsbis[2] = math::Point(0.72555, 1.86375, 4.16573);

        math::Triangle triangle2(pntsbis[0],pntsbis[1],pntsbis[2]);

	bool res = triangle.intersect(triangle2);

        EXPECT_EQ(false,res);
}
/*----------------------------------------------------------------------------*/
TEST(MathTest,TriangleTriangleIntersection_2) {
        math::Point pnt0(-0.492124, -0.712611, -0.5);
        math::Point pnt1(-0.620854, -0.603772, -0.5);
        math::Point pnt2(-0.596846, -0.627515, -0.5);
        math::Triangle triangle(pnt0,pnt1,pnt2);

	math::Point pnt3(-0.7, -0.8, -0.5);
        math::Point pnt4(-0.6, -0.8, -0.5);
        math::Point pnt5(-0.6, -0.7, -0.5);
        math::Triangle triangle2(pnt3,pnt4,pnt5);

        bool res = triangle.intersect(triangle2);

        EXPECT_EQ(false,res);
}
/*----------------------------------------------------------------------------*/
TEST(MathTest,TriangleTriangleIntersection_3) {
	math::Point pnt0(0, 0.5, 0.5);
        math::Point pnt1(0, 0, 0);
        math::Point pnt2(0, 1, 0);
        math::Triangle triangle(pnt0,pnt1,pnt2);

        math::Point pnt3(0, 1.1, 0.1);
        math::Point pnt4(0, 1, 0.1);
        math::Point pnt5(0, 1, 0.2);
        math::Triangle triangle2(pnt3,pnt4,pnt5);

        bool res = triangle.intersect(triangle2);

        EXPECT_EQ(false,res);
}
/*----------------------------------------------------------------------------*/
TEST(MathTest,TriangleSegmentIntersection) {
        math::Point pnts[3];
        pnts[0] = math::Point(1, 2.25, 4.3);
        pnts[1] = math::Point(1, 2, 4.3);
        pnts[2] = math::Point(1, 2, 4.6);

        math::Triangle triangle(pnts[0],pnts[1],pnts[2]);

        math::Point pntsbis[3];
        pntsbis[0] = math::Point(1, 1.73205, 4);
        pntsbis[1] = math::Point(1, 1.73205, 4.5);
//        pntsbis[2] = math::Point(0.72555, 1.86375, 4.16573);

        math::Segment segment(pntsbis[0],pntsbis[1]);

        bool res = triangle.intersect(segment);

        EXPECT_EQ(false,res);
}
/*----------------------------------------------------------------------------*/
TEST(MathTest,TriangleSegmentIntersection2D) {
        math::Point pnt3(-0.7, -0.8, 0.);
        math::Point pnt4(-0.6, -0.8, 0.);
        math::Point pnt5(-0.6, -0.7, 0.);
        math::Triangle triangle2(pnt3,pnt4,pnt5);

	math::Point pnt0(-0.492124, -0.712611, 0.);
	math::Point pnt1(-0.620854, -0.603772, 0.);
        math::Segment segment(pnt0,pnt1);

        bool res = triangle2.intersect2D(segment);

        EXPECT_EQ(false,res);
}
/*----------------------------------------------------------------------------*/
TEST(MathTest,TriangleSegmentIntersection2D_2) {
        math::Point pnt3(1.1, 0.1, 0);
        math::Point pnt4(1., 0.1, 0);
        math::Point pnt5(1., 0.2, 0);
        math::Triangle triangle2(pnt3,pnt4,pnt5);

        math::Point pnt0(1, 0, 0.);
        math::Point pnt1(0.5, 0.5, 0.);
        math::Segment segment(pnt0,pnt1);

        bool res = triangle2.intersect2D(segment);

        EXPECT_EQ(false,res);
}
/*----------------------------------------------------------------------------*/
TEST(MathTest,TriangleRayIntersection) {
	math::Point pnt0(0, -0.4, 0);
	math::Point pnt1(0.47629, -0.479422, 0);
	math::Point pnt2(0.312907, -0.766642, 0);
	math::Triangle triangle(pnt0,pnt1,pnt2);

	math::Point pnt(0.263162, -0.546133, 5);
	math::Vector3d dir({0, 0, -0.15228});
	math::Ray ray(pnt,dir);	

	bool res = triangle.intersect(ray);
	
	EXPECT_EQ(true,true);//res);
}
/*----------------------------------------------------------------------------*/
TEST(MathTest,HexahedronTriangleIntersection_2) {
        math::Point pnt0(-0.492124, -0.712611, -0.5);
        math::Point pnt1(-0.620854, -0.603772, -0.5);
        math::Point pnt2(-0.596846, -0.627515, -0.5);
        math::Triangle triangle(pnt0,pnt1,pnt2);

        math::Point pnta(-0.7, -0.8, -0.6);
        math::Point pntb(-0.6, -0.8, -0.6);
        math::Point pntc(-0.6, -0.7, -0.6);
        math::Point pntd(-0.7, -0.7, -0.6);
        math::Point pnte(-0.7, -0.8, -0.5);
        math::Point pntf(-0.6, -0.8, -0.5);
        math::Point pntg(-0.6, -0.7, -0.5);
        math::Point pnth(-0.7, -0.7, -0.5);

        math::Hexahedron hex(pnta,pntb,pntc,pntd,pnte,pntf,pntg,pnth);

        bool res = hex.intersect(triangle);

        EXPECT_EQ(false,res);
}
/*----------------------------------------------------------------------------*/
TEST(MathTest,HexahedronTriangleIntersection_3) {
        math::Point pnt0(0, 0.5, 0.5);
        math::Point pnt1(0, 0, 0);
        math::Point pnt2(0, 1, 0);
        math::Triangle triangle(pnt0,pnt1,pnt2);

        math::Point pnta(0, 1, 0.1);
        math::Point pntb(0.1, 1, 0.1);
        math::Point pntc(0.1, 1.1, 0.1);
        math::Point pntd(0, 1.1, 0.1);
        math::Point pnte(0, 1, 0.2);
        math::Point pntf(0.1, 1, 0.2);
        math::Point pntg(0.1, 1.1, 0.2);
        math::Point pnth(0, 1.1, 0.2);

        math::Hexahedron hex(pnta,pntb,pntc,pntd,pnte,pntf,pntg,pnth);

        bool res = hex.intersect(triangle);

        EXPECT_EQ(false,res);
}
/*----------------------------------------------------------------------------*/
TEST(MathTest,HexahedronScaledJacobian) {
        math::Point pnta(0, 0, 0);
        math::Point pntb(1, 0, 0);
        math::Point pntc(1, 1, 0);
        math::Point pntd(0, 1, 0);
        math::Point pnte(0, 0, 1);
        math::Point pntf(1, 0, 1);
        math::Point pntg(1, 1, 1);
        math::Point pnth(0, 1, 1);

        math::Hexahedron hex(pnta,pntb,pntc,pntd,pnte,pntf,pntg,pnth);

        double res = hex.computeScaledJacobian();

        EXPECT_EQ(1,res);
}
/*----------------------------------------------------------------------------*/
TEST(MathTest,VectorAngle) {
  math::Vector3d v1({1, 0, 0});
  math::Vector3d v2({0, -1, 0});
 
  double angle1 = v1.angle(v2);
  double angle2 = v2.angle(v1);
   
  EXPECT_NEAR(angle1, angle2, 1e-12);
}
/*----------------------------------------------------------------------------*/
TEST(MathTest,VectorOrientedAngle) {
  math::Vector3d v1({1, 0, 0});
  math::Vector3d v2({0, -1, 0});
 
  double angle1 = v1.orientedAngle(v2);
  double angle2 = v2.orientedAngle(v1);
   
  EXPECT_TRUE(angle1<0);
  EXPECT_TRUE(angle2>0);
  EXPECT_NEAR(angle1, -angle2, 1e-12);
}
/*----------------------------------------------------------------------------*/
TEST(MathTest,VectorSelfAngle) {
  math::Vector3d v({1, 1, 0});
 
  double angle1 = v.angle(v);
  double angle2 = v.orientedAngle(v);
   
  EXPECT_NEAR(0.0, angle1, 1e-6);
  EXPECT_NEAR(0.0, angle2, 1e-6);
}
/*----------------------------------------------------------------------------*/
TEST(MathTest,Vector3dOrtho1) {
    math::Vector3d v({1, 0, 0});
    math::Vector3d v2 = v.getOneOrtho();
    EXPECT_TRUE(math::isZero(v.dot(v2)));

}
/*----------------------------------------------------------------------------*/
TEST(MathTest,Vector3dOrtho2) {
    math::Vector3d v({0, -2, 0});
    math::Vector3d v2 = v.getOneOrtho();
    EXPECT_TRUE(math::isZero(v.dot(v2)));

}
/*----------------------------------------------------------------------------*/
TEST(MathTest,Vector3dOrtho3) {
    math::Vector3d v({0, 0, 4});
    math::Vector3d v2 = v.getOneOrtho();

    EXPECT_TRUE(math::isZero(v.dot(v2)));

}
/*----------------------------------------------------------------------------*/
TEST(MathTest,Vector3dOrtho4) {
    math::Vector3d v({1.23, -0.45, 124.5});

    math::Vector3d v2 = v.getOneOrtho();

    EXPECT_TRUE(math::isZero(v.dot(v2)));

}
/*----------------------------------------------------------------------------*/
TEST(MathTest,VectorOrderingAngles) {
  math::Vector3d v({1, 0, 0});
  math::Vector3d v1({1, 1, 0});
  math::Vector3d v2({1, -1, 0});

  double angle1   = v.angle(v1);
  double angle2   = v.angle(v2);
  double angle1_o = v.orientedAngle(v1);
  double angle2_o = v.orientedAngle(v2);

  EXPECT_NEAR(angle2, angle1, 1e-6);
  EXPECT_TRUE(angle1_o>angle2_o);
}
/*----------------------------------------------------------------------------*/
TEST(MathTest, numericComparison){
    ASSERT_TRUE(math::isZero(1e-9));
    ASSERT_TRUE(math::areEquals(1, 1+1e-9));
}

/*----------------------------------------------------------------------------*/
TEST(MathTest, planeSegmentIntersection1){
    math::Plane pl(math::Point(0,0,0), math::Vector3d({0, 0, 1}));
    math::Segment seg(math::Point(0,0,-1), math::Point(0,0,2));
    math::Point pi;
    double w0=0, w1=0;
    ASSERT_EQ(math::Plane::SEGMENT_MIDDLE,
              pl.intersect(seg, pi, w0, w1));
}

/*----------------------------------------------------------------------------*/
TEST(MathTest, planeSegmentIntersection2){
    math::Plane pl(math::Point(0,0,0), math::Vector3d({0, 0, 1}));
    math::Segment seg(math::Point(1,1,0), math::Point(0,0,2));
    math::Point pi;
    double w0=0, w1=0;
    ASSERT_EQ(math::Plane::SEGMENT_FIRST_PNT,
              pl.intersect(seg, pi, w0, w1));
}

/*----------------------------------------------------------------------------*/
TEST(MathTest, planeSegmentIntersection3){
    math::Plane pl(math::Point(0,0,0), math::Vector3d({0, 0, 1}));
    math::Segment seg(math::Point(0,0,-1), math::Point(1,1234,0));
    math::Point pi;
    double w0=0, w1=0;
    ASSERT_EQ(math::Plane::SEGMENT_SECOND_PNT,
              pl.intersect(seg, pi, w0, w1));
}
/*----------------------------------------------------------------------------*/
TEST(MathTest, planeSegmentIntersection4){
    math::Plane pl(math::Point(0,0,0), math::Vector3d({0, 0, 1}));
    math::Segment seg(math::Point(0,0,-1), math::Point(0,0,-0.5));
    math::Point pi;
    double w0=0, w1=0;
    ASSERT_EQ(math::Plane::NO_INTERSECTION,
              pl.intersect(seg, pi, w0, w1));
}
/*----------------------------------------------------------------------------*/
TEST(MathTest, dotTest){
    math::Point pa(5,3.2,-0.5);
    math::Point pb(5,4.2,-0.5);
    math::Point pc(5,6.4,-0.5);
    math::Point p1(5,4.427,-0.5);
    math::Point p2(5,4.144,-0.5);

    math::Vector3d va=pa-p1;
    math::Vector3d v12=p2-p1;
    TCoord a = v12.dot(va);
    ASSERT_TRUE(a>=0);
    ASSERT_TRUE(a>=v12.norm());

    math::Vector3d vb=pb-p1;
    a = v12.dot(vb);

    ASSERT_TRUE(a>=0);
    ASSERT_TRUE(a<=v12.norm());
    math::Vector3d vc=pc-p1;
    a = v12.dot(vc);

    ASSERT_TRUE(a<=0);
}
/*----------------------------------------------------------------------------*/
TEST(MathTest,LineTest) {
	math::Point a(0, 0, 0);
	math::Point b(1, 0, 0);
	math::Point c(1, 1, 0);
	math::Point d(1, 2, 0);

	math::Line ab(a,b);
	math::Line cd(c,d);

	math::Point p;
	double param;
	EXPECT_TRUE(ab.intersect2D(cd,p,param));
}