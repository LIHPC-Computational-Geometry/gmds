//
// Created by rochec on 05/12/23.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/math/BezierHex.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(BezierHexClass, testBasic_1)
{
	Array3D<math::Point> ctrl_pts(2,2, 2);

	ctrl_pts(0,0,0) = {0,0,0};
	ctrl_pts(0,1,0) = {0,1,0};
	ctrl_pts(1,0,0) = {1,0,0};
	ctrl_pts(1,1,0) = {1,1,0};
	ctrl_pts(0,0,1) = {0,0,1};
	ctrl_pts(0,1,1) = {0,1,1};
	ctrl_pts(1,0,1) = {1,0,1};
	ctrl_pts(1,1,1) = {1,1,1};

	math::BezierHex bezier_hex(ctrl_pts);

	math::Point p = bezier_hex(0.5, 0.5, 0.5);

	ASSERT_NEAR(p.X(), 0.5, pow(10,-6));
	ASSERT_NEAR(p.Y(), 0.5, pow(10,-6));
	ASSERT_NEAR(p.Z(), 0.5, pow(10,-6));

	Array3D<math::Point> pts = bezier_hex.getDiscretization(10,10,10);

	ASSERT_NEAR(pts(5,5,5).X(), 0.5, pow(10,-6));
	ASSERT_NEAR(pts(5,5,5).Y(), 0.5, pow(10,-6));
	ASSERT_NEAR(pts(5,5,5).Z(), 0.5, pow(10,-6));
}
/*----------------------------------------------------------------------------*/
TEST(BezierHexClass, testBasic_2)
{
	Array3D<math::Point> ctrl_pts(2,2,3);

	ctrl_pts(0,0,0) = {0,0,0};
	ctrl_pts(0,1,0) = {0,1,0};
	ctrl_pts(1,0,0) = {1,0,0};
	ctrl_pts(1,1,0) = {1,1,0};
	ctrl_pts(0,0,1) = {0,0,0.5};
	ctrl_pts(0,1,1) = {0,1,0.5};
	ctrl_pts(1,0,1) = {1,0,0.5};
	ctrl_pts(1,1,1) = {1,1,0.5};
	ctrl_pts(0,0,2) = {0,0,1};
	ctrl_pts(0,1,2) = {0,1,1};
	ctrl_pts(1,0,2) = {1,0,1};
	ctrl_pts(1,1,2) = {1,1,1};

	math::BezierHex bezier_hex(ctrl_pts);

	math::Point p = bezier_hex(0.5, 0.5, 0.5);

	ASSERT_NEAR(p.X(), 0.5, pow(10,-6));
	ASSERT_NEAR(p.Y(), 0.5, pow(10,-6));
	ASSERT_NEAR(p.Z(), 0.5, pow(10,-6));

	Array3D<math::Point> pts = bezier_hex.getDiscretization(10,10,10);

	ASSERT_NEAR(pts(5,5,5).X(), 0.5, pow(10,-6));
	ASSERT_NEAR(pts(5,5,5).Y(), 0.5, pow(10,-6));
	ASSERT_NEAR(pts(5,5,5).Z(), 0.5, pow(10,-6));
}
/*----------------------------------------------------------------------------*/