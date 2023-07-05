//
// Created by rochec on 28/10/2022.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/math/BezierSurface.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(BezierSurfaceClass, testBasic_1)
{
	Array2D<math::Point> ctrl_pts(2,2);

	ctrl_pts(0,0) = {0,0,0};
	ctrl_pts(0,1) = {0,1,0};
	ctrl_pts(1,0) = {1,0,0};
	ctrl_pts(1,1) = {1,1,0};

	math::BezierSurface bezier_surface(ctrl_pts);

	math::Point p = bezier_surface(0.5, 0.5);

	ASSERT_NEAR(p.X(), 0.5, pow(10,-6));
	ASSERT_NEAR(p.Y(), 0.5, pow(10,-6));

	Array2D<math::Point> pts = bezier_surface.getDiscretization(10,10);

	ASSERT_NEAR(pts(5,5).X(), 0.5, pow(10,-6));
	ASSERT_NEAR(pts(5,5).Y(), 0.5, pow(10,-6));
}
/*----------------------------------------------------------------------------*/
TEST(BezierSurfaceClass, testBasic_2)
{
	Array2D<math::Point> ctrl_pts(3,2);

	ctrl_pts(0,0) = {0,0,4};
	ctrl_pts(0,1) = {0,0,0};

	ctrl_pts(1,0) = {0.5,2,4};
	ctrl_pts(1,1) = {0.5,8,0};

	ctrl_pts(2,0) = {1,0,4};
	ctrl_pts(2,1) = {1,0,0};


	math::BezierSurface bezier_surface(ctrl_pts);
	Array2D<math::Point> pts = bezier_surface.getDiscretization(40,40);

	// Unit tests
	{
		ASSERT_NEAR(pts(0, 0).X(), 0, pow(10,-6));
		ASSERT_NEAR(pts(0, 0).Y(), 0, pow(10,-6));
		ASSERT_NEAR(pts(0, 0).Z(), 4, pow(10,-6));
		ASSERT_NEAR(pts(0, 10).X(), 0, pow(10,-6));
		ASSERT_NEAR(pts(0, 10).Y(), 0, pow(10,-6));
		ASSERT_NEAR(pts(0, 10).Z(), 3, pow(10,-6));
		ASSERT_NEAR(pts(0, 20).X(), 0, pow(10,-6));
		ASSERT_NEAR(pts(0, 20).Y(), 0, pow(10,-6));
		ASSERT_NEAR(pts(0, 20).Z(), 2, pow(10,-6));
		ASSERT_NEAR(pts(0, 30).X(), 0, pow(10,-6));
		ASSERT_NEAR(pts(0, 30).Y(), 0, pow(10,-6));
		ASSERT_NEAR(pts(0, 30).Z(), 1, pow(10,-6));
		ASSERT_NEAR(pts(0, 40).X(), 0, pow(10,-6));
		ASSERT_NEAR(pts(0, 40).Y(), 0, pow(10,-6));
		ASSERT_NEAR(pts(0, 40).Z(), 0, pow(10,-6));
		ASSERT_NEAR(pts(10, 0).X(), 0.25, pow(10,-6));
		ASSERT_NEAR(pts(10, 0).Y(), 0.75, pow(10,-6));
		ASSERT_NEAR(pts(10, 0).Z(), 4, pow(10,-6));
		ASSERT_NEAR(pts(10, 10).X(), 0.25, pow(10,-6));
		ASSERT_NEAR(pts(10, 10).Y(), 1.3125, pow(10,-6));
		ASSERT_NEAR(pts(10, 10).Z(), 3, pow(10,-6));
		ASSERT_NEAR(pts(10, 20).X(), 0.25, pow(10,-6));
		ASSERT_NEAR(pts(10, 20).Y(), 1.875, pow(10,-6));
		ASSERT_NEAR(pts(10, 20).Z(), 2, pow(10,-6));
		ASSERT_NEAR(pts(10, 30).X(), 0.25, pow(10,-6));
		ASSERT_NEAR(pts(10, 30).Y(), 2.4375, pow(10,-6));
		ASSERT_NEAR(pts(10, 30).Z(), 1, pow(10,-6));
		ASSERT_NEAR(pts(10, 40).X(), 0.25, pow(10,-6));
		ASSERT_NEAR(pts(10, 40).Y(), 3, pow(10,-6));
		ASSERT_NEAR(pts(10, 40).Z(), 0, pow(10,-6));
		ASSERT_NEAR(pts(20, 0).X(), 0.5, pow(10,-6));
		ASSERT_NEAR(pts(20, 0).Y(), 1, pow(10,-6));
		ASSERT_NEAR(pts(20, 0).Z(), 4, pow(10,-6));
		ASSERT_NEAR(pts(20, 10).X(), 0.5, pow(10,-6));
		ASSERT_NEAR(pts(20, 10).Y(), 1.75, pow(10,-6));
		ASSERT_NEAR(pts(20, 10).Z(), 3, pow(10,-6));
		ASSERT_NEAR(pts(20, 20).X(), 0.5, pow(10,-6));
		ASSERT_NEAR(pts(20, 20).Y(), 2.5, pow(10,-6));
		ASSERT_NEAR(pts(20, 20).Z(), 2, pow(10,-6));
		ASSERT_NEAR(pts(20, 30).X(), 0.5, pow(10,-6));
		ASSERT_NEAR(pts(20, 30).Y(), 3.25, pow(10,-6));
		ASSERT_NEAR(pts(20, 30).Z(), 1, pow(10,-6));
		ASSERT_NEAR(pts(20, 40).X(), 0.5, pow(10,-6));
		ASSERT_NEAR(pts(20, 40).Y(), 4, pow(10,-6));
		ASSERT_NEAR(pts(20, 40).Z(), 0, pow(10,-6));
		ASSERT_NEAR(pts(30, 0).X(), 0.75, pow(10,-6));
		ASSERT_NEAR(pts(30, 0).Y(), 0.75, pow(10,-6));
		ASSERT_NEAR(pts(30, 0).Z(), 4, pow(10,-6));
		ASSERT_NEAR(pts(30, 10).X(), 0.75, pow(10,-6));
		ASSERT_NEAR(pts(30, 10).Y(), 1.3125, pow(10,-6));
		ASSERT_NEAR(pts(30, 10).Z(), 3, pow(10,-6));
		ASSERT_NEAR(pts(30, 20).X(), 0.75, pow(10,-6));
		ASSERT_NEAR(pts(30, 20).Y(), 1.875, pow(10,-6));
		ASSERT_NEAR(pts(30, 20).Z(), 2, pow(10,-6));
		ASSERT_NEAR(pts(30, 30).X(), 0.75, pow(10,-6));
		ASSERT_NEAR(pts(30, 30).Y(), 2.4375, pow(10,-6));
		ASSERT_NEAR(pts(30, 30).Z(), 1, pow(10,-6));
		ASSERT_NEAR(pts(30, 40).X(), 0.75, pow(10,-6));
		ASSERT_NEAR(pts(30, 40).Y(), 3, pow(10,-6));
		ASSERT_NEAR(pts(30, 40).Z(), 0, pow(10,-6));
		ASSERT_NEAR(pts(40, 0).X(), 1, pow(10,-6));
		ASSERT_NEAR(pts(40, 0).Y(), 0, pow(10,-6));
		ASSERT_NEAR(pts(40, 0).Z(), 4, pow(10,-6));
		ASSERT_NEAR(pts(40, 10).X(), 1, pow(10,-6));
		ASSERT_NEAR(pts(40, 10).Y(), 0, pow(10,-6));
		ASSERT_NEAR(pts(40, 10).Z(), 3, pow(10,-6));
		ASSERT_NEAR(pts(40, 20).X(), 1, pow(10,-6));
		ASSERT_NEAR(pts(40, 20).Y(), 0, pow(10,-6));
		ASSERT_NEAR(pts(40, 20).Z(), 2, pow(10,-6));
		ASSERT_NEAR(pts(40, 30).X(), 1, pow(10,-6));
		ASSERT_NEAR(pts(40, 30).Y(), 0, pow(10,-6));
		ASSERT_NEAR(pts(40, 30).Z(), 1, pow(10,-6));
		ASSERT_NEAR(pts(40, 40).X(), 1, pow(10,-6));
		ASSERT_NEAR(pts(40, 40).Y(), 0, pow(10,-6));
		ASSERT_NEAR(pts(40, 40).Z(), 0, pow(10,-6));
	}

}
/*----------------------------------------------------------------------------*/