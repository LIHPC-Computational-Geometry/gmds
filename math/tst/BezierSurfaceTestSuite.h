//
// Created by rochec on 28/10/2022.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/math/BezierSurface.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(BezierSurfaceClass, testBasic)
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