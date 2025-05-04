//
// Created by rochec on 28/08/24.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/math/BezierCurve.h>
#include <fstream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(BezierCurveClass, testBasic)
{
	std::vector<math::Point> ctrl_pts;

	ctrl_pts.push_back({0,0,0});
	ctrl_pts.push_back({0.33,0.3,0});
	ctrl_pts.push_back({0.66,-0.5,0});
	ctrl_pts.push_back({1,0,0});

	math::BezierCurve bezier_curve(ctrl_pts);

	math::Point p = bezier_curve(0.2);

	ASSERT_NEAR(p.X(), 0.19808, pow(10,-6));
	ASSERT_NEAR(p.Y(), 0.067199999999999996, pow(10,-6));

}
/*----------------------------------------------------------------------------*/
