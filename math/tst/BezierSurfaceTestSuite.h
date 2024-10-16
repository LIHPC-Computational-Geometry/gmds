//
// Created by rochec on 28/10/2022.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/math/BezierSurface.h>
#include <fstream>
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
TEST(BezierSurfaceClass, testBasic_Manuscrit)
{

	Array2D<math::Point> ctrl_pts_1(3,3);

	ctrl_pts_1(0,0) = {0,0,0};
	ctrl_pts_1(0,1) = {0.65,0.35,0};
	ctrl_pts_1(0,2) = {1,0,0};
	ctrl_pts_1(1,0) = {0,-0.5,0};
	ctrl_pts_1(1,1) = {0.5,-0.6,0};
	ctrl_pts_1(1,2) = {1,-0.5,0};
	ctrl_pts_1(2,0) = {0,-1,0};
	ctrl_pts_1(2,1) = {0.65,-0.8,0};
	ctrl_pts_1(2,2) = {1,-1,0};

	math::BezierSurface bezier_surface_1(ctrl_pts_1);


	/*
	Array2D<math::Point> ctrl_pts_1(3,4);
	ctrl_pts_1(0,0) = {1,0,0};
	ctrl_pts_1(0,1) = {1.35,-0.35,0};
	ctrl_pts_1(0,2) = {1.6,0.1,0};
	ctrl_pts_1(0,3) = {2.0,0,0};

	ctrl_pts_1(1,0) = {1,-0.5,0};
	ctrl_pts_1(1,1) = {1.35,-0.75,0};
	ctrl_pts_1(1,2) = {1.7,-0.6,0};
	ctrl_pts_1(1,3) = {2,-0.5,0};

	ctrl_pts_1(2,0) = {1,-1,0};
	ctrl_pts_1(2,1) = {1.35,-1.2,0};
	ctrl_pts_1(2,2) = {1.7,-1.2,0};
	ctrl_pts_1(2,3) = {2,-1,0};
	math::BezierSurface bezier_surface_1(ctrl_pts_1);
	 */

	int Nx(10);
	int Ny(10);
	Array2D<math::Point> pts = bezier_surface_1.getDiscretization(10,10);

	std::ofstream stream= std::ofstream("BezierFace_1.table", std::ios::out);
	stream.precision(15);

	for (auto i=0;i<Nx;i++)
	{
		for (auto j=0;j<Ny;j++)
		{
			stream << pts(i,j).X() << " " << pts(i,j).Y() << "\n";
			stream << pts(i+1,j).X() << " " << pts(i+1,j).Y() << "\n";
			stream << pts(i+1,j+1).X() << " " << pts(i+1,j+1).Y() << "\n";
			stream << pts(i,j+1).X() << " " << pts(i,j+1).Y() << "\n";
			stream << pts(i,j).X() << " " << pts(i,j).Y() << "\n";
			stream << "\n";
		}
	}
	stream.close();

}
/*----------------------------------------------------------------------------*/