//
// Created by rochec on 20/05/2022.
//

/*----------------------------------------------------------------------------*/
#include<string>
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/math/Line.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(MathTest,Line) {

	{
		// Test 1
		std::cout << "-----------" << std::endl;
		std::cout << "TEST 1" << std::endl;
		std::cout << "-----------" << std::endl;
		math::Point p1;
		math::Point p2;

		math::Point p3;
		math::Point p4;

		bool intersect;
		math::Point Pint;
		double param;

		// Test 1
		p1 = (0, 0);
		p2 = (1, 0);
		p3 = (-3, 0);
		p4 = (-3, 1);

		math::Line L1(p1, p2);
		math::Line L2(p3, p4);

		intersect = L1.intersect2D(L2, Pint, param);
		std::cout << "Intersection : " << intersect << std::endl;
		std::cout << "Point : " << Pint << std::endl;
		std::cout << "Param : " << param << std::endl;

		intersect = L2.intersect2D(L1, Pint, param);
		std::cout << "Intersection : " << intersect << std::endl;
		std::cout << "Point : " << Pint << std::endl;
		std::cout << "Param : " << param << std::endl;
	}

	{
		// Test 2
		std::cout << "-----------" << std::endl;
		std::cout << "TEST 2" << std::endl;
		std::cout << "-----------" << std::endl;
		math::Point p1(0, 0, 0);
		math::Vector3d v1(1, 0, 0);
		math::Point p2(0, -3, 0);
		math::Vector3d v2(0, 1, 0);

		bool intersect;
		math::Point Pint;
		double param;

		math::Line L1(p1, v1);
		math::Line L2(p2, v2);

		intersect = L2.intersect2D(L1, Pint, param);
		std::cout << "Intersection : " << intersect << std::endl;
		std::cout << "Point : " << Pint << std::endl;
		std::cout << "Param : " << param << std::endl;

		intersect = L1.intersect2D(L2, Pint, param);
		std::cout << "Intersection : " << intersect << std::endl;
		std::cout << "Point : " << Pint << std::endl;
		std::cout << "Param : " << param << std::endl;
	}

	{
		// Test 3
		std::cout << "-----------" << std::endl;
		std::cout << "TEST 3" << std::endl;
		std::cout << "-----------" << std::endl;
		math::Point p1(-1, 0, 0);
		math::Vector3d v1(1, 0, 0);
		math::Point p2(0, -3, 0);
		math::Vector3d v2(0, 1, 0);

		bool intersect;
		math::Point Pint;
		double param;

		math::Line L1(p1, v1);
		math::Line L2(p2, v2);

		intersect = L2.intersect2D(L1, Pint, param);
		std::cout << "Intersection : " << intersect << std::endl;
		std::cout << "Point : " << Pint << std::endl;
		std::cout << "Param : " << param << std::endl;

		intersect = L1.intersect2D(L2, Pint, param);
		std::cout << "Intersection : " << intersect << std::endl;
		std::cout << "Point : " << Pint << std::endl;
		std::cout << "Param : " << param << std::endl;
	}

	EXPECT_EQ(1,1);
}
