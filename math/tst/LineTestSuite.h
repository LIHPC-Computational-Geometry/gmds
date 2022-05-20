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

	math::Point p1;
	math::Point p2;

	math::Point p3;
	math::Point p4;

	math::Vector3d v1;
	math::Vector3d v2;

	bool intersect;
	math::Point Pint;
	double param;

	// Test 1
	p1 = (0,0);
	p2 = (1,0);
	p3 = (-3,0);
	p4 = (-3, 1);

	math::Line L1(p1, p2);
	math::Line L2(p3, p4);

	intersect = L1.intersect2D(L2, Pint, param);
	std::cout << "Intersection : " << intersect << std::endl;
	std::cout << "Point : " << Pint << std::endl;
	std::cout << "Param : " << param << std::endl;

	EXPECT_EQ(1,1);
}
