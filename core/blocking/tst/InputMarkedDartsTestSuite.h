/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
//#include <filesystem>
//#include <iostream>
//#include <vector>
/*----------------------------------------------------------------------------*/
#include <gmds/blocking/InputMarkedDarts.h>
#include <gmds/blocking/Blocking.h>
/*----------------------------------------------------------------------------*/
TEST(InputMarkedDartsTestSuite, dummytest)
{
	ASSERT_EQ(0, 0);
}
/*----------------------------------------------------------------------------*/
TEST(InputMarkedDartsTestSuite, pillow3d)
{
	gmds::blocking::Blocking bl;
	bl.createGrid3d(gmds::math::Point(0,0,0), gmds::math::Point(2,2,2), 2,2,2);

	gmds::blocking::LCC_3::size_type mark = bl.lcc()->get_new_mark();
	gmds::blocking::InputMarkedDarts imd;
	int nbMarkedDarts = imd.pillow_mark_first_cells3d(bl.lcc(), mark, 5);
	bl.lcc()->free_mark(mark);

	ASSERT_EQ(160, nbMarkedDarts);
}
/*----------------------------------------------------------------------------*/
TEST(InputMarkedDartsTestSuite, insert3d)
{
	gmds::blocking::Blocking bl;
	bl.createGrid3d(gmds::math::Point(0,0,0), gmds::math::Point(2,2,2), 2,2,2);

	gmds::blocking::LCC_3::size_type mark = bl.lcc()->get_new_mark();
	gmds::blocking::InputMarkedDarts imd;
	int nbMarkedDarts = imd.insertsheet_mark_first_cells3d(bl.lcc(), mark, 5);
	bl.lcc()->free_mark(mark);

	ASSERT_EQ(40, nbMarkedDarts);
}
/*----------------------------------------------------------------------------*/
TEST(InputMarkedDartsTestSuite, insert3d_intersect_3x2x3)
{
	gmds::blocking::Blocking bl;
	bl.createGrid3d(gmds::math::Point(0,0,0), gmds::math::Point(4,4,4), 3,2,3);

	gmds::blocking::LCC_3::size_type mark = bl.lcc()->get_new_mark();
	gmds::blocking::InputMarkedDarts imd;
	int nbMarkedDarts = imd.insertsheet_mark_intersect_3d(bl.lcc(), mark, 3, 2, 3);
	bl.lcc()->free_mark(mark);

	ASSERT_EQ(96, nbMarkedDarts);
}
/*----------------------------------------------------------------------------*/
TEST(InputMarkedDartsTestSuite, insert3d_intersect_4x4x4)
{
	gmds::blocking::Blocking bl;
	bl.createGrid3d(gmds::math::Point(0,0,0), gmds::math::Point(4,4,4), 4,4,4);

	gmds::blocking::LCC_3::size_type mark = bl.lcc()->get_new_mark();
	gmds::blocking::InputMarkedDarts imd;
	int nbMarkedDarts = imd.insertsheet_mark_intersect_3d(bl.lcc(), mark, 4, 4, 4);
	bl.lcc()->free_mark(mark);

	ASSERT_EQ(320, nbMarkedDarts);
}
/*----------------------------------------------------------------------------*/