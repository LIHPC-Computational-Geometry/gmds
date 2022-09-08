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