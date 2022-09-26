/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <filesystem>
//#include <iostream>
/*----------------------------------------------------------------------------*/
#include <gmds/blocking/WriterDartsVTK.h>
/*----------------------------------------------------------------------------*/
TEST(WriterDartsVTKTestSuite, dummytest)
{
	ASSERT_EQ(0, 0);
}
/*----------------------------------------------------------------------------*/
TEST(WriterDartsVTKTestSuite, grid2d)
{
	gmds::blocking::Blocking bl;
	bl.createGrid2d(gmds::math::Point(0,0,0), gmds::math::Point(2,2, 0), 2,3);

	gmds::blocking::WriterDartsVTK writer;
	writer.setBl(&bl);
	gmds::blocking::WriterDartsVTK::STATUS st = writer.execute("darts_grid2d.vtk");
	ASSERT_EQ(gmds::blocking::WriterDartsVTK::SUCCESS, st);
}
/*----------------------------------------------------------------------------*/
TEST(WriterDartsVTKTestSuite, grid3d)
{
	gmds::blocking::Blocking bl;
	bl.createGrid3d(gmds::math::Point(0,0,0), gmds::math::Point(2,3,4), 2,3,4);

	gmds::blocking::WriterDartsVTK writer;
	writer.setBl(&bl);
	gmds::blocking::WriterDartsVTK::STATUS st = writer.execute("darts_grid3d.vtk");
	ASSERT_EQ(gmds::blocking::WriterDartsVTK::SUCCESS, st);
}
/*----------------------------------------------------------------------------*/