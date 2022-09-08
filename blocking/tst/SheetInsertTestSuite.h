/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
//#include <filesystem>
//#include <iostream>
//#include <vector>
/*----------------------------------------------------------------------------*/
#include <gmds/blocking/SheetInsert.h>

#include <gmds/blocking/Blocking.h>
/*----------------------------------------------------------------------------*/
TEST(SheetInsertTestSuite, dummytest)
{
	ASSERT_EQ(0, 0);
}
/*----------------------------------------------------------------------------*/
TEST(SheetInsertTestSuite, instanciate)
{
	gmds::blocking::Blocking bl;
	bl.createGrid2d();

	gmds::blocking::SheetInsert is;
	is.setBl(&bl);

	gmds::blocking::SheetInsert::STATUS status = is.execute();
	ASSERT_EQ(status, gmds::blocking::SheetInsert::NOT_YET_IMPLEMENTED);
}
/*----------------------------------------------------------------------------*/
TEST(SheetInsertTestSuite, pillow3d)
{
	gmds::blocking::Blocking bl;
//	bl.createGrid3d();
	bl.createGrid3d(gmds::math::Point(0,0,0), gmds::math::Point(2,2,2), 2,2,2);

	gmds::blocking::SheetInsert is;
	is.setBl(&bl);

	gmds::blocking::SheetInsert::STATUS status = is.pillow();
	ASSERT_EQ(status, gmds::blocking::SheetInsert::NOT_YET_IMPLEMENTED);

	bl.writeVTKFile("pillow3d.vtk");
	bl.writeMokaFile("pillow3d.moka");
}
/*----------------------------------------------------------------------------*/
TEST(SheetInsertTestSuite, insert3d)
{
	gmds::blocking::Blocking bl;
	//	bl.createGrid3d();
	bl.createGrid3d(gmds::math::Point(0,0,0), gmds::math::Point(2,2,2), 2,2,2);

	gmds::blocking::SheetInsert is;
	is.setBl(&bl);

	gmds::blocking::SheetInsert::STATUS status = is.execute();
	ASSERT_EQ(status, gmds::blocking::SheetInsert::NOT_YET_IMPLEMENTED);

	bl.writeVTKFile("insertsheet3d.vtk");
	bl.writeMokaFile("insertsheet3d.moka");
}
/*----------------------------------------------------------------------------*/