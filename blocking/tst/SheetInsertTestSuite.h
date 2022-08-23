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