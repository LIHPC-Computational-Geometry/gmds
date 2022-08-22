/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
//#include <filesystem>
//#include <iostream>
//#include <vector>
/*----------------------------------------------------------------------------*/
#include <gmds/blocking/InsertSheet.h>

#include <gmds/blocking/Blocking.h>
/*----------------------------------------------------------------------------*/
TEST(InsertSheetTestSuite, dummytest)
{
	ASSERT_EQ(0, 0);
}
/*----------------------------------------------------------------------------*/
TEST(InsertSheetTestSuite, instanciate)
{
	gmds::blocking::Blocking bl;
	bl.createGrid2d();

	gmds::blocking::InsertSheet is;
	is.setBl(&bl);

	gmds::blocking::InsertSheet::STATUS status = is.execute();
	ASSERT_EQ(status, gmds::blocking::InsertSheet::NOT_YET_IMPLEMENTED);
}
/*----------------------------------------------------------------------------*/