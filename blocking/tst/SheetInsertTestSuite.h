/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
//#include <filesystem>
//#include <iostream>
//#include <vector>
/*----------------------------------------------------------------------------*/
#include <gmds/blocking/Blocking.h>
#include <gmds/blocking/InputMarkedDarts.h>
#include <gmds/blocking/SheetInsert.h>
#include <gmds/blocking/WriterDartsVTK.h>
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

//	gmds::blocking::SheetInsert::STATUS status = is.execute();
//	ASSERT_EQ(status, gmds::blocking::SheetInsert::NOT_YET_IMPLEMENTED);
}
/*----------------------------------------------------------------------------*/
TEST(SheetInsertTestSuite, pillow3d)
{
	gmds::blocking::Blocking bl;
//	bl.createGrid3d();
	bl.createGrid3d(gmds::math::Point(0,0,0), gmds::math::Point(2,2,2), 2,2,2);

	gmds::blocking::SheetInsert is;
	is.setBl(&bl);

	gmds::blocking::LCC_3::size_type mark = is.lcc()->get_new_mark();
	gmds::blocking::InputMarkedDarts imd;
	imd.pillow_mark_first_cells3d(bl.lcc(), mark, 5);

	gmds::blocking::SheetInsert::STATUS status = is.pillow(mark);
	is.lcc()->free_mark(mark);

	ASSERT_EQ(status, gmds::blocking::SheetInsert::NOT_YET_IMPLEMENTED);

	bl.writeVTKFile("pillow3d.vtk");
	bl.writeMokaFile("pillow3d.moka");

	gmds::blocking::WriterDartsVTK writer;
	writer.setBl(&bl);
	gmds::blocking::WriterDartsVTK::STATUS st = writer.execute("darts_pillow3d.vtk");
	ASSERT_EQ(gmds::blocking::WriterDartsVTK::SUCCESS, st);
}
/*----------------------------------------------------------------------------*/
TEST(SheetInsertTestSuite, insert3d)
{
	gmds::blocking::Blocking bl;
	bl.createGrid3d(gmds::math::Point(0,0,0), gmds::math::Point(2,2,2), 2,2,2);

	gmds::blocking::SheetInsert is;
	is.setBl(&bl);

	gmds::blocking::LCC_3::size_type mark = is.lcc()->get_new_mark();
	gmds::blocking::InputMarkedDarts imd;
	imd.insertsheet_mark_first_cells3d(bl.lcc(), mark, 5);

	gmds::blocking::SheetInsert::STATUS status = is.execute(mark);
	is.lcc()->free_mark(mark);

	ASSERT_EQ(status, gmds::blocking::SheetInsert::NOT_YET_IMPLEMENTED);

	bl.writeVTKFile("insertsheet3d_firstcells.vtk");
	bl.writeMokaFile("insertsheet3d_firstcells.moka");

	gmds::blocking::WriterDartsVTK writer;
	writer.setBl(&bl);
	gmds::blocking::WriterDartsVTK::STATUS st = writer.execute("darts_insertsheet3d_firstcells.vtk");
	ASSERT_EQ(gmds::blocking::WriterDartsVTK::SUCCESS, st);
}
/*----------------------------------------------------------------------------*/
TEST(SheetInsertTestSuite, intersect3d)
{
	gmds::blocking::Blocking bl;
	bl.createGrid3d(gmds::math::Point(0,0,0), gmds::math::Point(4,4,4), 4,4,4);

	gmds::blocking::SheetInsert is;
	is.setBl(&bl);

	gmds::blocking::LCC_3::size_type mark = is.lcc()->get_new_mark();
	gmds::blocking::InputMarkedDarts imd;
	imd.insertsheet_mark_intersect_3d(bl.lcc(), mark, 4, 4, 4);

	gmds::blocking::SheetInsert::STATUS status = is.execute(mark);
	is.lcc()->free_mark(mark);

	ASSERT_EQ(status, gmds::blocking::SheetInsert::NOT_YET_IMPLEMENTED);

	bl.writeVTKFile("insertsheet3d_autointersect.vtk");
	bl.writeMokaFile("insertsheet3d_autointersect.moka");

	gmds::blocking::WriterDartsVTK writer;
	writer.setBl(&bl);
	gmds::blocking::WriterDartsVTK::STATUS st = writer.execute("darts_insertsheet3d_autointersect.vtk");
	ASSERT_EQ(gmds::blocking::WriterDartsVTK::SUCCESS, st);
}
/*----------------------------------------------------------------------------*/