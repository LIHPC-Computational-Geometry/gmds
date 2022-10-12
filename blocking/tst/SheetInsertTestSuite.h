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

#include <unit_test_config.h>
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

	gmds::blocking::WriterDartsVTK writer_tmp;
	writer_tmp.setBl(&bl);
	gmds::blocking::WriterDartsVTK::STATUS st_tmp = writer_tmp.execute("darts_tmp_pillow3d.vtk", mark);

	gmds::blocking::SheetInsert::STATUS status = is.pillow(mark);
	is.lcc()->free_mark(mark);

	ASSERT_EQ(status, gmds::blocking::SheetInsert::SUCCESS);

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
	bl.createGrid3d(gmds::math::Point(0,0,0), gmds::math::Point(2,4,4), 2,4,4);

	gmds::blocking::SheetInsert is;
	is.setBl(&bl);

	gmds::blocking::LCC_3::size_type mark = is.lcc()->get_new_mark();
	gmds::blocking::InputMarkedDarts imd;
	imd.insertsheet_mark_first_cells3d(bl.lcc(), mark, 16);

	gmds::blocking::WriterDartsVTK writer_tmp;
	writer_tmp.setBl(&bl);
	gmds::blocking::WriterDartsVTK::STATUS st_tmp = writer_tmp.execute("darts_tmp_insertsheet3d_firstcells.vtk", mark);

	gmds::blocking::SheetInsert::STATUS status = is.execute(mark);
	is.lcc()->free_mark(mark);

	ASSERT_EQ(status, gmds::blocking::SheetInsert::SUCCESS);

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

	bl.writeVTKFile("grid3d_insertsheet3d_autointersect.vtk");

	gmds::blocking::SheetInsert is;
	is.setBl(&bl);

	gmds::blocking::LCC_3::size_type mark = is.lcc()->get_new_mark();
	gmds::blocking::InputMarkedDarts imd;
	imd.insertsheet_mark_intersect_3d(bl.lcc(), mark, 4, 4, 4);

	gmds::blocking::WriterDartsVTK writer_tmp;
	writer_tmp.setBl(&bl);
	gmds::blocking::WriterDartsVTK::STATUS st_tmp = writer_tmp.execute("darts_tmp_insertsheet3d_autointersect.vtk", mark);

	gmds::blocking::SheetInsert::STATUS status = is.execute(mark);
	is.lcc()->free_mark(mark);

	ASSERT_EQ(status, gmds::blocking::SheetInsert::SUCCESS);

	bl.writeVTKFile("insertsheet3d_autointersect.vtk");
	bl.writeMokaFile("insertsheet3d_autointersect.moka");

	gmds::blocking::WriterDartsVTK writer;
	writer.setBl(&bl);
	gmds::blocking::WriterDartsVTK::STATUS st = writer.execute("darts_insertsheet3d_autointersect.vtk");
	ASSERT_EQ(gmds::blocking::WriterDartsVTK::SUCCESS, st);
}
/*----------------------------------------------------------------------------*/
TEST(SheetInsertTestSuite, DISABLED_insert3d_unstruct)
{
	gmds::blocking::Blocking bl;
	bl.createGrid3d(gmds::math::Point(0,0,0), gmds::math::Point(3,3,3), 1,1,1);

	gmds::blocking::SheetInsert is;
	is.setBl(&bl);

	gmds::blocking::LCC_3::size_type mark = is.lcc()->get_new_mark();
	gmds::blocking::InputMarkedDarts imd;
	imd.pillow_mark_first_cells3d(bl.lcc(), mark, 1);

	gmds::blocking::SheetInsert::STATUS status = is.pillow(mark);
	is.lcc()->free_mark(mark);

	ASSERT_EQ(status, gmds::blocking::SheetInsert::SUCCESS);

	bl.writeVTKFile("pillow3d_unstruct.vtk");
	bl.writeMokaFile("pillow3d_unstruct.moka");

	gmds::blocking::WriterDartsVTK writer;
	writer.setBl(&bl);
	gmds::blocking::WriterDartsVTK::STATUS st = writer.execute("darts_insertsheet3d_firstcells.vtk");
	ASSERT_EQ(gmds::blocking::WriterDartsVTK::SUCCESS, st);

	gmds::blocking::LCC_3::size_type mark_insert = is.lcc()->get_new_mark();
	imd.insertsheet_mark_first_cells3d(bl.lcc(), mark_insert, 3);
	is.execute(mark_insert);

	bl.writeVTKFile("insertsheet3d_unstruct.vtk");
	bl.writeMokaFile("insertsheet3d_unstruct.moka");
}
/*----------------------------------------------------------------------------*/
TEST(SheetInsertTestSuite, pillow_ogrid_3d)
{
	gmds::blocking::Blocking bl;

	std::string dir(TEST_SAMPLES_DIR);
	bl.readVTKFile(dir + "/ogrid.vtk");
	ASSERT_EQ(bl.nbBlocks(), 7);

	gmds::blocking::SheetInsert is;
	is.setBl(&bl);

	gmds::blocking::LCC_3::size_type mark = is.lcc()->get_new_mark();
	gmds::blocking::InputMarkedDarts imd;
	imd.pillow_mark_first_cells3d(bl.lcc(), mark, 3);

	gmds::blocking::WriterDartsVTK writer_tmp;
	writer_tmp.setBl(&bl);
	gmds::blocking::WriterDartsVTK::STATUS st_tmp = writer_tmp.execute("darts_tmp_pillow3d_ogrid.vtk", mark);

	gmds::blocking::SheetInsert::STATUS status_insert = is.pillow(mark);
	is.lcc()->free_mark(mark);

	ASSERT_EQ(19, bl.nbBlocks());
	ASSERT_EQ(status_insert, gmds::blocking::SheetInsert::SUCCESS);

	bl.writeVTKFile("pillow3d_ogrid.vtk");
	bl.writeMokaFile("pillow3d_ogrid.moka");

	gmds::blocking::WriterDartsVTK writer;
	writer.setBl(&bl);
	gmds::blocking::WriterDartsVTK::STATUS status_write = writer.execute("darts_pillow3d_ogrid.vtk");
	ASSERT_EQ(gmds::blocking::WriterDartsVTK::SUCCESS, status_write);
}
/*----------------------------------------------------------------------------*/
TEST(SheetInsertTestSuite, insertsheet_ogrid_3d)
{
	gmds::blocking::Blocking bl;

	std::string dir(TEST_SAMPLES_DIR);
	bl.readVTKFile(dir + "/ogrid.vtk");
	ASSERT_EQ(bl.nbBlocks(), 7);

	gmds::blocking::SheetInsert is;
	is.setBl(&bl);

	gmds::blocking::LCC_3::size_type mark = is.lcc()->get_new_mark();
	gmds::blocking::InputMarkedDarts imd;
	imd.insertsheet_mark_first_cells3d(bl.lcc(), mark, 3);

	gmds::blocking::WriterDartsVTK writer_tmp;
	writer_tmp.setBl(&bl);
	gmds::blocking::WriterDartsVTK::STATUS st_tmp = writer_tmp.execute("darts_tmp_insertsheet3d_ogrid.vtk", mark);

	gmds::blocking::SheetInsert::STATUS status_insert = is.execute(mark);
	is.lcc()->free_mark(mark);

	ASSERT_EQ(17, bl.nbBlocks());

	ASSERT_EQ(status_insert, gmds::blocking::SheetInsert::SUCCESS);

	bl.writeVTKFile("insertsheet3d_ogrid.vtk");
	bl.writeMokaFile("insertsheet3d_ogrid.moka");

	gmds::blocking::WriterDartsVTK writer;
	writer.setBl(&bl);
	gmds::blocking::WriterDartsVTK::STATUS status_write = writer.execute("darts_insertsheet3d_ogrid.vtk");
	ASSERT_EQ(gmds::blocking::WriterDartsVTK::SUCCESS, status_write);
}
/*----------------------------------------------------------------------------*/
TEST(SheetInsertTestSuite, insert3d_screenshots)
{
	gmds::blocking::Blocking bl;
	bl.createGrid3d(gmds::math::Point(0,0,0), gmds::math::Point(1,1,1), 1,1,2);

	gmds::blocking::SheetInsert is;
	is.setBl(&bl);

	bl.writeVTKFile("grid3d_insert3d_screenshot.vtk");

	gmds::blocking::LCC_3::size_type mark = is.lcc()->get_new_mark();
	gmds::blocking::InputMarkedDarts imd;
	imd.insertsheet_mark_first_cells3d(bl.lcc(), mark, 1);

	gmds::blocking::WriterDartsVTK writer_tmp;
	writer_tmp.setBl(&bl);
	gmds::blocking::WriterDartsVTK::STATUS st_tmp = writer_tmp.execute("darts_tmp_insertsheet3d_screenshot.vtk", mark);

	gmds::blocking::SheetInsert::STATUS status = is.execute(mark);
	is.lcc()->free_mark(mark);

	ASSERT_EQ(status, gmds::blocking::SheetInsert::SUCCESS);

	bl.writeVTKFile("insertsheet3d_screenshot.vtk");
	bl.writeMokaFile("insertsheet3d_screenshot.moka");

	gmds::blocking::WriterDartsVTK writer;
	writer.setBl(&bl);
	gmds::blocking::WriterDartsVTK::STATUS st = writer.execute("darts_insertsheet3d_screenshot.vtk");
	ASSERT_EQ(gmds::blocking::WriterDartsVTK::SUCCESS, st);
}
/*----------------------------------------------------------------------------*/