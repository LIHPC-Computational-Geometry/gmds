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
TEST(WriterDartsVTKTestSuite, screenshot_2d)
{
	gmds::blocking::Blocking bl;

	gmds::blocking::Point pt0(0.825, 0.825, 0);
	gmds::blocking::Point pt1(-0.4125, -0.4125, 0);
	gmds::blocking::Point pt2(4.175,  0.825, 0);
	gmds::blocking::Point pt3(5.825,  -0.4125, 0);
	gmds::blocking::Point pt4(5.825,  4.175, 0);
	gmds::blocking::Point pt5(0.825,  4.175, 0);
	gmds::blocking::Point pt6(-0.4125, 3.96875, 0);
	gmds::blocking::Point pt7(9.175,  -0.4125, 0);
	gmds::blocking::Point pt8(10.4125, -0.4125, 0);
	gmds::blocking::Point pt9(9.175,  4.175, 0);
	gmds::blocking::Point pt10(10.4125, 5, 0);
	gmds::blocking::Point pt11(8.175, 5, 0);
	gmds::blocking::Point pt12(5, 5, 0);
	gmds::blocking::Point pt13(5.825, 5.825, 0);
	gmds::blocking::Point pt14(4.175, 9.175, 0);
	gmds::blocking::Point pt15(5.825, 10.4125, 0);
	gmds::blocking::Point pt16(0.825, 9.175, 0);
	gmds::blocking::Point pt17(-0.4125, 10.4125, 0);
	gmds::blocking::Point pt18(9.175, 5.825, 0);
	gmds::blocking::Point pt19(9.175, 10.4125, 0);
	gmds::blocking::Point pt20(10.4125, 10.4125, 0);

	gmds::blocking::Dart_handle dh0  = bl.lcc()->make_quadrangle(pt0 , pt2 , pt12, pt5);
	gmds::blocking::Dart_handle dh1  = bl.lcc()->make_quadrangle(pt3 , pt7 , pt9 , pt4);
	gmds::blocking::Dart_handle dh2  = bl.lcc()->make_quadrangle(pt16, pt5 , pt12, pt14);
	gmds::blocking::Dart_handle dh3  = bl.lcc()->make_quadrangle(pt15, pt13, pt18, pt19);
	gmds::blocking::Dart_handle dh4  = bl.lcc()->make_quadrangle(pt1 , pt3 , pt2 , pt0);
	gmds::blocking::Dart_handle dh5  = bl.lcc()->make_quadrangle(pt2 , pt3 , pt4 , pt12);
	gmds::blocking::Dart_handle dh6  = bl.lcc()->make_quadrangle(pt1 , pt0 , pt5 , pt6);
	gmds::blocking::Dart_handle dh7  = bl.lcc()->make_quadrangle(pt7 , pt8 , pt10, pt9);
	gmds::blocking::Dart_handle dh8  = bl.lcc()->make_quadrangle(pt10, pt18, pt11, pt9);
	gmds::blocking::Dart_handle dh9  = bl.lcc()->make_quadrangle(pt4 , pt9 , pt11, pt12);
	gmds::blocking::Dart_handle dh10 = bl.lcc()->make_quadrangle(pt12, pt13, pt15, pt14);
	gmds::blocking::Dart_handle dh11 = bl.lcc()->make_quadrangle(pt16, pt14, pt15, pt17);
	gmds::blocking::Dart_handle dh12 = bl.lcc()->make_quadrangle(pt6 , pt5 , pt16, pt17);
	gmds::blocking::Dart_handle dh13 = bl.lcc()->make_quadrangle(pt12, pt11, pt18, pt13);
	gmds::blocking::Dart_handle dh14 = bl.lcc()->make_quadrangle(pt18, pt10, pt20, pt19);

	bl.lcc()->sew<2>(dh0, bl.lcc()->alpha(dh4, 1, 0, 1));
	bl.lcc()->sew<2>(bl.lcc()->alpha(dh0, 1), bl.lcc()->alpha(dh6, 0, 1));
	bl.lcc()->sew<2>(bl.lcc()->alpha(dh0, 1, 0, 1), bl.lcc()->alpha(dh2, 0, 1));
	bl.lcc()->sew<2>(bl.lcc()->alpha(dh0, 0, 1), bl.lcc()->alpha(dh5, 1));

	bl.lcc()->sew<2>(bl.lcc()->alpha(dh1, 1), bl.lcc()->alpha(dh5, 0, 1));
	bl.lcc()->sew<2>(bl.lcc()->alpha(dh1, 0, 1), bl.lcc()->alpha(dh7, 1));
	bl.lcc()->sew<2>(bl.lcc()->alpha(dh1, 1, 0, 1), dh9);

	bl.lcc()->sew<2>(bl.lcc()->alpha(dh2, 0, 1, 1), bl.lcc()->alpha(dh12, 0, 1));
	bl.lcc()->sew<2>(bl.lcc()->alpha(dh2, 0, 1, 0, 1), bl.lcc()->alpha(dh10, 1));
	bl.lcc()->sew<2>(bl.lcc()->alpha(dh2, 1), dh11);

	bl.lcc()->sew<2>(bl.lcc()->alpha(dh3, 0, 1), bl.lcc()->alpha(dh13, 1, 0, 1));
	bl.lcc()->sew<2>(bl.lcc()->alpha(dh3, 0, 1, 0, 1), bl.lcc()->alpha(dh14, 1));
	bl.lcc()->sew<2>(dh3, bl.lcc()->alpha(dh10, 0, 1, 0));

	bl.lcc()->sew<2>(bl.lcc()->alpha(dh4, 0, 1), bl.lcc()->alpha(dh5, 0));
	bl.lcc()->sew<2>(bl.lcc()->alpha(dh4, 1), dh6);

	bl.lcc()->sew<2>(bl.lcc()->alpha(dh5, 0, 1, 0, 1), bl.lcc()->alpha(dh9, 1));

	bl.lcc()->sew<2>(bl.lcc()->alpha(dh6, 1, 0, 1), dh12);

	bl.lcc()->sew<2>(bl.lcc()->alpha(dh7, 1, 0, 1), bl.lcc()->alpha(dh8, 1, 0));

	bl.lcc()->sew<2>(bl.lcc()->alpha(dh8, 1, 0, 1), bl.lcc()->alpha(dh9, 0, 1));
	bl.lcc()->sew<2>(bl.lcc()->alpha(dh8, 0, 1), bl.lcc()->alpha(dh13, 0, 1, 0));
	bl.lcc()->sew<2>(dh8, bl.lcc()->alpha(dh14, 0));

	bl.lcc()->sew<2>(bl.lcc()->alpha(dh9, 1, 0, 1), dh13);

	bl.lcc()->sew<2>(bl.lcc()->alpha(dh10, 1, 0, 1), bl.lcc()->alpha(dh11, 0, 1));

	bl.lcc()->sew<2>(bl.lcc()->alpha(dh11, 1), bl.lcc()->alpha(dh12, 0, 1, 0, 1));

	gmds::blocking::WriterDartsVTK writer;
	writer.setBl(&bl);
	gmds::blocking::WriterDartsVTK::STATUS st = writer.execute("darts_screenshot_2d.vtk");
	ASSERT_EQ(gmds::blocking::WriterDartsVTK::SUCCESS, st);
}
/*----------------------------------------------------------------------------*/
