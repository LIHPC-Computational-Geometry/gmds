/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <vector>
/*----------------------------------------------------------------------------*/
//#include <gmds/blocking/Blocking.h>
/*----------------------------------------------------------------------------*/
//#include <CGAL/draw_linear_cell_complex.h>
#include <CGAL/Generalized_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
/*----------------------------------------------------------------------------*/
TEST(CGALTestSuite, using_cgal)
{
	// This is taken from CGAL Generalized_map/gmap_4_simple_example.cpp

	typedef CGAL::Generalized_map<4> GMap_4;
	typedef GMap_4::Dart_handle Dart_handle;

	GMap_4 gm;
	Dart_handle d1 = gm.make_combinatorial_tetrahedron();
	Dart_handle d2 = gm.make_combinatorial_tetrahedron();
	ASSERT_TRUE(gm.is_valid());
	gm.sew<4>(d1,d2);
//	gm.display_characteristics(std::cout);
	ASSERT_TRUE(gm.is_valid());
}
/*----------------------------------------------------------------------------*/
TEST(CGALTestSuite, using_lcc_gmaps)
{
	// This is taken from CGAL Linear_cell_complex/linear_cell_complex_4.cpp

	typedef CGAL::Linear_cell_complex_for_generalized_map<4,5> LCC_4;
	typedef LCC_4::Dart_handle                                 Dart_handle;
	typedef LCC_4::Point                                       Point;
	typedef LCC_4::Vector                                      Vector;
	typedef LCC_4::FT                                          FT;

	LCC_4 lcc;
	// Create two tetrahedra.
	FT p1[5]={ 0, 0, 0, 0, 0}; std::vector<FT> v1(p1,p1+5);
	FT p2[5]={ 0, 2, 0, 0, 0}; std::vector<FT> v2(p2,p2+5);
	FT p3[5]={ 0, 1, 2, 0, 0}; std::vector<FT> v3(p3,p3+5);
	FT p4[5]={ 2, 1, 0, 0, 0}; std::vector<FT> v4(p4,p4+5);
	FT p5[5]={-1, 0, 0, 0, 0}; std::vector<FT> v5(p5,p5+5);
	FT p6[5]={-1, 2, 0, 0, 0}; std::vector<FT> v6(p6,p6+5);
	FT p7[5]={-1, 1, 2, 0, 0}; std::vector<FT> v7(p7,p7+5);
	FT p8[5]={-3, 1, 2, 0, 0}; std::vector<FT> v8(p8,p8+5);
	Dart_handle d1 = lcc.make_tetrahedron(Point(5, v1.begin(), v1.end()),
													  Point(5, v2.begin(), v2.end()),
													  Point(5, v3.begin(), v3.end()),
													  Point(5, v4.begin(), v4.end()));
	Dart_handle d2 = lcc.make_tetrahedron(Point(5, v5.begin(), v5.end()),
													  Point(5, v6.begin(), v6.end()),
													  Point(5, v7.begin(), v7.end()),
													  Point(5, v8.begin(), v8.end()));

	ASSERT_TRUE(lcc.is_valid());

	lcc.sew<4>(d1,d2);
//	lcc.display_characteristics(std::cout);

	// Add one vertex on the middle of the edge containing dart d1.
	Dart_handle d3 = lcc.insert_barycenter_in_cell<1>(d1);
	ASSERT_TRUE(lcc.is_valid());
//	lcc.display_characteristics(std::cout);

	// Add one edge to cut the face containing dart d3 in two.
	Dart_handle d4 = lcc.insert_cell_1_in_cell_2(d3, lcc.alpha(d1, 1, 0, 1));
	ASSERT_TRUE(lcc.is_valid());
//	lcc.display_characteristics(std::cout);

	// We use removal operations to get back to the initial configuration.
	lcc.remove_cell<1>(d4);
	ASSERT_TRUE(lcc.is_valid());
	lcc.remove_cell<0>(d3);
	ASSERT_TRUE(lcc.is_valid());
	lcc.unsew<4>(d1);
	ASSERT_TRUE(lcc.is_valid());
//	lcc.display_characteristics(std::cout);
}
/*----------------------------------------------------------------------------*/
TEST(CGALTestSuite, DISABLED_draw)
{
	// This is taken from CGAL Linear_cell_complex/draw_linear_cell_complex.cpp

	typedef CGAL::Linear_cell_complex_for_generalized_map<3> LCC;
	typedef LCC::Dart_handle Dart_handle;
	typedef LCC::Point Point;

	LCC lcc;
	Dart_handle dh1=
		lcc.make_hexahedron(Point(0,0,0), Point(5,0,0),
								  Point(5,5,0), Point(0,5,0),
								  Point(0,5,4), Point(0,0,4),
								  Point(5,0,4), Point(5,5,4));
	Dart_handle dh2=
		lcc.make_hexahedron(Point(5,0,0), Point(10,0,0),
								  Point(10,5,0), Point(5,5,0),
								  Point(5,5,4), Point(5,0,4),
								  Point(10,0,4), Point(10,5,4));
//	lcc.sew<3>(lcc.beta(dh1, 1, 1, 2), lcc.beta(dh2, 2));
	ASSERT_TRUE(lcc.is_valid());
//	CGAL::draw(lcc);
}
/*----------------------------------------------------------------------------*/