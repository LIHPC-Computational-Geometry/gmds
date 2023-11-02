/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <vector>
/*----------------------------------------------------------------------------*/
//#include <gmds/blocking/Blocking.h>
/*----------------------------------------------------------------------------*/
#if WITH_CGAL_QT5
#include <CGAL/draw_linear_cell_complex.h>
#endif
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
#if WITH_CGAL_QT5
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
#endif
/*----------------------------------------------------------------------------*/
TEST(CGALTestSuite, using_lcc_gmaps_attributes)
{
	typedef CGAL::Linear_cell_complex_for_generalized_map<3> LCC_3;
	typedef LCC_3::size_type size_type;
	typedef LCC_3::Dart_handle Dart_handle;
	typedef LCC_3::Point Point;

	LCC_3 lcc;
	Dart_handle dh1=
	   lcc.make_hexahedron(Point(0,0,0), Point(5,0,0),
	                       Point(5,5,0), Point(0,5,0),
	                       Point(0,5,4), Point(0,0,4),
	                       Point(5,0,4), Point(5,5,4));

	typedef LCC_3::Attribute_handle<3>::type region_attribute;

//	for (LCC_3::Attribute_range<3>::type::iterator it = lcc.attributes<3>().begin(), itend = lcc.attributes<2>().end();
//             it != itend; ++it) {
//		id = lcc.info_of_attribute<2>(it)[0];
//	}

//	LCC_3::CellAttribute ca;

//	lcc.create_attribute<3>(int);

//	ASSERT_EQ(lcc.number_of_attributes<0>(), 1);
//	ASSERT_EQ(lcc.number_of_attributes<3>(), 8);
//	lcc.is_attribute_used<3>();

	// 1) Create all 2-attributes and associated them to darts.
//	for (LCC_3::Dart_range::iterator
//			  it=lcc.darts().begin(), itend=lcc.darts().end();
//		  it!=itend; ++it)
//	{
//		if ( lcc.attribute<2>(it)==nullptr ) {
//			lcc.set_attribute<2>(it, lcc.create_attribute<2>());
//		}
//	}

//	// Set the "color" of all vertices of the first cube to 1.
//	for (LCC_3::One_dart_per_incident_cell_range<0, 3>::iterator
//			  it=lcc.one_dart_per_incident_cell<0,3>(dh1).begin(),
//			  itend=lcc.one_dart_per_incident_cell<0,3>(dh1).end(); it!=itend; ++it) {
////		lcc.info<0>(it)=1;
//	}
}
/*----------------------------------------------------------------------------*/
TEST(CGALTestSuite, using_lcc_gmaps_marks)
{
	typedef CGAL::Linear_cell_complex_for_generalized_map<3> LCC_3;
	typedef LCC_3::size_type size_type;
	typedef LCC_3::Dart_handle Dart_handle;
	typedef LCC_3::Point Point;

	LCC_3 lcc;
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

	LCC_3::size_type mark = lcc.get_new_mark();
	lcc.mark(dh1, mark);

	ASSERT_EQ(lcc.number_of_used_marks(), 1);
	ASSERT_EQ(lcc.number_of_marked_darts(mark), 1);


}
/*----------------------------------------------------------------------------*/
TEST(CGALTestSuite, using_links)
{
	typedef CGAL::Generalized_map<3> GMap_3;
	typedef GMap_3::Dart_handle Dart_handle;

	GMap_3 gm;
	Dart_handle d1 = gm.create_dart();
	Dart_handle d2 = gm.create_dart();

	ASSERT_EQ(d1, gm.alpha(d1, 0));
	ASSERT_EQ(d2, gm.alpha(d2, 0));
	gm.link_alpha<0>(d1, d2);
	ASSERT_EQ(d2, gm.alpha(d1, 0));
	ASSERT_EQ(d1, gm.alpha(d2, 0));
}
/*----------------------------------------------------------------------------*/
TEST(CGALTestSuite, vertex_position)
{
	typedef CGAL::Linear_cell_complex_for_generalized_map<3> LCC_3;
	typedef LCC_3::size_type size_type;
	typedef LCC_3::Dart_handle Dart_handle;
	typedef LCC_3::Point Point;

	LCC_3 lcc;
	LCC_3::Vertex_attribute_handle v0 = lcc.create_vertex_attribute( LCC_3::Point (0.5, 0, 0));
	LCC_3::Vertex_attribute_handle v1 = lcc.create_vertex_attribute( LCC_3::Point (1, 0, 0));
	LCC_3::Vertex_attribute_handle v2 = lcc.create_vertex_attribute( LCC_3::Point (2, 0, 0));
	LCC_3::Vertex_attribute_handle v3 = lcc.create_vertex_attribute( LCC_3::Point (3, 0, 0));

	Dart_handle d0 = lcc.create_dart(v0);
	Dart_handle d1 = lcc.create_dart(v1);
	Dart_handle d2 = lcc.create_dart(v2);
	Dart_handle d3 = lcc.create_dart(v3);

	ASSERT_TRUE(lcc.is_valid());
	ASSERT_EQ(4, lcc.number_of_vertex_attributes());

//	std::cout<<lcc.vertex_attribute(d0)->point()<<std::endl;
	//	std::cout<<lcc.vertex_attribute(d1)->point()<<std::endl;
	//std::cout<<lcc.vertex_attribute(d2)->point()<<std::endl;
	//std::cout<<lcc.vertex_attribute(d3)->point()<<std::endl;

	std::set<Dart_handle> darts;
	darts.insert(d0);
	darts.insert(d1);

	lcc.link_alpha<1>(d0, d1);

	ASSERT_TRUE(darts.find(d0) != darts.end());
	ASSERT_TRUE(darts.find(d1) != darts.end());
	ASSERT_NE(d0, d1);

	//std::cout<<lcc.vertex_attribute(d0)->point()<<std::endl;
	//std::cout<<lcc.vertex_attribute(d1)->point()<<std::endl;
	//std::cout<<lcc.vertex_attribute(d2)->point()<<std::endl;
	//std::cout<<lcc.vertex_attribute(d3)->point()<<std::endl;
	ASSERT_EQ(lcc.vertex_attribute(d0)->point(), lcc.vertex_attribute(d1)->point());
	ASSERT_EQ(3, lcc.number_of_vertex_attributes());


	//std::cout<<lcc.vertex_attribute(d0)->point()<<std::endl;
	//std::cout<<lcc.vertex_attribute(d1)->point()<<std::endl;

	lcc.link_alpha<2>(d0, d2);
	ASSERT_TRUE(lcc.is_valid());
	ASSERT_EQ(lcc.vertex_attribute(d0)->point(), lcc.vertex_attribute(d2)->point());

//	lcc.correct_invalid_attributes();
	ASSERT_EQ(2, lcc.number_of_vertex_attributes());

	ASSERT_TRUE(lcc.is_valid());
	lcc.link_alpha<3>(d0, d3);
	ASSERT_FALSE(lcc.is_valid());

	ASSERT_EQ(lcc.vertex_attribute(d0)->point(), lcc.vertex_attribute(d3)->point());
	ASSERT_EQ(1, lcc.number_of_vertex_attributes());

	//std::cout<< "automatic " << lcc.are_attributes_automatically_managed()<<std::endl;
	//std::cout<<"attribute "<<lcc.is_attribute_used<0>(v0)<<" "<<lcc.is_valid_attribute<0>(v0)<<std::endl;
	//std::cout<<"attribute "<<lcc.is_attribute_used<0>(v1)<<" "<<lcc.is_valid_attribute<0>(v1)<<std::endl;
	//std::cout<<"attribute "<<lcc.is_attribute_used<0>(v2)<<" "<<lcc.is_valid_attribute<0>(v2)<<std::endl;
	//std::cout<<"attribute "<<lcc.is_attribute_used<0>(v3)<<" "<<lcc.is_valid_attribute<0>(v3)<<std::endl;
	lcc.correct_invalid_attributes();
	//std::cout<<"attribute "<<lcc.is_attribute_used<0>(v0)<<" "<<lcc.is_valid_attribute<0>(v0)<<std::endl;
	//std::cout<<"attribute "<<lcc.is_attribute_used<0>(v1)<<" "<<lcc.is_valid_attribute<0>(v1)<<std::endl;
	//std::cout<<"attribute "<<lcc.is_attribute_used<0>(v2)<<" "<<lcc.is_valid_attribute<0>(v2)<<std::endl;
	//std::cout<<"attribute "<<lcc.is_attribute_used<0>(v3)<<" "<<lcc.is_valid_attribute<0>(v3)<<std::endl;

	//std::cout<<"nbref "<<v0->get_nb_refs()<<" "<<v1->get_nb_refs()<<" "<<v2->get_nb_refs()<<" "<<v3->get_nb_refs()<<std::endl;

	lcc.unlink_alpha(d0, 3);
	lcc.correct_invalid_attributes();
	ASSERT_TRUE(lcc.is_valid());
	ASSERT_EQ(2, lcc.number_of_vertex_attributes());
	ASSERT_EQ(lcc.vertex_attribute(d0)->point(), lcc.vertex_attribute(d3)->point());
	//std::cout<<lcc.vertex_attribute(d0)->point()<<" "<<lcc.vertex_attribute(d3)->point()<<std::endl;
	v0->point() = LCC_3::Point (7,7,7);
	//std::cout<<lcc.vertex_attribute(d0)->point()<<" "<<lcc.vertex_attribute(d3)->point()<<std::endl;
}
/*----------------------------------------------------------------------------*/
TEST(CGALTestSuite, vertex_attributes)
{
	typedef CGAL::Linear_cell_complex_for_generalized_map<3> LCC_3;
	typedef LCC_3::size_type size_type;
	typedef LCC_3::Dart_handle Dart_handle;
	typedef LCC_3::Point Point;

	LCC_3 lcc;
	LCC_3::Vertex_attribute_handle v0 = lcc.create_vertex_attribute(LCC_3::Point(0.5, 0, 0));
	LCC_3::Vertex_attribute_handle v1 = lcc.create_vertex_attribute(LCC_3::Point(1, 0, 0));

	Dart_handle d0 = lcc.create_dart(v0);
	Dart_handle d1 = lcc.create_dart(v1);

	//std::cout<<"attribute "<<lcc.is_attribute_used<0>(v0)<<" "<<lcc.is_valid_attribute<0>(v0)<<std::endl;
	//std::cout<<"attribute "<<lcc.is_attribute_used<0>(v1)<<" "<<lcc.is_valid_attribute<0>(v1)<<std::endl;

	lcc.link_alpha<1>(d0, d1);
	//std::cout<<"attribute "<<lcc.is_attribute_used<0>(v0)<<" "<<lcc.is_valid_attribute<0>(v0)<<std::endl;
	//std::cout<<"attribute "<<lcc.is_attribute_used<0>(v1)<<" "<<lcc.is_valid_attribute<0>(v1)<<std::endl;
	ASSERT_EQ(1, lcc.number_of_vertex_attributes());

	lcc.set_vertex_attribute_of_dart(d1, v1);
	//std::cout<<"attribute "<<lcc.is_attribute_used<0>(v0)<<" "<<lcc.is_valid_attribute<0>(v0)<<std::endl;
	//std::cout<<"attribute "<<lcc.is_attribute_used<0>(v1)<<" "<<lcc.is_valid_attribute<0>(v1)<<std::endl;


	ASSERT_EQ(1, lcc.number_of_vertex_attributes());
}
/*----------------------------------------------------------------------------*/