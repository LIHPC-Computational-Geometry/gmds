/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <filesystem>
//#include <iostream>
//#include <vector>
/*----------------------------------------------------------------------------*/
#include <gmds/blocking/Blocking.h>
/*----------------------------------------------------------------------------*/
TEST(BlockingTestSuite, dummytest)
{
	ASSERT_EQ(0, 0);
}
/*----------------------------------------------------------------------------*/
TEST(BlockingTestSuite, createGrid)
{
	gmds::blocking::Blocking bl2d;
	bl2d.createGrid2d(gmds::math::Point(0,0,0), gmds::math::Point(1,1,1), 2,3);
	ASSERT_EQ(bl2d.nbVertices(), 12);
	ASSERT_EQ(bl2d.nbEdges(), 17);
	ASSERT_EQ(bl2d.nbFaces(), 6);
//	ASSERT_EQ(bl2d.nbBlocks(), 0);

	gmds::blocking::Blocking bl3d;
	bl3d.createGrid3d(gmds::math::Point(0,0,0), gmds::math::Point(1,1,1), 2,3,4);
	ASSERT_EQ(bl3d.nbVertices(), 60);
	ASSERT_EQ(bl3d.nbEdges(), 133);
	ASSERT_EQ(bl3d.nbFaces(), 98);
	ASSERT_EQ(bl3d.nbBlocks(), 24);
}
/*----------------------------------------------------------------------------*/
TEST(BlockingTestSuite, writeMoka)
{
	gmds::blocking::Blocking bl2d;
	bl2d.createGrid2d();
	ASSERT_EQ(bl2d.nbVertices(), 16);
	ASSERT_EQ(bl2d.nbEdges(), 24);
	ASSERT_EQ(bl2d.nbFaces(), 9);
	//	ASSERT_EQ(bl2d.nbBlocks(), 0);
	bl2d.writeMokaFile("grid2d.moka");

	gmds::blocking::Blocking bl3d;
	bl3d.createGrid3d();
	ASSERT_EQ(bl3d.nbVertices(), 64);
	ASSERT_EQ(bl3d.nbEdges(), 144);
	ASSERT_EQ(bl3d.nbFaces(), 108);
	ASSERT_EQ(bl3d.nbBlocks(), 27);
	bl3d.writeMokaFile("grid3d.moka");
//#if __cplusplus >= 201703L
#if __cpp_lib_filesystem
	std::filesystem::exists(std::filesystem::path ("grid3d.moka"));
#endif // C++17
}
/*----------------------------------------------------------------------------*/
TEST(BlockingTestSuite, writeVTK_2d)
{
	gmds::blocking::Blocking bl2d;
	bl2d.createGrid2d();
	bl2d.writeVTKFile("grid2d.vtk");
	ASSERT_EQ(bl2d.nbVertices(), 16);
}
/*----------------------------------------------------------------------------*/
TEST(BlockingTestSuite, writeVTK_3d)
{
	gmds::blocking::Blocking bl3d;
	bl3d.createGrid3d();
	bl3d.writeVTKFile("grid3d.vtk");
	ASSERT_EQ(bl3d.nbVertices(), 64);
}
/*----------------------------------------------------------------------------*/
TEST(BlockingTestSuite, unsew_and_vertices_move)
{
	typedef CGAL::Linear_cell_complex_for_generalized_map<3> LCC_3;
	typedef LCC_3::size_type size_type;
	typedef LCC_3::Dart_handle Dart_handle;
	typedef LCC_3::Point Point;

	gmds::blocking::Blocking bl;

	Dart_handle dh1=
		bl.lcc()->make_hexahedron(Point(0,0,0), Point(5,0,0),
										  Point(5,5,0), Point(0,5,0),
										  Point(0,5,4), Point(0,0,4),
										  Point(5,0,4), Point(5,5,4));
	Dart_handle dh2=
		bl.lcc()->make_hexahedron(Point(5,0,0), Point(10,0,0),
										  Point(10,5,0), Point(5,5,0),
										  Point(5,5,4), Point(5,0,4),
										  Point(10,0,4), Point(10,5,4));

	bl.writeVTKFile("unsew_a.vtk");
	bl.lcc()->sew<3>(bl.lcc()->alpha(dh1, 1,0,1,2), bl.lcc()->alpha(dh2,2));
	bl.writeVTKFile("unsew_b.vtk");
	bl.lcc()->unsew<3>(bl.lcc()->alpha(dh1, 1,0,1,2));
	bl.writeVTKFile("unsew_c.vtk");

}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/