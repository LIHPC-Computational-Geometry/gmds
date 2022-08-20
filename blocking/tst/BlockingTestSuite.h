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
/*----------------------------------------------------------------------------*/