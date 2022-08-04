/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
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

	gmds::blocking::Blocking bl3d;
	bl3d.createGrid3d(gmds::math::Point(0,0,0), gmds::math::Point(1,1,1), 2,3,4);
	ASSERT_EQ(bl3d.nbVertices(), 60);
	ASSERT_EQ(bl3d.nbEdges(), 133);
	ASSERT_EQ(bl3d.nbFaces(), 98);
	ASSERT_EQ(bl3d.nbBlocks(), 24);
}
/*----------------------------------------------------------------------------*/
TEST(BlockingTestSuite, write)
{
	gmds::blocking::Blocking bl2d;
	bl2d.createGrid2d();
	bl2d.writeMokaFile("grid2d.moka");
	ASSERT_EQ(bl2d.nbVertices(), 16);

	gmds::blocking::Blocking bl3d;
	bl3d.createGrid3d();
	bl3d.writeMokaFile("grid3d.moka");
	ASSERT_EQ(bl3d.nbVertices(), 64);
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/