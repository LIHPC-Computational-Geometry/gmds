/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <filesystem>
//#include <iostream>
//#include <vector>
/*----------------------------------------------------------------------------*/
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKWriter.h"
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
TEST(BlockingTestSuite, readVTK_3d)
{
	gmds::blocking::Blocking bl3d;
	bl3d.readVTKFile("grid3d.vtk");

//	ASSERT_EQ(bl3d.nbVertices(), 64);
//	ASSERT_EQ(bl3d.nbEdges(), 144);
//	ASSERT_EQ(bl3d.nbFaces(), 108);
	ASSERT_EQ(bl3d.nbBlocks(), 27);
}
/*----------------------------------------------------------------------------*/
TEST(BlockingTestSuite, createBlocks_3d)
{
	gmds::Mesh m3d(gmds::MeshModel(gmds::DIM3| gmds::R | gmds::N | gmds::R2N | gmds::N2R));

	gmds::Node n0 = m3d.newNode(0,0,0);
	gmds::Node n1 = m3d.newNode(1,0,0);
	gmds::Node n2 = m3d.newNode(1,1,0);
	gmds::Node n3 = m3d.newNode(0,1,0);
	gmds::Node n4 = m3d.newNode(0,0,1);
	gmds::Node n5 = m3d.newNode(1,0,1);
	gmds::Node n6 = m3d.newNode(1,1,1);
	gmds::Node n7 = m3d.newNode(0,1,1);

	m3d.newHex(n0,n1,n2,n3,n4,n5,n6,n7);

	gmds::Node n8 = m3d.newNode(2,0,0);
	gmds::Node n9 = m3d.newNode(2,1,0);
	gmds::Node n10 = m3d.newNode(2,0,1);
	gmds::Node n11 = m3d.newNode(2,1,1);

	m3d.newHex(n1,n8,n9,n2,n5,n10,n11,n6);

	gmds::Node n12 = m3d.newNode(1,0,2);
	gmds::Node n13 = m3d.newNode(2,0,2);
	gmds::Node n14 = m3d.newNode(2,1,2);
	gmds::Node n15 = m3d.newNode(1,1,2);

	m3d.newHex(n5,n10,n11,n6,n12,n13,n14,n15);

	gmds::blocking::Blocking bl3d;
	bl3d.createBlocks3dFromMesh(m3d);
	//ASSERT_EQ(bl3d.nbVertices(), 16);
}
/*----------------------------------------------------------------------------*/
TEST(BlockingTestSuite, unsew_and_vertices_move)
{
	gmds::blocking::Blocking bl;

	gmds::blocking::Dart_handle dh1=
		bl.lcc()->make_hexahedron(gmds::blocking::Point(0,0,0), gmds::blocking::Point(5,0,0),
										  gmds::blocking::Point(5,5,0), gmds::blocking::Point(0,5,0),
										  gmds::blocking::Point(0,5,4), gmds::blocking::Point(0,0,4),
										  gmds::blocking::Point(5,0,4), gmds::blocking::Point(5,5,4));
	gmds::blocking::Dart_handle dh2=
		bl.lcc()->make_hexahedron(gmds::blocking::Point(5,0,0), gmds::blocking::Point(10,0,0),
										  gmds::blocking::Point(10,5,0), gmds::blocking::Point(5,5,0),
										  gmds::blocking::Point(5,5,4), gmds::blocking::Point(5,0,4),
										  gmds::blocking::Point(10,0,4), gmds::blocking::Point(10,5,4));

	bl.writeVTKFile("unsew_a.vtk");
	bl.lcc()->sew<3>(bl.lcc()->alpha(dh1, 1,0,1,2), bl.lcc()->alpha(dh2,2));
	bl.writeVTKFile("unsew_b.vtk");
	bl.lcc()->unsew<3>(bl.lcc()->alpha(dh1, 1,0,1,2));
	bl.writeVTKFile("unsew_c.vtk");
}
/*----------------------------------------------------------------------------*/