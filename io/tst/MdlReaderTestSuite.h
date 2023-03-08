/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/io/MdlReader.h>
#include <gtest/gtest.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(MdlReaderTestSuite, testMdlReader1)
{
	// WE WRITE
	Mesh m2(MeshModel(DIM2|E|N|E2N));

	std::string dir(TEST_SAMPLES_DIR);

	MdlReader reader2(m2);
	reader2.read(dir + "/test4Mdl.mdl");
	reader2.createVariablesFromGroups();

	ASSERT_EQ(19,m2.getNbNodes());
	ASSERT_EQ(16,m2.getNbEdges());
	ASSERT_EQ(2,m2.getNbGroups<Node>());
	ASSERT_EQ(m2.getGroup<Node>(0)->name(), "FER");
	ASSERT_EQ(m2.getGroup<Node>(1)->name(), "ALU");

}
/*----------------------------------------------------------------------------*/
TEST(MdlReaderTestSuite, testMdlReader2)
{
	// WE WRITE
	Mesh m1(MeshModel(DIM2|E|N|E2N));

	std::string dir(TEST_SAMPLES_DIR);

	MdlReader reader1(m1,"ALU");
	reader1.read(dir + "/test4Mdl.mdl");
	reader1.createVariablesFromGroups();

	ASSERT_EQ(10,m1.getNbNodes());
	ASSERT_EQ(8,m1.getNbEdges());
	ASSERT_EQ(1,m1.getNbGroups<Node>());
	ASSERT_EQ(m1.getGroup<Node>(0)->name(), "ALU");

}