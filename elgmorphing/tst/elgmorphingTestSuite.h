/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
//#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/MdlReader.h>
#include <gtest/gtest.h>
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
TEST(elgmorphingTestSuite, dummytest)
{
	ASSERT_EQ(0, 0);
}
/*----------------------------------------------------------------------------*/
TEST(elgmorphingTestSuite, testGrid2D)
{
	gmds::Mesh m2(gmds::MeshModel(gmds::DIM2|gmds::E|gmds::N|gmds::E2N));

	std::string dir(TEST_SAMPLES_DIR);

	gmds::MdlReader reader2(m2);
	reader2.read(dir + "/test4Mdl.mdl");
	reader2.createVariablesFromGroups();

	ASSERT_TRUE(true);
}
/*----------------------------------------------------------------------------*/