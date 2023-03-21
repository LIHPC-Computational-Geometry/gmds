/*----------------------------------------------------------------------------*/
#include <gmds/elgmorphing/ElgMorphing.h>
#include <gmds/ig/Mesh.h>
//#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/MdlReader.h>
#include <gmds/math/Point.h>

#include <gtest/gtest.h>
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
TEST(elgmorphingTestSuite, dummytest)
{
	ASSERT_EQ(0, 0);
}
/*----------------------------------------------------------------------------*/
TEST(elgmorphingTestSuite, bbtransform2d)
{
	gmds::math::Point pt(1,1,1);
	gmds::TCoord minXYZ_orig[3] = {-1,-2,-3};
	gmds::TCoord maxXYZ_orig[3] = {2,3,4};
	gmds::TCoord minXYZ_dest[3] = {10,12,13};
	gmds::TCoord maxXYZ_dest[3] = {11,15,15};

	gmds::elgmorphing::ElgMorphing elg;
	gmds::math::Point ptnew = elg.bbtransform2d(pt, minXYZ_orig, maxXYZ_orig, minXYZ_dest, maxXYZ_dest);

	ASSERT_DOUBLE_EQ(ptnew.X(), 10.+1.*(2./3.));
	ASSERT_DOUBLE_EQ(ptnew.Y(), 12.+3.*(3./5.));
//	ASSERT_DOUBLE_EQ(ptnew.Z(), 13.+2.*(4./7.));
}
/*----------------------------------------------------------------------------*/
TEST(elgmorphingTestSuite, boundingbox2d)
{
	gmds::Mesh m_dest(gmds::MeshModel(gmds::DIM2|gmds::E|gmds::N|gmds::E2N));

	std::string dir(TEST_SAMPLES_DIR);

	gmds::MdlReader reader_dest(m_dest);
	reader_dest.read(dir + "/test4Mdl.mdl");
	reader_dest.createVariablesFromGroups();

	gmds::elgmorphing::ElgMorphing elg;
	gmds::TCoord minXYZ_dest[3];
	gmds::TCoord maxXYZ_dest[3];
	elg.computeBoundingBox(&m_dest, std::string("ALU"), minXYZ_dest, maxXYZ_dest);

	gmds::math::Point ptmin(minXYZ_dest[0],minXYZ_dest[1],minXYZ_dest[2]);
	gmds::math::Point ptmax(maxXYZ_dest[0],maxXYZ_dest[1],maxXYZ_dest[2]);

	ASSERT_EQ(gmds::math::Point(0.,-1.,0.), ptmin);
	ASSERT_EQ(gmds::math::Point(4.,2.,0.), ptmax);
}
/*----------------------------------------------------------------------------*/
TEST(elgmorphingTestSuite, morphing2d)
{
	gmds::Mesh m_dest(gmds::MeshModel(gmds::DIM2|gmds::E|gmds::N|gmds::E2N));

	std::string dir(TEST_SAMPLES_DIR);

	gmds::MdlReader reader_dest(m_dest);
	reader_dest.read(dir + "/test4Mdl.mdl");
	reader_dest.createVariablesFromGroups();

	gmds::elgmorphing::ElgMorphing elg;
	gmds::TCoord minXYZ_dest[3];
	gmds::TCoord maxXYZ_dest[3];
	elg.computeBoundingBox(&m_dest, std::string("ALU"), minXYZ_dest, maxXYZ_dest);

	gmds::math::Point ptmin(minXYZ_dest[0],minXYZ_dest[1],minXYZ_dest[2]);
	gmds::math::Point ptmax(maxXYZ_dest[0],maxXYZ_dest[1],maxXYZ_dest[2]);

	ASSERT_EQ(gmds::math::Point(0.,-1.,0.), ptmin);
	ASSERT_EQ(gmds::math::Point(4.,2.,0.), ptmax);
}
/*----------------------------------------------------------------------------*/