//
// Created by ledouxf on 1/22/19.
//

/*----------------------------------------------------------------------------*/
#include <cstdlib>
#include <gmds/quality/HexQuality.h>
#include <gtest/gtest.h>
#include <iostream>
#ifdef _WIN32
#include <io.h>
#else
#include <dirent.h>
#endif // _WIN32
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::quality;
/*----------------------------------------------------------------------------*/
TEST(HexQualityTestClass, test1)
{
	HexQuality he({0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1});

	ASSERT_NEAR(he.p[0][0], 0, 0.001);
	ASSERT_NEAR(he.p[0][1], 0, 0.001);
	ASSERT_NEAR(he.p[0][2], 0, 0.001);

	ASSERT_NEAR(he.volume(), 1, 0.01);
	ASSERT_NEAR(he.diagonal(), 1, 0.01);
	ASSERT_NEAR(he.jacobian(), 1, 0.01);
	ASSERT_NEAR(he.scaledJacobian(), 1, 0.01);
	ASSERT_NEAR(he.edgeRatio(), 1, 0.01);
	ASSERT_NEAR(he.maximumEdgeRatio(), 1, 0.01);
	ASSERT_NEAR(he.meanAspectFrobenius(), 8, 0.01);
	ASSERT_NEAR(he.maximumAspectFrobenius(), 1, 0.01);
	ASSERT_NEAR(he.stretch(), 1, 0.01);
	ASSERT_NEAR(he.shear(), 1, 0.01);
	ASSERT_NEAR(he.shape(), 1, 0.01);
	ASSERT_NEAR(he.skew(), 0, 0.01);
	ASSERT_NEAR(he.oddy(), 0, 0.01);
	ASSERT_NEAR(he.taper(), 0, 0.01);

}

/*----------------------------------------------------------------------------*/
TEST(HexQualityTestClass, test2)
{
	HexQuality he = HexQuality::build(math::Point(0,0,0),
												 math::Point(2,0,0),
												 math::Point(2,1,0),
												 math::Point(0,1,0),
												 math::Point(0,0,1),
												 math::Point(2,0,1),
												 math::Point(2,1,1),
												 math::Point(0,1,1));



		ASSERT_NEAR(he.p[0][0], 0, 0.001);
		ASSERT_NEAR(he.p[0][1], 0, 0.001);
		ASSERT_NEAR(he.p[0][2], 0, 0.001);

		ASSERT_NEAR(he.volume(), 2, 0.01);
		ASSERT_NEAR(he.diagonal(), 1, 0.01);
		ASSERT_NEAR(he.jacobian(), 2, 0.01);
		ASSERT_NEAR(he.scaledJacobian(), 1, 0.01);
		ASSERT_NEAR(he.edgeRatio(), 2, 0.01);
		ASSERT_NEAR(he.maximumEdgeRatio(), 2, 0.01);
		ASSERT_NEAR(he.meanAspectFrobenius(), 9.797, 0.01);
		ASSERT_NEAR(he.maximumAspectFrobenius(), 1.2247, 0.01);
		ASSERT_NEAR(he.stretch(), 0.707, 0.01);
		ASSERT_NEAR(he.shear(), 2, 0.01);
		ASSERT_NEAR(he.shape(), 0.7937, 0.01);
		ASSERT_NEAR(he.skew(), 0, 0.01);
		ASSERT_NEAR(he.oddy(), 2.3811, 0.01);
		ASSERT_NEAR(he.taper(), 0, 0.01);

}
/*----------------------------------------------------------------------------*/
TEST(HexQualityTestClass, test3)
{
	   HexQuality he = HexQuality::build(math::Point(188,28,29),
	                                     math::Point(191,28,29),
	                                     math::Point(191,31,26),
	                                     math::Point(188,31,26),
	                                     math::Point(188,27,28),
	                                     math::Point(191,27,28),
	                                     math::Point(191,29,25),
	                                     math::Point(188,29,25));

	   ASSERT_NEAR(he.skew(), 0.106, 0.01);

}
