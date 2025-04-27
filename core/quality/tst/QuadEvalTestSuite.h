//
// Created by ledouxf on 1/22/19.
//

/*----------------------------------------------------------------------------*/
#include <cstdlib>
#include <gmds/quality/QuadQuality.h>
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
TEST(QuadQualityTestClass, test1)
{
	QuadQuality qe({0,0,0,
					 1,0,0,
					 1,1,0,
					 0,1,0});


	ASSERT_NEAR(qe.p[0][0],0,0.001);
	ASSERT_NEAR(qe.p[0][1],0,0.001);
	ASSERT_NEAR(qe.p[0][2],0,0.001);

	ASSERT_NEAR(qe.jacobian(),1,0.01);
	ASSERT_NEAR(qe.scaledJacobian(),1,0.01);
	ASSERT_NEAR(qe.aspectRatio(),1,0.01);
	ASSERT_NEAR(qe.signedArea(),1,0.01);
	ASSERT_NEAR(qe.angleDeviation(),0,0.01);
	ASSERT_NEAR(qe.skew(),0,0.01);

	QuadQuality qe2 = QuadQuality::build(math::Point(0,0,0),
	                               math::Point(2,0,0),
	                               math::Point(2,1,0),
	                               math::Point(0,1,0));


	ASSERT_NEAR(qe2.p[0][0],0,0.001);
	ASSERT_NEAR(qe2.p[0][1],0,0.001);
	ASSERT_NEAR(qe2.p[0][2],0,0.001);

	ASSERT_NEAR(qe2.jacobian(),2,0.01);
	ASSERT_NEAR(qe2.scaledJacobian(),1,0.01);
	ASSERT_NEAR(qe2.aspectRatio(),6,0.01);
	ASSERT_NEAR(qe2.signedArea(),2,0.01);
	ASSERT_NEAR(qe2.angleDeviation(),0,0.01);
	ASSERT_NEAR(qe2.skew(),0,0.01);

	QuadQuality qe3 = QuadQuality::build(math::Point(0,0,0),
	                                     math::Point(2,0,0),
	                                     math::Point(2,1,0),
	                                     math::Point(1,1,0));
	ASSERT_NEAR(qe3.skew(),0.447,0.01);

}
