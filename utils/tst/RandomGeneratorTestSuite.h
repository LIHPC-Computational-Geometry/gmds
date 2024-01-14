#ifndef GMDS_RANDOM_GENERATOR_TESTSUITE_H
#	define GMDS_RANDOM_GENERATOR_TESTSUITE_H

#	include "gmds/utils/RandomGenerator.h"
#	include "gtest/gtest.h"
#	include <unit_test_config.h>

TEST(RandomGeneratorTestSuite, ValueInRange)
{
	gmds::RandomGenerator randomGen;
	randomGen.init();

	for (int i = 0; i < 1000; ++i) {
		double value = randomGen.value();
		EXPECT_GE(value, 0.0);
		EXPECT_LE(value, 1.0);
	}
}

TEST(RandomGeneratorTestSuite, Initialization)
{
	gmds::RandomGenerator randomGen1;
	gmds::RandomGenerator randomGen2;

	randomGen1.init();
	randomGen2.init();

	// Test if two instances have different seed values
	EXPECT_NE(randomGen1.value(), randomGen2.value());
}
