#ifndef GMDS_COMMON_TYPES_TESTSUITE_H
#define GMDS_COMMON_TYPES_TESTSUITE_H

#include "gtest/gtest.h"
#include <gmds/utils/CommonTypes.h>
#include <unit_test_config.h>
using namespace gmds;

// Test for getCommonBut
TEST(CommonTypesTestSuite, GetCommonButTest)
{
	std::vector<gmds::TCellID> set1 = {1, 2, 3, 4, 5};
	std::vector<gmds::TCellID> set2 = {3, 4, 5, 6, 7};
	gmds::TCellID but = 4;

	std::vector<gmds::TCellID> result = gmds::getCommonBut(set1, set2, but);

	ASSERT_EQ(result.size(), 2);
	ASSERT_TRUE(std::find(result.begin(), result.end(), 2) != result.end());
	ASSERT_TRUE(std::find(result.begin(), result.end(), 3) != result.end());
	ASSERT_TRUE(std::find(result.begin(), result.end(), 4) == result.end());
}

TEST(CommonTypesTestSuite, KeepFilterTest)
{
	std::vector<gmds::TCellID> set = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4};
	gmds::TInt nb = 2;

	std::vector<gmds::TCellID> result = gmds::keepFilter(set, nb);

	ASSERT_EQ(result.size(), 2);
	ASSERT_TRUE(std::find(result.begin(), result.end(), 3) != result.end());
	ASSERT_TRUE(std::find(result.begin(), result.end(), 4) != result.end());
}

#endif     // GMDS_COMMON_TYPES_TEST_SUITE_H
