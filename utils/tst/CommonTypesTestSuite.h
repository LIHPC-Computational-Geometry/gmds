#ifndef COMMON_TYPES_TEST_SUITE_H
#define COMMON_TYPES_TEST_SUITE_H

#include "gmds/utils/CommonTypes.h"
#include "gtest/gtest.h"
#include <unit_test_config.h>

// Test for getCommonBut
TEST(CommonTypesTestSuite, GetCommonButTest)
{
	std::vector<TCellID> set1 = {1, 2, 3, 4, 5};
	std::vector<TCellID> set2 = {3, 4, 5, 6, 7};
	TCellID but = 2;

	std::vector<TCellID> result = gmds::getCommonBut(set1, set2, but);

	ASSERT_EQ(result.size(), 2);
	ASSERT_TRUE(std::find(result.begin(), result.end(), 3) != result.end());
	ASSERT_TRUE(std::find(result.begin(), result.end(), 5) != result.end());
}

// Test for keepFilter
TEST(CommonTypesTestSuite, KeepFilterTest)
{
	std::vector<TCellID> set = {1, 2, 2, 3, 3, 3, 4, 4, 4, 4};
	TInt nb = 2;

	std::vector<TCellID> result = gmds::keepFilter(set, nb);

	ASSERT_EQ(result.size(), 2);
	ASSERT_TRUE(std::find(result.begin(), result.end(), 2) != result.end());
	ASSERT_TRUE(std::find(result.begin(), result.end(), 3) != result.end());
}

// Ajoutez d'autres tests pour les fonctions de CommonTypes.cpp

#endif     // COMMON_TYPES_TEST_SUITE_H
