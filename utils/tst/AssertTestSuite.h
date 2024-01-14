#ifndef GMDS_ASSERT_TEST_H
#define GMDS_ASSERT_TEST_H

#include "gmds/utils/Assert.h"
#include "gtest/gtest.h"
using namespace gmds;

TEST(AssertTestSuite, SetAndGetAssertMode)
{
	// Test setAssertMode and getAssertMode functions
	setAssertMode(ASSERT_MODE_THROW);
	EXPECT_EQ(getAssertMode(), ASSERT_MODE_THROW);

	setAssertMode(ASSERT_MODE_ABORT);
	EXPECT_EQ(getAssertMode(), ASSERT_MODE_ABORT);
}

TEST(AssertTestSuite, AssertFailedThrow)
{
	// Test assertFailed with ASSERT_MODE_THROW
	setAssertMode(ASSERT_MODE_THROW);

	EXPECT_THROW(assertFailed("1 == 2", "test_file.cpp", 42), GMDSException);
	// Add more tests based on your assertion logic
}

TEST(AssertTestSuite, AssertFailedAbort)
{
	// Test assertFailed with ASSERT_MODE_ABORT
	setAssertMode(ASSERT_MODE_ABORT);

	// Replace the exit(0) test with output to a file or other mechanism
	ASSERT_EXIT(assertFailed("1 == 2", "test_file.cpp", 42), ::testing::ExitedWithCode(0), "");
	// Add more tests based on your assertion logic
}

TEST(AssertTestSuite, AssertRangeFailedThrow)
{
	// Test assertRangeFailed with ASSERT_MODE_THROW
	setAssertMode(ASSERT_MODE_THROW);

	EXPECT_THROW(assertRangeFailed(5.0, 0.0, 2.0, "test_file.cpp", 42), GMDSException);
}

TEST(AssertTestSuite, AssertRangeFailedAbort)
{
	// Test assertRangeFailed with ASSERT_MODE_ABORT
	setAssertMode(ASSERT_MODE_ABORT);

	// Replace the exit(0) test with output to a file or other mechanism
	ASSERT_EXIT(assertRangeFailed(5.0, 0.0, 2.0, "test_file.cpp", 42), ::testing::ExitedWithCode(0), "");
}

TEST(AssertTestSuite, SetAssertModeAndAssertFailed)
{
	// Test setting assertion mode and calling assertFailed
	setAssertMode(ASSERT_MODE_THROW);
	GMDS_ASSERT(1 == 1);     // Should not throw

	setAssertMode(ASSERT_MODE_ABORT);
	// Replace the exit(0) test with output to a file or other mechanism
	ASSERT_EXIT(GMDS_ASSERT(1 == 2), ::testing::ExitedWithCode(0), "");
}

TEST(AssertTestSuite, SetAssertModeAndRangeAssert)
{
	// Test setting assertion mode and calling range assertion
	setAssertMode(ASSERT_MODE_THROW);
	GMDS_RANGE_ASSERT(1.0, 0.0, 2.0);     // Should not throw

	setAssertMode(ASSERT_MODE_ABORT);
	// Replace the exit(0) test with output to a file or other mechanism
	ASSERT_EXIT(GMDS_RANGE_ASSERT(5.0, 0.0, 2.0), ::testing::ExitedWithCode(0), "");
}

#endif     // GMDS_ASSERT_TEST_H
