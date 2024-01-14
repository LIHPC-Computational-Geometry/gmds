#ifndef LOG_TEST_SUITE_H
#	define LOG_TEST_SUITE_H

#	include "gmds/utils/Log.h"
#	include "gmds/utils/LogStream.h"
#	include "gtest/gtest.h"
#	include <unit_test_config.h>

using namespace gmds;

TEST(LogTestSuite, ReportingLevelTest)
{
	LogLevel initialLevel = Log::reportingLevel();

	// Test setting and getting reporting level
	Log::reportingLevel() = LOG_DEBUG;
	EXPECT_EQ(Log::reportingLevel(), LOG_DEBUG);

	// Reset to the initial level
	Log::reportingLevel() = initialLevel;
}

TEST(LogTestSuite, AddStreamTest)
{
	Log log;
	LogStream stream1(LOG_INFO);
	LogStream stream2(LOG_WARNING);

	// Add streams to the log
	log.addStream(stream1);
	log.addStream(stream2);

	// Check if streams are added correctly
	ASSERT_EQ(log.out_.size(), 2);
	EXPECT_EQ(log.out_[0].level(), LOG_INFO);
	EXPECT_EQ(log.out_[1].level(), LOG_WARNING);
}

TEST(LogTestSuite, ClearTest)
{
	Log log;
	LogStream stream1(LOG_INFO);
	LogStream stream2(LOG_WARNING);

	// Add streams to the log
	log.addStream(stream1);
	log.addStream(stream2);

	// Clear the log
	log.clear();

	// Check if the log is cleared
	EXPECT_TRUE(log.out_.empty());
}

TEST(LogTestSuite, OperatorStreamTest)
{
	Log log;
	LogStream stream1(LOG_INFO);
	LogStream stream2(LOG_WARNING);

	// Add streams to the log
	log.addStream(stream1);
	log.addStream(stream2);

	// Test the operator<< with log
	log << "Test message";

	// Verify if the message is written to the streams
	for (const auto &outStream : log.out_) {
		std::ostringstream oss;
		oss << "Test message";
		EXPECT_EQ(outStream.str(), oss.str());
	}
}