/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Files containing the different test suites to launch
#include "ArrayTestSuite.h"
#include "AssertTestSuite.h"
#include "BitVectorTestSuite.h"
#include "CommonTypesTestSuite"
#include "OrientedGraphTestSuite.h"
#include "ParamTestSuite.h"
#include "RandomGeneratorTest.h"
#include "UtilsTestSuite.h"
/*----------------------------------------------------------------------------*/
int
main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
/*----------------------------------------------------------------------------*/
