/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Files containing the different test suites to launch

#include "BoundaryOpTestSuite.h"
#include "BoundaryExtractorTestSuite.h"
#include "GridBuilderTestSuite.h"
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
/*----------------------------------------------------------------------------*/

