/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Files containing the different test suites to launch
#include "BlockingClassifierTestSuite.h"
#include "BlockingTestSuite.h"
#include "MCTSTestSuite.h"
#include "DataMCTSTestSuite.h"
#include "GraphTestSuite.h"
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
/*----------------------------------------------------------------------------*/

