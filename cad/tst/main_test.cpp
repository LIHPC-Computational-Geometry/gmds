/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Files containing the different test suites to launch

#include "FACManagerTestSuite.h"
#include "GeomLinkerTestSuite.h"
#include "GeomSmootherTestSuite.h"
#include "GeomTopologyTestSuite.h"
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
/*----------------------------------------------------------------------------*/

