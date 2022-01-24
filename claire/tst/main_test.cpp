/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Files containing the different test suites to launch

#include "ClaireTestSuite.h"
#include "LevelSet2DTestSuite.h"
#include "DistanceMapTestSuite.h"
#include "GradientComputation2DTestSuite.h"
#include "GradientComputation3DTestSuite.h"
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
/*----------------------------------------------------------------------------*/

