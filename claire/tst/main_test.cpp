/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Files containing the different test suites to launch

#include "ClaireTestSuite.h"
#include "LevelSetTestSuite.h"
#include "DistanceMapTestSuite.h"
#include "GradientComputationTestSuite.h"
#include "PointFollowingVectorField2DTestSuite.h"
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
/*----------------------------------------------------------------------------*/

