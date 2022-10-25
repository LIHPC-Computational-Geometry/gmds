/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Files containing the different test suites to launch

#include "ClaireTestSuite.h"
#include "LevelSetTestSuite.h"
#include "DistanceMapTestSuite.h"
#include "GradientComputationTestSuite.h"
#include "PointFollowingVectorFieldTestSuite.h"
#include "AeroPipelineTestSuite.h"
#include "AeroBoundariesTestSuite.h"
#include "UtilsTestSuite.h"
#include "SU2WriterTestSuite.h"
#include "AeroMeshQualityTestSuite.h"
#include "RefinementBetaTestSuite.h"
#include "SmoothLineSweepingTestSuite.h"
#include "FastLocalizeTestSuite.h"
#include "DiffusionEquation2DTestSuite.h"
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
/*----------------------------------------------------------------------------*/

