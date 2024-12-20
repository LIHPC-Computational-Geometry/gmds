/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Files containing the different test suites to launch

#include "AeroTestSuite.h"
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
#include "NodeNeighbourhoodOnFront_3DTestSuite.h"
#include "MeshAlignmentTestSuite.h"
#include "Front_3DTestSuite.h"
#include "Blocking3DTestSuite.h"
#include "TransfiniteInterpolation_3DTestSuite.h"
#include "ComputeBezierDegree_3DTestSuite.h"
#include "ComputeBezierCurveCtrlPtstoInterpolateCurveTestSuite.h"
#include "RefinementBetaBlock3DTestSuite.h"
#include "ComputeBezierCtrlPtstoInterpolateSurfaceTestSuite.h"
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
/*----------------------------------------------------------------------------*/

