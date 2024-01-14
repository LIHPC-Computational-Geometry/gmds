/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Files containing the different test suites to launch

#include "BezierSurfaceTestSuite.h"
#include "ChartTestSuite.h"
#include "Cross2DTestSuite.h"
#include "CrossTestSuite.h"
#include "DiscretizationScheme1DTestSuite.h"
#include "MathTestSuite.h"
#include "PointTestSuite.h"
#include "QuaternionTestSuite.h"
#include "OrientationTestSuite.h"
#include "TransfiniteInterpolationTestSuite.h"
#include "BezierCurveTestSuite.h"
#include "BezierTriangleTestSuite.h"
#include "AxisAngleRotationTestSuite.h"
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
/*----------------------------------------------------------------------------*/
