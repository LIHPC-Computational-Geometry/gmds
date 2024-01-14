/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Files containing the different test suites to launch

#include "BooleanMarkTestSuite.h"
#include "Blocking2DTestSuite.h"
#include "EdgeTestSuite.h"
#include "MeshDoctorTestSuite.h"
#include "MeshTestSuite.h"
#include "CellTestSuite.h"
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
/*----------------------------------------------------------------------------*/

