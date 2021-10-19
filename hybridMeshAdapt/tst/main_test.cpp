/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Files containing the different test suites to launch

//#include "SimplexMeshTestSuite.h"
//#include "SimplexMeshTestAdd_NTt.h"
//#include "SimplexMeshTestReorient.h"
//#include "SimplexMeshTestDelete_NTt.h"
#include "SimplexReadAndWriteTestSuite.h"

/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
/*----------------------------------------------------------------------------*/
