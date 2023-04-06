/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Files containing the different test suites to launch

//#include "SimplexMeshTestSuite.h"
//#include "SimplexMeshTestAdd_NTt.h"
//#include "SimplexMeshTestReorient.h"
//#include "SimplexMeshTestDelete_NTt.h"
#include "ModeleCAD0.h"
//#include "ModeleCAD1.h"
//#include "ModeleCAD2.h"
//#include "ModeleCAD3.h"
//#include "ModeleCAD4.h"
//#include "ModeleCAD5.h"

/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
/*----------------------------------------------------------------------------*/
