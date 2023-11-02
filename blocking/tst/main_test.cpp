/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Files containing the different test suites to launch
#include "BlockingTestSuite.h"
#include "CurvedBlockingTestSuite.h"
#include "CurvedBlockingClassificationTestSuite.h"
#include "CGALTestSuite.h"
#include "InputMarkedDartsTestSuite.h"
#include "SheetInsertTestSuite.h"
#include "SheetOperationTestSuite.h"
#include "WriterDartsVTKTestSuite.h"
#include "ExecutionActionsTestSuite.h"
#ifdef USE_CGNS
#include "CGNSWriterTestSuite.h"
#endif
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
/*----------------------------------------------------------------------------*/

