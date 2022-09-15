/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Files containing the different test suites to launch
#include "BlockingTestSuite.h"
#include "CGALTestSuite.h"
#include "InputMarkedDartsTestSuite.h"
#include "SheetInsertTestSuite.h"
#include "WriterDartsVTKTestSuite.h"
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
/*----------------------------------------------------------------------------*/

