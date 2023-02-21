/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Files containing the different test suites to launch
#include "ReaderTestSuite.h"
#include "WriterTestSuite.h"
#ifdef WITH_LIMA
#include "LimaWriteAndReadTestSuite.h"
#endif//WITH_LIMA
//#include "MeshBWriterTestSuite.h"
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
/*----------------------------------------------------------------------------*/

