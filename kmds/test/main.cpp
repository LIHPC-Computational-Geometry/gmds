/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Files containing the different test suites to launch
#include "ConnectivityPerfTest.h"
#include "ConnectivityTest.h"
#include "ContainersTest.h"
#include "EContainerTest.h"
#include "FaceTest.h"
#include "FakeTypesTest.h"
#include "FContainersTest.h"
#include "GraphTest.h"
#include "GrowingViewTest.h"
#include "InitToolsTest.h"
#include "KokkosUnorderedMapTest.h"
#include "MeshTest.h"
#include "NContainerTest.h"
#include "RContainersTest.h"
#include "RegionTest.h"
#include "VariableTest.h"
#include "WriterTest.h"
/*----------------------------------------------------------------------------*/
int
main(int argc, char** argv)
{
        ::testing::InitGoogleTest(&argc, argv);
        return RUN_ALL_TESTS();
}
/*----------------------------------------------------------------------------*/
