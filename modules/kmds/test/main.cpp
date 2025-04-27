/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Kokkos headers
#include <Kokkos_Core.hpp>
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
	Kokkos::InitializationSettings kargs;
	kargs.set_num_threads(3);
	Kokkos::initialize(kargs);

	::testing::InitGoogleTest(&argc, argv);
	int ret = RUN_ALL_TESTS();

	Kokkos::finalize();

	return ret;
}
/*----------------------------------------------------------------------------*/
