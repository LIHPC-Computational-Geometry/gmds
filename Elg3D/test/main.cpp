/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
// Kokkos headers
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// Files containing the different test suites to launch
#include "AssignCellsTest.h"
#include "BadPillowingTest.h"
#include "BoundingBoxGeomAssociationTest.h"
#include "DefeaturingTest.h"
#include "DummyTest.h"
#include "EigenTest.h"
#include "ExtractGeomTest.h"
#include "FracPresEnforcementTest.h"
#include "FracPresTest.h"
#include "GCOTest.h"
#include "InihTest.h"
#include "InitDataTest.h"
#include "InterfaceNodesPosTest.h"
#include "InterfaceNodesPosSmoothVFTest.h"
#include "ManifoldDetectionTest.h"
#include "MaterialAssignmentTest.h"
#include "MaterialGradientComputationTest.h"
#include "MaterialInterfacesTest.h"
#include "MoveToNewPosTest.h"
#include "OptimizationSmoothTest.h"
#include "PillowingTest.h"
#include "Refinement.h"
#include "SmartLaplacianTest.h"
#include "SubsetProblemTest.h"
#include "ToolsTest.h"

#include <Eigen/Core>
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv) {

	Kokkos::InitializationSettings kargs;
	kargs.set_num_threads(1);
	Kokkos::initialize(kargs);

  ::testing::InitGoogleTest(&argc, argv);

  Eigen::initParallel();

  int ret = RUN_ALL_TESTS();

  Kokkos::finalize();
  return ret;
}
/*----------------------------------------------------------------------------*/

