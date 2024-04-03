#include <gtest/gtest.h>
#include <gmds/milp/milp.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

TEST(DummyTestClass, aaa)
{
	gmds::milp::milp m;
	ASSERT_EQ(0,0);
}
