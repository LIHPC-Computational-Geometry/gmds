//
// Created by rochec on 14/01/2022.
//

#include <gmds/claire/DistanceMap.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gtest/gtest.h>
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

TEST(DistanceMapTestClass, DistanceMap_Test1)
{

	DistanceMap distmap ;
	DistanceMap::STATUS result = distmap.execute();

	ASSERT_EQ(DistanceMap::SUCCESS, result);
}


