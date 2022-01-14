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
	distmap.add(0,1);
	distmap.add(2,2);
	distmap.add(0,3);
	distmap.add(1,4);
	distmap.add(1.2,5);
	std::cout<<distmap<<std::endl;
	ASSERT_EQ(distmap(0).size(), 2);
}


