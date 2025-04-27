//
// Created by rochec on 14/01/2022.
//

#include <gmds/aero/DistanceMap.h>
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

	ASSERT_EQ(distmap(0).size(), 2);

	// Test des méthodes .add et .getNbrIds
	int nbr_ids;
	distmap.getNbrIds(0, nbr_ids);
	ASSERT_EQ(nbr_ids, 2);

	// Test de la méthode remove
	distmap.remove(0, 1);

	distmap.getNbrIds(0, nbr_ids);
	ASSERT_EQ(nbr_ids, 1);

	// Test de la méthode check
	distmap.add(2,2);

	bool test = distmap.check() ;
	ASSERT_EQ(test, false);

	distmap.remove(2, 2);
	test = distmap.check() ;
	ASSERT_EQ(test, true);

	// Test de la méthode getAndRemoveFirst
	double MinDist ;
	TCellID MinId ;
	distmap.getAndRemoveFirst(MinDist, MinId);
 	ASSERT_EQ(MinDist, 0);
	ASSERT_EQ(MinId, 3);

	// Test de la méthode update
	distmap.update(2, 1.2, 2);
 	distmap.getNbrIds(1.2, nbr_ids);
	ASSERT_EQ(nbr_ids, 2);

	// Test des méthodes .getNbrIdsTotal
	distmap.getNbrIdsTotal(nbr_ids);
	ASSERT_EQ(nbr_ids, 3);

	// Test de la méthode .isEmpty
	ASSERT_EQ(distmap.isEmpty(), false);
}


