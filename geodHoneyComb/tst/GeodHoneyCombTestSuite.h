//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/geodHoneyComb/GeodHexMesher.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(GeodHoneyCombestClass, test_geode1)
{
    GeodHexMesher ghm(1);
    ASSERT_EQ(GeodHexMesher::GEOD_FAILURE, ghm.execute());

}
/*----------------------------------------------------------------------------*/
