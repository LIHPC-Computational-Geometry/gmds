//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/geodHoneyComb/GeodHexMesher.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(GeodHoneyCombestClass, test_geode1)
{
    GeodHexMesher ghm;
    ghm.setCenter(math::Point(0,0,0));
    ghm.setRadius(1);
    math::DiscretizationScheme1DUniform d(math::Point(0,0,0),
                                          math::Point(1,0,0),
                                          10);
    math::DiscretizationScheme1DGeometric d10(0.4,math::Point(0,0,0),
                                          math::Point(1,0,0),
                                          10);
    auto layer0 = std::make_pair(0., &d);
    auto layer1 = std::make_pair(0.5, &d);
    auto layer2 = std::make_pair(0.8,&d);
    auto layer3 = std::make_pair(0.9,&d);
    std::map<double, math::DiscretizationScheme1D*> layers;
    layers.insert(layer1);
    layers.insert(layer2);
    layers.insert(layer3);
    layers.insert(layer0);
    ghm.setLayerData(layers);
    ASSERT_EQ(GeodHexMesher::GEOD_SUCCESS, ghm.execute());

    std::unique_ptr<Mesh> m = ghm.getMesh();
    IGMeshIOService ioService(m.get());
    VTKWriter writer(&ioService);
    writer.setCellOptions(N|R);
    writer.setDataOptions(N|R);
    writer.write("geode.vtk");

}
/*----------------------------------------------------------------------------*/
