//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/geodHoneyComb/RegularIcosahedron.h>
#include <gmds/utils/Array.h>
#include <gmds/math/TransfiniteInterpolation.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(RegularIcosahedronTestClass, test_unit_icosahedron)
{
    math::Point center(0,0,0);
    RegularIcosahedron ico(center,1);
    std::unique_ptr<Mesh> m = ico.getRepresentation();


    ASSERT_EQ(12, m->getNbNodes());
    ASSERT_EQ(20, m->getNbFaces());

    for(auto n_id:m->nodes()){
        math::Point p = m->get<Node>(n_id).point();
        ASSERT_FLOAT_EQ(1,p.distance(center));
    }
}
/*----------------------------------------------------------------------------*/
TEST(RegularIcosahedronTestClass, test_unit_icosahedron_subdivide)
{
    math::Point center(0,0,0);
    RegularIcosahedron ico(center,1,20);
    std::unique_ptr<Mesh> m = ico.getRepresentation();


//    ASSERT_EQ(12, m->getNbNodes());
    ASSERT_EQ(20, m->getNbFaces());

    for(auto n_id:m->nodes()){
        math::Point p = m->get<Node>(n_id).point();
    //    ASSERT_FLOAT_EQ(1,p.distance(center));
    }

    IGMeshIOService ioService(m.get());
    VTKWriter writer(&ioService);
    writer.setCellOptions(N|F);
    writer.setDataOptions(N|F);
    writer.write("icosahedron_subdvide.vtk");
}
/*----------------------------------------------------------------------------*/
TEST(RegularIcosahedronTestClass, test_icosahedron_length_3)
{
    math::Point center(0,0,0);
    RegularIcosahedron ico(center,3);
    std::unique_ptr<Mesh> m = ico.getRepresentation();



    ASSERT_EQ(12, m->getNbNodes());
    ASSERT_EQ(20, m->getNbFaces());

    for(auto n_id:m->nodes()){
        math::Point p = m->get<Node>(n_id).point();
        ASSERT_FLOAT_EQ(3,p.distance(center));
    }

}
/*----------------------------------------------------------------------------*/
TEST(RegularIcosahedronTestClass, test_icosahedron_length_2_center555)
{
    math::Point center(5,5,5);
    RegularIcosahedron ico(center,2);
    std::unique_ptr<Mesh> m = ico.getRepresentation();

    ASSERT_EQ(12, m->getNbNodes());
    ASSERT_EQ(20, m->getNbFaces());

    for(auto n_id:m->nodes()){
        math::Point p = m->get<Node>(n_id).point();
        ASSERT_FLOAT_EQ(2,p.distance(center));
    }

    math::Point p1 = m->get<Node>(1).point();
    math::Point p2 = m->get<Node>(2).point();

    ASSERT_FLOAT_EQ(0,center.distance(0.5*(p1+p2)));

}
/*----------------------------------------------------------------------------*/
TEST(RegularIcosahedronTestClass, test_icosahedron_length_2_center555_dualized)
{
    math::Point center(5,5,5);
    RegularIcosahedron ico(center,2);
    ico.performQuadDualization();
    std::unique_ptr<Mesh> m = ico.getRepresentation();


    for(auto n_id:m->nodes()){
        math::Point p = m->get<Node>(n_id).point();
        ASSERT_FLOAT_EQ(2,p.distance(center));
    }

    math::Point p1 = m->get<Node>(1).point();
    math::Point p2 = m->get<Node>(2).point();

    ASSERT_FLOAT_EQ(0,center.distance(0.5*(p1+p2)));

    IGMeshIOService ioService(m.get());
    VTKWriter writer(&ioService);
    writer.setCellOptions(N|F);
    writer.setDataOptions(N|F);
    writer.write("icosahedron.vtk");
}
/*----------------------------------------------------------------------------*/
