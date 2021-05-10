//
// Created by ledouxf on 1/22/19.
//

/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/MeshBWriter.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/math/Point.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(MeshClass, testMeshBWriter_TriQuad)
{
    // WE WRITE
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::F2N));


    gmds::Node n0 = m.newNode(0,0,0);
    gmds::Node n1 = m.newNode(1,0,0);
    gmds::Node n2 = m.newNode(1,1,0);
    gmds::Node n3 = m.newNode(0,1,0);

    gmds::Node n4 = m.newNode(2,0,0);
    gmds::Node n5 = m.newNode(2,1,0);
    gmds::Node n6 = m.newNode(3,0,0);
    gmds::Node n7 = m.newNode(3,1,0);

    gmds::Face f0 = m.newTriangle(n0,n1,n2);
    gmds::Face f1 = m.newTriangle(n0,n2,n3);
    gmds::Face f2 = m.newQuad(n1,n2,n5,n4);
    gmds::Face f3 = m.newQuad(n4,n5,n7,n6);

    gmds::IGMeshIOService ioService(&m);
    gmds::MeshBWriter writer(&ioService);
    writer.setCellOptions(gmds::N|gmds::F);
    writer.setDataOptions(gmds::N|gmds::F);
    writer.write("toto_2D.mesh");
    ASSERT_EQ(0, 0);

}
/*----------------------------------------------------------------------------*/
TEST(MeshClass, testMeshBWriter_Tetra)
{
    // WE WRITE
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::R|gmds::N|gmds::R2N));


    gmds::Node n0 = m.newNode(0,0,0);
    gmds::Node n1 = m.newNode(1,0,0);
    gmds::Node n2 = m.newNode(1,1,0);
    gmds::Node n3 = m.newNode(0,1,0);

    gmds::Node n4 = m.newNode(0,0,1);
    gmds::Node n5 = m.newNode(1,0,1);
    gmds::Node n6 = m.newNode(1,1,1);
    gmds::Node n7 = m.newNode(0,1,1);

    m.newTet(n0,n1,n3,n4);
    m.newTet(n3,n2,n1,n6);
    m.newTet(n3,n1,n4,n6);

    gmds::IGMeshIOService ioService(&m);
    gmds::MeshBWriter writer(&ioService);
    writer.setCellOptions(gmds::N|gmds::R);
    writer.setDataOptions(gmds::N|gmds::R);
    writer.write("toto_3D.mesh");
    ASSERT_EQ(0, 0);

}