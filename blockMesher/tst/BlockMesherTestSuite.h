//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/blockMesher/BlockMesher.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(BlockMesherTestClass, test_1_block_non_classified)
{
    Mesh blocking(MeshModel(DIM3|R|N|E|F|R2N|F2N|E2N|R2F|R2E|F2E));
    Node n0 = blocking.newNode(0,0,0);
    Node n1 = blocking.newNode(1,0,0);
    Node n2 = blocking.newNode(1,1,0);
    Node n3 = blocking.newNode(0,1,0);
    Node n4 = blocking.newNode(0,0,1);
    Node n5 = blocking.newNode(1,0,1);
    Node n6 = blocking.newNode(1,1,1);
    Node n7 = blocking.newNode(0,1,1);
    blocking.newHex(n0,n1,n2,n3,n4,n5,n6,n7);

    MeshDoctor doc(&blocking);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    BlockMesher bm(&blocking,NULL);
    ASSERT_EQ(BlockMesher::SUCCESS,bm.execute(4));

    IGMeshIOService ioService(bm.mesh());

    VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::F);
    vtkWriter.setDataOptions(gmds::N|gmds::F);
    vtkWriter.write("block1_mesh.vtk");
}
/*----------------------------------------------------------------------------*/
TEST(BlockMesherTestClass, test_2_blocks_non_classified)
{
    Mesh blocking(MeshModel(DIM3|R|N|E|F|R2N|F2N|E2N|R2F|R2E|F2E));
    Node n0 = blocking.newNode(0,0,0);
    Node n1 = blocking.newNode(1,0,0);
    Node n2 = blocking.newNode(1,1,0);
    Node n3 = blocking.newNode(0,1,0);
    Node n4 = blocking.newNode(0,0,1);
    Node n5 = blocking.newNode(1,0,1);
    Node n6 = blocking.newNode(1,1,1);
    Node n7 = blocking.newNode(0,1,1);
    Node n8 = blocking.newNode(0,0,2);
    Node n9 = blocking.newNode(1,0,2);
    Node n10= blocking.newNode(1,1,2);
    Node n11= blocking.newNode(0,1,2);
    blocking.newHex(n0,n1,n2,n3,n4,n5,n6,n7);
    blocking.newHex(n4,n5,n6,n7,n8,n9,n10,n11);

    MeshDoctor doc(&blocking);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    BlockMesher bm(&blocking,NULL);
    ASSERT_EQ(BlockMesher::SUCCESS,bm.execute(4));

    IGMeshIOService ioService(bm.mesh());

    VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::F);
    vtkWriter.setDataOptions(gmds::N|gmds::F);
    vtkWriter.write("block2_mesh.vtk");

}
/*----------------------------------------------------------------------------*/
