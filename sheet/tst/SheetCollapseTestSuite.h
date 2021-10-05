/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 1/25/19.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Mesh.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/sheet/Collapse2D.h>
#include <gmds/sheet/Collapse3D.h>
#include <gmds/cad/GeomMeshLinker.h>
#include <gmds/cad/FACManager.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(SheetCollapseClass, DISABLED_testSheetCollapse_2D)
{
    Mesh m(gmds::MeshModel(DIM2|N|E|F|F2N|N2F|E2N));

    GridBuilder gb(&m, 2);
    gb.execute(4,1.,4,1.);
    ASSERT_EQ(m.getNbFaces(),9);
    MeshDoctor doc(&m);
    doc.buildBoundaryCells();
    doc.updateUpwardConnectivity();


    cad::FACManager manager;
    cad::GeomMeshLinker linker;

    manager.initAndLinkFrom2DMesh(&m,&linker);
    //Here, edge and node of m are connected to the geometry
    ASSERT_EQ(m.getNbEdges(),12);
    Collapse2D op(&m,&linker);
    ASSERT_TRUE(op.isValid());

    op.execute(0,1);
    ASSERT_EQ(m.getNbFaces(),6);
    op.execute(4,8);
    ASSERT_EQ(m.getNbFaces(),4);

    op.execute(2,3);
    ASSERT_EQ(m.getNbFaces(),2);
}
/*----------------------------------------------------------------------------*/
TEST(SheetCollapseClass, DISABLED_testSheetCollapse_3D)
{
    Mesh m(gmds::MeshModel(DIM3|N|R|R2N));

    GridBuilder gb(&m, 3);
    gb.execute(4,1.,4,1.,4,1.);
    ASSERT_EQ(m.getNbRegions(),27);

    Collapse3D op(&m);
    ASSERT_TRUE(op.isValid());

    op.execute(0,1);
    ASSERT_EQ(m.getNbRegions(),18);
}
/*----------------------------------------------------------------------------*/
void func(int* a, int* b){
    *a=*a+1;
    *b=*b-1;
    ASSERT_TRUE(*a<*b);
}
TEST(SheetCool, f){
    int x=1;
    int y=4;
    ASSERT_TRUE(x+y==5);
    func(&x,&y);
    ASSERT_TRUE(x+y==5);


}