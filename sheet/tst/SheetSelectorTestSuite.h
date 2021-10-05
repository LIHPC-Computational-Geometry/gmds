/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 1/25/19.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Mesh.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/sheet/Selector2D.h>
#include <gmds/sheet/Selector3D.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(SelectorClass, testSelection3D)
{
    Mesh m(gmds::MeshModel(DIM3|N|R|R2N));

    GridBuilder gb(&m, 3);
    gb.execute(4,1.,4,1.,4,1.);
    ASSERT_EQ(m.getNbRegions(),27);

    Selector3D select(&m);
    ASSERT_TRUE(select.isValid());

    select.execute(0,1);
    ASSERT_EQ(select.getSheetCells().size(),9);
}
/*----------------------------------------------------------------------------*/
TEST(SelectorClass, testSelection2D)
{
    Mesh m(MeshModel(DIM3|N|E|F
                     |F2N |N2F | E2N));

    GridBuilder gb(&m, 2);
    gb.execute(4,1.,4,1.);
    ASSERT_EQ(m.getNbFaces(),9);

    MeshDoctor doc(&m);
    doc.updateUpwardConnectivity();
    doc.buildBoundaryCells();
    Selector2D select(&m);
    ASSERT_TRUE(select.isValid());

    select.execute(0,1);
    ASSERT_EQ(select.getSheetCells().size(),3);
}