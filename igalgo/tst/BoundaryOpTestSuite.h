/*----------------------------------------------------------------------------*/
//
// Created by ledouxf on 1/25/19.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <iostream>
#include <gmds/io/IGMeshIOService.h>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
TEST(BoundaryOpClass, testBoundary2D)
{
    // WE WRITE
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|
    gmds::F2N|gmds::N2F|gmds::E|gmds::E2N));

    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/triangulated_quad.vtk";
    gmds::IGMeshIOService ioService(&m);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::F);
    vtkReader.read(vtk_file);

    gmds::MeshDoctor doc(&m);
    doc.updateUpwardConnectivity();
    doc.buildBoundaryCells();

    gmds::BoundaryOperator2D bop(&m);

    ASSERT_TRUE(bop.isValid());
    int mark_edge_curve = m.newMark<gmds::Edge>();
    int mark_node_curve = m.newMark<gmds::Node>();
    bop.markCellsOnCurves(mark_edge_curve,mark_node_curve);
    int e=0,n=0;
    for(auto i:m.nodes()){
        if(m.isMarked<gmds::Node>(i,mark_node_curve))
            n++;
    }
    for(auto i:m.edges()){
        if(m.isMarked<gmds::Edge>(i,mark_edge_curve))
            e++;
    }
    ASSERT_TRUE(n!=0);
    ASSERT_TRUE(e!=0);

    m.unmarkAll<gmds::Edge>(mark_edge_curve);
    m.unmarkAll<gmds::Node>(mark_node_curve);
    m.freeMark <gmds::Edge>(mark_edge_curve);
    m.freeMark <gmds::Node>(mark_node_curve);

}