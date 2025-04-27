/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/polyblock/PolyBlock2D.h>
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(PolyBlock2DTestSuite, testA)
{
    // WE WRITE
    Mesh m(MeshModel(DIM3|F|N|N2E|E2F|F2N|N2F|E|E2N|F2E));

    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/HolesInSquare1.vtk";
    IGMeshIOService ioService(&m);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::F);
    vtkReader.read(vtk_file);

    gmds::MeshDoctor doc(&m);
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();
  //  doc.buildBoundaryCells();

  doc.orient2DFaces();
    gmds::BoundaryOperator2D bop(&m);

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


    PolyBlock2D pb(&m, mark_node_curve, mark_edge_curve);
    pb.execute();
    ASSERT_EQ(1,1);

    VTKWriter w(&ioService);
    w.setCellOptions(gmds::N|gmds::F);
    w.setDataOptions(gmds::N|gmds::F);
    w.write("polycut2d.vtk");

    m.unmarkAll<gmds::Edge>(mark_edge_curve);
    m.unmarkAll<gmds::Node>(mark_node_curve);
    m.freeMark <gmds::Edge>(mark_edge_curve);
    m.freeMark <gmds::Node>(mark_node_curve);


}
