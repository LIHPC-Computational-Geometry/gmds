/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/Blocking2D.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
/*----------------------------------------------------------------------------*/
TEST(GridBuildOpClass, test2D)
{
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::F|gmds::F2N));

    gmds::GridBuilder gb(&m,2);

    ASSERT_TRUE(gb.isValid());

    gb.execute(3,1.0, 4, 1.0);

    ASSERT_EQ(m.getNbNodes(),12);
    ASSERT_EQ(m.getNbFaces(),6);

}
/*----------------------------------------------------------------------------*/
TEST(GridBuildOpClass, test3D)
{
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::N|gmds::R|gmds::R2N));

    gmds::GridBuilder gb(&m,3);

    ASSERT_TRUE(gb.isValid());

    gb.execute(3,1.0, 4, 1.0, 3, 2.0);

    ASSERT_EQ(m.getNbNodes(),36);
    ASSERT_EQ(m.getNbRegions(),12);

}
/*----------------------------------------------------------------------------*/
TEST(GridTestSuite, test_blocking2D_output)
{
    Blocking2D m;
    Node n1 = m.newBlockCorner(0,0);
    Node n2 = m.newBlockCorner(1,0);
    Node n3 = m.newBlockCorner(1,1);
    Node n4=  m.newBlockCorner(0,1);

    Blocking2D::Block b1 = m.newBlock(n1,n2,n3,n4);

    Node n5 = m.newBlockCorner(2,0,0);
    Node n6 = m.newBlockCorner(2,1.5,0);
    Blocking2D::Block b2 = m.newBlock(n2,n5,n6,n3);
    b1.setNbDiscretizationI(10);
    b1.setNbDiscretizationJ(10);
    b2.setNbDiscretizationI(10);
    b2.setNbDiscretizationJ(10);

    m.initializeGridPoints();

    IGMeshIOService ios(&m);
    VTKWriter writer(&ios);
    writer.setCellOptions(N|F);
    writer.setDataOptions(N|F);
    writer.write("blocking2D_sample.vtk");

}
/*----------------------------------------------------------------------------*/
