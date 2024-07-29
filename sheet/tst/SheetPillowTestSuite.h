/*----------------------------------------------------------------------------*/
#include <gmds/cadfac/FACManager.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/sheet/Pillow2D.h>
#include <gmds/sheet/Pillow3D.h>
#include <gmds/smoothy/LaplacianSmoother3C.h>
#include <gtest/gtest.h>
#include <iostream>

using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(PillowOpClass, test3D_1)
{
    Mesh m(MeshModel(DIM3|N|R|R2N));

    Pillow3D pillow(&m);

    ASSERT_TRUE(pillow.isValid());

    GridBuilder gb(&m, 3);
    gb.execute(3,1.,3,1.,3,1.);
    ASSERT_EQ(m.getNbRegions(),8);


    std::vector<TCellID> shrink_set;
    shrink_set.push_back(1);
    pillow.execute(shrink_set);



    ASSERT_EQ(m.getNbRegions(),11);
}
/*----------------------------------------------------------------------------*/
TEST(PillowOpClass, test3D_2)
{
    Mesh m(MeshModel(DIM3|N|R|R2N));

    Pillow3D pillow(&m, false);

    ASSERT_TRUE(pillow.isValid());

    GridBuilder gb(&m, 3);
    gb.execute(3,1.,3,1.,3,1.);
    ASSERT_EQ(m.getNbRegions(),8);

    std::vector<TCellID> shrink_set;
    shrink_set.push_back(0);
    pillow.execute(shrink_set);

    ASSERT_EQ(m.getNbRegions(),14);
}
/*----------------------------------------------------------------------------*/
TEST(PillowOpClass, test3D_3)
{
    Mesh m(MeshModel(DIM3|N|R|R2N|N2R));

    GridBuilder gb(&m, 3);
    gb.execute(3,1.,3,1.,3,1.);
    ASSERT_EQ(m.getNbRegions(),8);

    Pillow3D pillow(&m);
    ASSERT_TRUE(pillow.isValid());

    std::vector<TCellID> shrink_set;
    shrink_set.push_back(0);
    pillow.execute(shrink_set);
    //check right N2R connectivity
    MeshDoctor doc(&m);
    doc.updateUpwardConnectivity();

    Node n = m.get<Node>(0);
    ASSERT_EQ(n.get<Region>().size(),1);
    n = m.get<Node>(1);
    ASSERT_EQ(n.get<Region>().size(),2);
    n = m.get<Node>(3);
    ASSERT_EQ(n.get<Region>().size(),2);
    n = m.get<Node>(9);
    ASSERT_EQ(n.get<Region>().size(),2);
    n = m.get<Node>(4);
    ASSERT_EQ(n.get<Region>().size(),3);
    n = m.get<Node>(10);
    ASSERT_EQ(n.get<Region>().size(),3);
    n = m.get<Node>(12);
    ASSERT_EQ(n.get<Region>().size(),3);

    ASSERT_EQ(m.getNbRegions(),11);

}/*----------------------------------------------------------------------------*/
TEST(PillowOpClass, test3D_4)
{
    Mesh m(MeshModel(DIM3|N|R|R2N));

    GridBuilder gb(&m, 3);
    gb.execute(3,1.,3,1.,3,1.);
    ASSERT_EQ(m.getNbRegions(),8);

    Pillow3D pillow(&m, false);
    ASSERT_TRUE(pillow.isValid());

    std::vector<TCellID> shrink_set;
    shrink_set.push_back(0);
    shrink_set.push_back(1);
    pillow.execute(shrink_set);

    ASSERT_EQ(m.getNbRegions(),18);
}
/*----------------------------------------------------------------------------*/
TEST(PillowOpClass, DISABLED_test3D_5)
{
    Mesh m(MeshModel(DIM3|N|R|R2N));

    GridBuilder gb(&m, 3);
    gb.execute(3,1.,3,1.,3,1.);
    ASSERT_EQ(m.getNbRegions(),8);

    Pillow3D pillow(&m, true);
    ASSERT_TRUE(pillow.isValid());

    std::vector<TCellID> shrink_set;
    shrink_set.push_back(1);
    shrink_set.push_back(3);
    shrink_set.push_back(2);
    pillow.execute(shrink_set);


    IGMeshIOService ioService(&m);
    std::string vtk_file2 = ("pillow_out.vtk");

    VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(N|R);
    vtkWriter.setDataOptions(N|R);
    vtkWriter.write(vtk_file2);
    ASSERT_EQ(m.getNbRegions(),22);
}
/*----------------------------------------------------------------------------*/
TEST(PillowOpClass, test3D_6)
{
    Mesh m(MeshModel(DIM3|N|R|R2N));

    GridBuilder gb(&m, 3);
    gb.execute(3,1.,3,1.,3,1.);
    ASSERT_EQ(m.getNbRegions(),8);

    Pillow3D pillow(&m);
    ASSERT_TRUE(pillow.isValid());


    std::vector<TCellID> shrink_set;
    shrink_set.push_back(0);
    shrink_set.push_back(1);
    shrink_set.push_back(2);
    shrink_set.push_back(3);
    shrink_set.push_back(4);
    shrink_set.push_back(5);
    shrink_set.push_back(6);
    pillow.execute(shrink_set);

    ASSERT_EQ(m.getNbRegions(),11);
}
/*----------------------------------------------------------------------------*/
TEST(PillowOpClass, test3D_7)
{
    Mesh m(MeshModel(DIM3|N|R|R2N));

    GridBuilder gb(&m, 3);
    gb.execute(3,1.,3,1.,3,1.);
    ASSERT_EQ(m.getNbRegions(),8);

    Pillow3D pillow(&m, false);//false
    ASSERT_TRUE(pillow.isValid());


    std::vector<TCellID> shrink_set;
    shrink_set.push_back(0);
    shrink_set.push_back(1);
    shrink_set.push_back(2);
    shrink_set.push_back(3);
    shrink_set.push_back(4);
    shrink_set.push_back(5);
    shrink_set.push_back(6);

    pillow.execute(shrink_set);

    ASSERT_EQ(m.getNbRegions(),32);
}
/*----------------------------------------------------------------------------*/
TEST(PillowOpClass, test2D_1)
{
    Mesh m(MeshModel(DIM3|N|E|F
                                 |F2N |N2F | E2N));

    GridBuilder gb(&m, 2);
    gb.execute(3,1.,3,1.);
    ASSERT_EQ(m.getNbFaces(),4);

    MeshDoctor doc(&m);
    doc.updateUpwardConnectivity();
    doc.buildBoundaryCells();
    //We have 8 boundary edges
    ASSERT_EQ(m.getNbEdges(),8);

    Pillow2D pillow(&m, false);//false
    ASSERT_TRUE(pillow.isValid());

    std::vector<TCellID> shrink_set;
    shrink_set.push_back(0);
    pillow.execute(shrink_set);

    ASSERT_EQ(m.getNbFaces(),8);
    //we check the number of boundary edges
    ASSERT_EQ(m.getNbEdges(),8);
}
/*----------------------------------------------------------------------------*/
TEST(PillowOpClass, test2D_2)
{
    Mesh m(MeshModel(DIM3|N|E|F
                     |F2N |N2F | E2N));

    GridBuilder gb(&m, 2);
    gb.execute(3,1.,3,1.);
    ASSERT_EQ(m.getNbFaces(),4);

    MeshDoctor doc(&m);
    doc.updateUpwardConnectivity();
    doc.buildBoundaryCells();
    Pillow2D pillow(&m, true);
    ASSERT_TRUE(pillow.isValid());

    std::vector<TCellID> shrink_set;
    shrink_set.push_back(0);
    pillow.execute(shrink_set);

    ASSERT_EQ(m.getNbFaces(),6);
}
/*----------------------------------------------------------------------------*/
TEST(PillowOpClass, test2D_3)
{
    Mesh m(MeshModel(DIM3|N|E|F
                     |F2N |N2F | E2N));

    GridBuilder gb(&m, 2);
    gb.execute(3,1.,3,1.);
    ASSERT_EQ(m.getNbFaces(),4);

    MeshDoctor doc(&m);
    doc.updateUpwardConnectivity();
    doc.buildBoundaryCells();
    Pillow2D pillow(&m, true);
    ASSERT_TRUE(pillow.isValid());

    std::vector<TCellID> shrink_set;
    shrink_set.push_back(0);
    shrink_set.push_back(1);
    pillow.execute(shrink_set);

    ASSERT_EQ(m.getNbFaces(),6);
}
/*----------------------------------------------------------------------------*/
TEST(PillowOpClass, test2D_4)
{
    Mesh m(MeshModel(DIM3|N|E|F
                     |F2N |N2F | E2N));

    GridBuilder gb(&m, 2);
    gb.execute(3,1.,3,1.);
    ASSERT_EQ(m.getNbFaces(),4);

    MeshDoctor doc(&m);
    doc.updateUpwardConnectivity();
    doc.buildBoundaryCells();
    Pillow2D pillow(&m, false);
    ASSERT_TRUE(pillow.isValid());

    std::vector<TCellID> shrink_set;
    shrink_set.push_back(0);
    shrink_set.push_back(1);
    pillow.execute(shrink_set);

    ASSERT_EQ(m.getNbFaces(),10);
}

/*----------------------------------------------------------------------------*/
TEST(PillowOpClass, test_pillow3D_with_geom)
{
    Mesh m_vol(gmds::MeshModel(DIM3|R|F|E|N|
                               R2N|R2F|R2E|
                               F2N|F2R|F2E|
                               E2F|E2N|N2E|N2R|N2F));

    gmds::GridBuilder gb(&m_vol, 3);
    gb.execute(3,1.,3,1.,3,1.);
    ASSERT_EQ(m_vol.getNbRegions(),8);


    gmds::Pillow3D pillow(&m_vol, true);
    ASSERT_TRUE(pillow.isValid());

    std::vector<TCellID> shrink_set;
    shrink_set.push_back(2);
    shrink_set.push_back(1);
    shrink_set.push_back(3);
    pillow.execute(shrink_set);

    MeshDoctor doc(&m_vol);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    cad::FACManager manager;
    cad::GeomMeshLinker linker;

    manager.initAndLinkFrom3DMesh(&m_vol,&linker);

    smoothy::LaplacianSmoother3C smoother(&m_vol, &linker);
	 smoother.setNbIterations(10);
    smoother.smoothCurves();
    smoother.smoothSurfaces();
    smoother.smoothVolumes();

    std::string vtk_file2 = ("hexa_pillow_out.vtk");
    IGMeshIOService ioService(&m_vol);

    VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::R);
    vtkWriter.setDataOptions(gmds::N|gmds::R);
    vtkWriter.write(vtk_file2);
}