//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/cad/FACManager.h>
#include <gmds/cad/GeomMeshLinker.h>
#include <gmds/cad/GeomSmoother.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/igalgo/GridBuilder.h>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(GeomSmootherTestSuite, tet_in_cube)
{
    // WE WRITE
    Mesh m_vol(gmds::MeshModel(DIM3|R|F|E|N|
                               R2N|R2F|R2E|
                               F2N|F2R|F2E|
                               E2F|E2N|N2E|N2R|N2F));

    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/tet_in_box.vtk";
    IGMeshIOService ioService(&m_vol);

    VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::R);
    vtkReader.read(vtk_file);

    MeshDoctor doc(&m_vol);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    cad::FACManager manager;
    cad::GeomMeshLinker linker;

    manager.initAndLinkFrom3DMesh(&m_vol,&linker);

    cad::GeomSmoother smoother(&linker);

    //We perturb mesh node locations to test the smoothing algorithm
    for(auto n_id:m_vol.nodes()){
        if(n_id==8){
            m_vol.get<Node>(n_id).setPoint(math::Point(5,1,5));
        }
        else if(n_id==10){
            m_vol.get<Node>(n_id).setPoint(math::Point(-5,1,5));
        }
        else if(n_id==20){
            m_vol.get<Node>(n_id).setPoint(math::Point(-2,2,5));
        }
        else if (n_id==28){
            m_vol.get<Node>(n_id).setPoint(math::Point(1,1,3));
        }
    }



    smoother.smoothCurves();

    ASSERT_NEAR(m_vol.get<Node>(8).getPoint().Y(), 0, 0.01);
    ASSERT_NEAR(m_vol.get<Node>(10).getPoint().Y(), 0, 0.01);

    smoother.smoothSurfaces();

    ASSERT_NEAR(m_vol.get<Node>(20).getPoint().X(), -1, 0.01);
    ASSERT_NEAR(m_vol.get<Node>(20).getPoint().Y(),  1, 0.01);

    smoother.smoothVolumes();

    ASSERT_NEAR(m_vol.get<Node>(28).getPoint().X(),  2.08, 0.01);
    ASSERT_NEAR(m_vol.get<Node>(28).getPoint().Y(),  0.64, 0.01);
    ASSERT_NEAR(m_vol.get<Node>(28).getPoint().Z(),  1.97, 0.01);


}

/*----------------------------------------------------------------------------*/
TEST(GeomSmootherTestSuite, test3)
{
    Mesh m_vol(gmds::MeshModel(DIM3|R|F|E|N|
                               R2N|R2F|R2E|
                               F2N|F2R|F2E|
                               E2F|E2N|N2E|N2R|N2F));
    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/hexa.vtk";

    IGMeshIOService ioService(&m_vol);

    VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::R);
    vtkReader.read(vtk_file);

    m_vol.deleteRegion(0);
    m_vol.deleteRegion(1);
    m_vol.deleteRegion(2);
    m_vol.deleteRegion(3);

    MeshDoctor doc(&m_vol);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    cad::FACManager manager;
    cad::GeomMeshLinker linker;

    manager.initAndLinkFrom3DMesh(&m_vol,&linker);

    cad::GeomSmoother smoother(&linker);


    smoother.smoothCurves(10);

    smoother.smoothSurfaces(10);

    smoother.smoothVolumes(10);

    std::string vtk_file2 = ("test_samples/hexa_out.vtk");

    VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::R);
   // vtkWriter.setDataOptions(gmds::N|gmds::R);
    vtkWriter.write(vtk_file2);
    ASSERT_TRUE(true);
}
