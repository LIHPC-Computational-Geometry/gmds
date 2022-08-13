//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/cad/GeomMeshLinker.h>
#include <gmds/smoothy/AngleBasedQuadSmoother.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(AngleBasedQuadSmootherTestSuite, s24_quad_smoother)
{
    // WE WRITE
    Mesh m_vol(gmds::MeshModel(DIM3|R|F|E|N|
                               R2N|R2F|R2E|
                               F2N|F2R|F2E|
                               E2F|E2N|N2E|N2R|N2F));

    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/Notch/notch_hexa.vtk";
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

    smoothy::AngleBasedQuadSmoother smoother(&linker);
    ASSERT_TRUE(smoother.isValid());

    smoother.smooth(100);


    std::string vtk_file2 = dir+"/Notch/out.vtk";

    VTKWriter writer(&ioService);
    writer.setCellOptions(gmds::N|gmds::R);
    writer.write(vtk_file2);


}

