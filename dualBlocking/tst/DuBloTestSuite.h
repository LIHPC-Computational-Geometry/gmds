//
// Created by ledouxf on 1/22/19.
//

/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/MeditReader.h>
#include <gmds/io/IGMeshIOService.h>
#include <iostream>
#include <smoothy/inc/gmds/smoothy/LaplacianSmoother.h>
#include <gmds/dualBlocking/DualBlockingSession.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace db;
/*----------------------------------------------------------------------------*/
TEST(DuBloTestClass, test1)
{
    // WE WRITE
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::R|gmds::F|gmds::N));

    std::string vtk_file = ("test_samples/B1_tet.vtk");
    //std::string vtk_file = ("/ccc/home/cont001/ocre/calderans/Dev/gmds/test_samples/B1_tet.vtk");
    gmds::IGMeshIOService ioService(&m);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::R);
    vtkReader.setDataOptions(gmds::N);
    vtkReader.read(vtk_file);

    ASSERT_EQ(m.getNbNodes(),2658);
    ASSERT_EQ(m.getNbRegions(),9872);
}
/*----------------------------------------------------------------------------*/
TEST(DuBloTestClass, test_B7)
{
    std::string param_file("test_samples/duBloInputs/B7.vtk");

    Mesh mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                        R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    IGMeshIOService ioService(&mesh);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::R);
    vtkReader.setDataOptions(gmds::N|gmds::R);
    vtkReader.read(param_file);

    MeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();


    Mesh hmesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                         R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    MeshDoctor doch(&hmesh);
    doch.buildFacesAndR2F();
    doch.buildEdgesAndX2E();
    doch.updateUpwardConnectivity();


    math::Vector3d x(1,0,0);
    math::Vector3d y(0,1,0);
    math::Vector3d z(0,0,1);


    DualBlockingSession session(&mesh,&hmesh);

    ASSERT_TRUE(session.createSurface(12769,x,1));
    ASSERT_TRUE(session.createSurface(7353,y,2));
    ASSERT_TRUE(session.createSurface(7353,z,3));
    ASSERT_TRUE(session.createSurface(18526,z,4));

    ASSERT_TRUE(session.colorDual());

    session.createBlock();

    //Nb blocks in the final primal structure
    ASSERT_EQ(hmesh.getNbRegions(),4);
}
/*----------------------------------------------------------------------------*/
TEST(DuBloTestClass, test_B10)
{
    std::string param_file("test_samples/duBloInputs/B10.vtk");

    Mesh mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                        R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    IGMeshIOService ioService(&mesh);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::R);
    vtkReader.setDataOptions(gmds::N|gmds::R);
    vtkReader.read(param_file);

    MeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();


    Mesh hmesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                         R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    MeshDoctor doch(&hmesh);
    doch.buildFacesAndR2F();
    doch.buildEdgesAndX2E();
    doch.updateUpwardConnectivity();


    math::Vector3d x(1,0,0);
    math::Vector3d y(0,1,0);
    math::Vector3d z(0,0,1);


    DualBlockingSession session(&mesh,&hmesh);

    ASSERT_TRUE(session.createSurface(81912,x,1));
    ASSERT_TRUE(session.createSurface(81912,y,2));
    ASSERT_TRUE(session.createSurface(141649,x,3));
    ASSERT_TRUE(session.createSurface(141649,y,4));
    ASSERT_TRUE(session.createSurface(98280,y,5));
    ASSERT_TRUE(session.createSurface(32943,x,6));
    ASSERT_TRUE(session.createSurface(133927,y,7));
    ASSERT_TRUE(session.createSurface(59617,x,8));
    ASSERT_TRUE(session.createSurface(263729,z,9));
    ASSERT_TRUE(session.createBoundary(188762,10));

    session.refineSurfaceSheet();

    ASSERT_TRUE(session.colorDual());

    session.createBlock();

    //Nb blocks in the final primal structure
    ASSERT_EQ(hmesh.getNbRegions(),12);
}
/*----------------------------------------------------------------------------*/
TEST(DuBloTestClass, test_B8)
{
    std::string param_file("test_samples/duBloInputs/B8.vtk");

    Mesh mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                        R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    IGMeshIOService ioService(&mesh);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::R);
    vtkReader.setDataOptions(gmds::N|gmds::R);
    vtkReader.read(param_file);

    MeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();


    Mesh hmesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                         R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    MeshDoctor doch(&hmesh);
    doch.buildFacesAndR2F();
    doch.buildEdgesAndX2E();
    doch.updateUpwardConnectivity();


    math::Vector3d x(1,0,0);
    math::Vector3d y(0,1,0);
    math::Vector3d z(0,0,1);


    DualBlockingSession session(&mesh,&hmesh);

    ASSERT_TRUE(session.createSurface(24605,y,1));
    ASSERT_TRUE(session.createSurface(24605,z,2));
    ASSERT_TRUE(session.createSurface(62440,y,3));
    ASSERT_TRUE(session.createSurface(62440,z,4));
    ASSERT_TRUE(session.createSurface(16975,x,5));
    ASSERT_TRUE(session.createSurface(67297,x,6));
    ASSERT_TRUE(session.createBoundary(73944,7));

    session.refineSurfaceSheet();

    ASSERT_TRUE(session.colorDual());

    session.createBlock();

    //Nb blocks in the final primal structure
    ASSERT_EQ(hmesh.getNbRegions(),10);
}
/*----------------------------------------------------------------------------*/
TEST(DuBloTestClass, test_B40)
{
    std::string param_file("test_samples/duBloInputs/B40.vtk");

    Mesh mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                        R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    IGMeshIOService ioService(&mesh);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::R);
    vtkReader.setDataOptions(gmds::N|gmds::R);
    vtkReader.read(param_file);

    MeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();


    Mesh hmesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                         R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    MeshDoctor doch(&hmesh);
    doch.buildFacesAndR2F();
    doch.buildEdgesAndX2E();
    doch.updateUpwardConnectivity();


    math::Vector3d x(1,0,0);
    math::Vector3d y(0,1,0);
    math::Vector3d z(0,0,1);


    DualBlockingSession session(&mesh,&hmesh);

    ASSERT_TRUE(session.createSurface(31893,x,1));
    ASSERT_TRUE(session.createSurface(31893,y,2));
    ASSERT_TRUE(session.createSurface(32427,y,3));
    ASSERT_TRUE(session.createSurface(20693,x,4));
    ASSERT_TRUE(session.createSurface(12280,z,5));

    session.refineSurfaceSheet();

    ASSERT_TRUE(session.colorDual());

    session.createBlock();

    //Nb blocks in the final primal structure
    ASSERT_EQ(hmesh.getNbRegions(),3);
}
/*----------------------------------------------------------------------------*/
TEST(DuBloTestClass, test_B45)
{
    std::string param_file("test_samples/duBloInputs/B45.vtk");

    Mesh mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                        R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    IGMeshIOService ioService(&mesh);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::R);
    vtkReader.setDataOptions(gmds::N|gmds::R);
    vtkReader.read(param_file);

    MeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();


    Mesh hmesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                         R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    MeshDoctor doch(&hmesh);
    doch.buildFacesAndR2F();
    doch.buildEdgesAndX2E();
    doch.updateUpwardConnectivity();


    math::Vector3d x(1,0,0);
    math::Vector3d y(0,1,0);
    math::Vector3d z(0,0,1);


    DualBlockingSession session(&mesh,&hmesh);

    ASSERT_TRUE(session.createSurface(66344,x,1));
    ASSERT_TRUE(session.createSurface(66344,y,2));
    ASSERT_TRUE(session.createSurface(128894,x,3));
    ASSERT_TRUE(session.createSurface(146876,x,4));
    ASSERT_TRUE(session.createSurface(99890,y,5));
    ASSERT_TRUE(session.createSurface(142099,y,6));
    ASSERT_TRUE(session.createSurface(123163,z,7));

    session.refineSurfaceSheet();

    ASSERT_TRUE(session.colorDual());

    session.createBlock();

    //Nb blocks in the final primal structure
    ASSERT_EQ(hmesh.getNbRegions(),6);
}
/*----------------------------------------------------------------------------*/
TEST(DuBloTestClass, test_S7)
{
    std::string param_file("test_samples/duBloInputs/S7.vtk");

    Mesh mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                        R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    IGMeshIOService ioService(&mesh);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::R);
    vtkReader.setDataOptions(gmds::N|gmds::R);
    vtkReader.read(param_file);

    MeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();


    Mesh hmesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                         R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    MeshDoctor doch(&hmesh);
    doch.buildFacesAndR2F();
    doch.buildEdgesAndX2E();
    doch.updateUpwardConnectivity();


    math::Vector3d x(1,0,0);
    math::Vector3d y(0,1,0);
    math::Vector3d z(0,0,1);


    DualBlockingSession session(&mesh,&hmesh);

    ASSERT_TRUE(session.createSurface(115356,y,1));
    ASSERT_TRUE(session.createSurface(74918,x,2));
    ASSERT_TRUE(session.createSurface(61147,x,3));
    ASSERT_TRUE(session.createSurface(127654,x,4));
    ASSERT_TRUE(session.createSurface(119658,y,5));
    ASSERT_TRUE(session.createSurface(133739,y,6));
    ASSERT_TRUE(session.createSurface(128446,z,7));
    ASSERT_TRUE(session.createSurface(72460,z,8));
    ASSERT_TRUE(session.createSurface(132286,z,9));
    ASSERT_TRUE(session.createSurface(74954,z,10));
    ASSERT_TRUE(session.createSurface(115029,z,11));
    ASSERT_TRUE(session.createBoundary(95013,12));
    ASSERT_TRUE(session.createBoundary(103210,13));

    session.refineSurfaceSheet();

    ASSERT_TRUE(session.colorDual());

    session.createBlock();

    //Nb blocks in the final primal structure
    ASSERT_EQ(hmesh.getNbRegions(),51);
}
/*----------------------------------------------------------------------------*/
TEST(DuBloTestClass, test_S34)
{
    std::string param_file("test_samples/duBloInputs/S34.vtk");

    Mesh mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                        R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    IGMeshIOService ioService(&mesh);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::R);
    vtkReader.setDataOptions(gmds::N|gmds::R);
    vtkReader.read(param_file);

    MeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();


    Mesh hmesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                         R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    MeshDoctor doch(&hmesh);
    doch.buildFacesAndR2F();
    doch.buildEdgesAndX2E();
    doch.updateUpwardConnectivity();


    math::Vector3d x(1,0,0);
    math::Vector3d y(0,1,0);
    math::Vector3d z(0,0,1);


    DualBlockingSession session(&mesh,&hmesh);

    ASSERT_TRUE(session.createSurface(114611,z,1));
    ASSERT_TRUE(session.createSurface(114611,y,2));
    ASSERT_TRUE(session.createSurface(113011,x,3));
    ASSERT_TRUE(session.createSurface(98691,x,4));
    ASSERT_TRUE(session.createSurface(108064,z,5));
    ASSERT_TRUE(session.createSurface(118850,y,6));
    ASSERT_TRUE(session.createSurface(77644,y,7));
    ASSERT_TRUE(session.createSurface(120621,z,8));
    ASSERT_TRUE(session.createSurface(117077,y,9));
    ASSERT_TRUE(session.createBoundary(75624,10));

    session.refineSurfaceSheet();

    ASSERT_TRUE(session.colorDual());

    session.createBlock();

    //Nb blocks in the final primal structure
    ASSERT_EQ(hmesh.getNbRegions(),17);
}
/*----------------------------------------------------------------------------*/
TEST(DuBloTestClass, test_S38)
{
    std::string param_file("test_samples/duBloInputs/S38.vtk");

    Mesh mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                        R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    IGMeshIOService ioService(&mesh);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::R);
    vtkReader.setDataOptions(gmds::N|gmds::R);
    vtkReader.read(param_file);

    MeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();


    Mesh hmesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                         R2F | F2R | F2E | E2R | E2F | R2E | N2R | N2F | N2E));

    MeshDoctor doch(&hmesh);
    doch.buildFacesAndR2F();
    doch.buildEdgesAndX2E();
    doch.updateUpwardConnectivity();


    math::Vector3d x(1,0,0);
    math::Vector3d y(0,1,0);
    math::Vector3d z(0,0,1);


    DualBlockingSession session(&mesh,&hmesh);

    ASSERT_TRUE(session.createSurface(172121,z,1));
    ASSERT_TRUE(session.createSurface(172121,y,2));
    ASSERT_TRUE(session.createSurface(195020,y,3));
    ASSERT_TRUE(session.createSurface(188465,y,4));
    ASSERT_TRUE(session.createSurface(183340,x,5));
    ASSERT_TRUE(session.createSurface(183340,z,6));
    ASSERT_TRUE(session.createSurface(97514,x,7));
    ASSERT_TRUE(session.createSurface(194196,x,8));
    ASSERT_TRUE(session.createSurface(113566,z,9));
    ASSERT_TRUE(session.createSurface(89555,z,10));
    ASSERT_TRUE(session.createSurface(160493,x,11));
    ASSERT_TRUE(session.createSurface(67246,z,12));
    ASSERT_TRUE(session.createSurface(113313,x,13));
    ASSERT_TRUE(session.createSurface(109139,z,14));
    ASSERT_TRUE(session.createSurface(189776,z,15));
    ASSERT_TRUE(session.createBoundary(174552,16));
    ASSERT_TRUE(session.createBoundary(190834,17));
    ASSERT_TRUE(session.createBoundary(105337,18));

    session.refineSurfaceSheet();

    ASSERT_TRUE(session.colorDual());

    session.createBlock();

    //Nb blocks in the final primal structure

    //Not working debug needed
    //ASSERT_EQ(hmesh.getNbRegions(),17);
}