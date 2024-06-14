//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/cad/GeomMeshLinker.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKReader.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(GeomLinkerTestSuite, fromSurfMesh)
{

    // WE WRITE
    Mesh m_vol(gmds::MeshModel(DIM3|R|F|E|N|
                               R2N|R2F|R2E|
                               F2N|F2R|F2E|
                               E2F|E2N|N2E|N2R));

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
    VTKWriter vtkW(&ioService);
    vtkW.setCellOptions(gmds::N|gmds::R);
    vtkW.setDataOptions(gmds::N|gmds::R);
    vtkW.write("toto_link.vtk");

	 VTKWriter vtkE(&ioService);
	 vtkE.setCellOptions(gmds::N|gmds::E);
	 vtkE.setDataOptions(gmds::N|gmds::E);
	 vtkE.write("toto_link_e.vtk");

	 VTKWriter vtkF(&ioService);
	 vtkF.setCellOptions(gmds::N|gmds::F);
	 vtkF.setDataOptions(gmds::N|gmds::F);
	 vtkF.write("toto_link_f.vtk");

    ASSERT_EQ(cad::GeomMeshLinker::LinkPoint, linker.getGeomDim<Node>(1));
    ASSERT_EQ(cad::GeomMeshLinker::LinkPoint, linker.getGeomDim<Node>(2));
    ASSERT_EQ(cad::GeomMeshLinker::LinkPoint, linker.getGeomDim<Node>(4));
    ASSERT_EQ(cad::GeomMeshLinker::LinkPoint, linker.getGeomDim<Node>(7));


    ASSERT_EQ(cad::GeomMeshLinker::LinkCurve, linker.getGeomDim<Node>(9));
    ASSERT_EQ(cad::GeomMeshLinker::LinkCurve, linker.getGeomDim<Node>(15));
    ASSERT_EQ(cad::GeomMeshLinker::LinkCurve, linker.getGeomDim<Node>(18));
    ASSERT_EQ(cad::GeomMeshLinker::LinkCurve, linker.getGeomDim<Node>(19));

    ASSERT_EQ(cad::GeomMeshLinker::LinkSurface, linker.getGeomDim<Node>(20));
    ASSERT_EQ(cad::GeomMeshLinker::LinkSurface, linker.getGeomDim<Node>(24));


    ASSERT_EQ(2, linker.getGeomId<Node>(1));
    ASSERT_EQ(3, linker.getGeomId<Node>(2));
    ASSERT_EQ(5, linker.getGeomId<Node>(4));
    ASSERT_EQ(8, linker.getGeomId<Node>(7));

    ASSERT_EQ(4, linker.getGeomId<Node>(9));
    ASSERT_EQ(9, linker.getGeomId<Node>(15));
    ASSERT_EQ(11, linker.getGeomId<Node>(18));
    ASSERT_EQ(2, linker.getGeomId<Node>(19));

    ASSERT_EQ(4, linker.getGeomId<Node>(20));
    ASSERT_EQ(3, linker.getGeomId<Node>(24));

    ASSERT_EQ(cad::GeomMeshLinker::LinkPoint  , linker.getGeomInfo<Node>(4). first);
    ASSERT_EQ(cad::GeomMeshLinker::LinkCurve  , linker.getGeomInfo<Node>(18).first);
    ASSERT_EQ(cad::GeomMeshLinker::LinkSurface, linker.getGeomInfo<Node>(24).first);
    ASSERT_EQ(5 , linker.getGeomInfo<Node>(4). second);
    ASSERT_EQ(11, linker.getGeomInfo<Node>(18).second);
    ASSERT_EQ(3 , linker.getGeomInfo<Node>(24).second);

    Node n1 = m_vol.get<Node>(1);
    ASSERT_EQ(cad::GeomMeshLinker::LinkPoint, linker.getGeomDim(n1));
    ASSERT_EQ(2, linker.getGeomId(n1));

}
/*----------------------------------------------------------------------------*/
TEST(GeomLinkerTestSuite, linkUsingTypes){
	 Mesh m_vol(gmds::MeshModel(DIM3|R|F|E|N|
	                            R2N|R2F|R2E|
	                            F2N|F2R|F2E|
	                            E2F|E2N|N2E|N2R));

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
	 cad::GeomMeshLinker linker(&m_vol,&manager);

	 linker.linkToPoint(m_vol.get<Node>(1),2);
	 linker.linkToPoint(m_vol.get<Node>(2),3);
	 linker.linkToPoint(m_vol.get<Node>(4),5);
	 linker.linkToPoint(m_vol.get<Node>(7),8);

	 ASSERT_EQ(cad::GeomMeshLinker::LinkPoint, linker.getGeomDim<Node>(1));
	 ASSERT_EQ(cad::GeomMeshLinker::LinkPoint, linker.getGeomDim<Node>(2));
	 ASSERT_EQ(cad::GeomMeshLinker::LinkPoint, linker.getGeomDim<Node>(4));
	 ASSERT_EQ(cad::GeomMeshLinker::LinkPoint, linker.getGeomDim<Node>(7));
	 ASSERT_EQ(2, linker.getGeomId<Node>(1));
	 ASSERT_EQ(3, linker.getGeomId<Node>(2));
	 ASSERT_EQ(5, linker.getGeomId<Node>(4));
	 ASSERT_EQ(8, linker.getGeomId<Node>(7));

	 linker.linkToCurve(m_vol.get<Edge>(14),2);
	 linker.linkToCurve(m_vol.get<Edge>(67),7);
	 linker.linkToCurve(m_vol.get<Edge>(83),6);
	 linker.linkToCurve(m_vol.get<Edge>(109),1);

	 ASSERT_EQ(cad::GeomMeshLinker::LinkCurve, linker.getGeomDim<Edge>(14));
	 ASSERT_EQ(cad::GeomMeshLinker::LinkCurve, linker.getGeomDim<Edge>(67));
	 ASSERT_EQ(cad::GeomMeshLinker::LinkCurve, linker.getGeomDim<Edge>(83));
	 ASSERT_EQ(cad::GeomMeshLinker::LinkCurve, linker.getGeomDim<Edge>(109));
	 ASSERT_EQ(2, linker.getGeomId<Edge>(14));
	 ASSERT_EQ(7, linker.getGeomId<Edge>(67));
	 ASSERT_EQ(6, linker.getGeomId<Edge>(83));
	 ASSERT_EQ(1, linker.getGeomId<Edge>(109));

	 linker.linkToSurface(m_vol.get<Face>(90), 6);
	 linker.linkToSurface(m_vol.get<Face>(145), 1);

	 ASSERT_EQ(cad::GeomMeshLinker::LinkSurface, linker.getGeomDim<Face>(90));
	 ASSERT_EQ(cad::GeomMeshLinker::LinkSurface, linker.getGeomDim<Face>(145));
	 ASSERT_EQ(6, linker.getGeomId<Face>(90));
	 ASSERT_EQ(1, linker.getGeomId<Face>(145));

	 linker.writeVTKDebugMesh("toto_type_link_debug.vtk");
}

