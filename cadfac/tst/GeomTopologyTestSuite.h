//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/cad/GeomMeshLinker.h>
#include <gmds/igalgo/GridBuilder.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(GeomTopologyTestSuite, cube_topo)
{
    Mesh m_vol(gmds::MeshModel(DIM3|R|F|E|N|
                               R2N|R2F|R2E|
                               F2N|F2R|F2E|
                               E2F|E2N|N2E|N2R|N2F));

    GridBuilder gb(&m_vol, 3);
    gb.execute(3,1.,3,1.,3,1.);
    ASSERT_EQ(m_vol.getNbRegions(),8);
    MeshDoctor doc(&m_vol);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    cad::FACManager manager;
    cad::GeomMeshLinker linker;

    manager.initAndLinkFrom3DMesh(&m_vol,&linker);

//==================================
// tests on the points
    std::vector<cad::GeomPoint*> pnts;
    manager.getPoints(pnts);
    for(auto p:pnts){
        ASSERT_EQ(p->volumes().size(),1);
        ASSERT_EQ(p->surfaces().size(),3);
        ASSERT_EQ(p->curves().size(),3);
    }
//==================================
// tests on the curves
    std::vector<cad::GeomCurve*> curves;
    manager.getCurves(curves);
    for(auto c:curves){
        ASSERT_EQ(c->volumes().size(),1);
        std::vector<cad::GeomSurface*> surfs =c->surfaces();
        std::vector<cad::GeomPoint*> pts =c->points();
        ASSERT_EQ(c->surfaces().size(),2);
        ASSERT_EQ(c->points().size(),2);
    }
//==================================
// tests on the surfaces
    std::vector<cad::GeomSurface*> surfs;
    manager.getSurfaces(surfs);
    for(auto s:surfs){
        ASSERT_EQ(s->volumes().size(),1);
        ASSERT_EQ(s->curves().size(),4);
        ASSERT_EQ(s->points().size(),4);
    }
//==================================
// tests on the volumes
    std::vector<cad::GeomVolume*> vols;
    manager.getVolumes(vols);
    for(auto v:vols){
        ASSERT_EQ(v->surfaces().size(),6);
        ASSERT_EQ(v->curves().size(),12);
        ASSERT_EQ(v->points().size(),8);
    }
}
/*----------------------------------------------------------------------------*/
TEST(GeomTopologyTestSuite, cube_convex_curve)
{
    Mesh m_vol(gmds::MeshModel(DIM3|R|F|E|N|
                               R2N|R2F|R2E|
                               F2N|F2R|F2E|
                               E2F|E2N|N2E|N2R|N2F));

    GridBuilder gb(&m_vol, 3);
    gb.execute(3,1.,3,1.,3,1.);
    ASSERT_EQ(m_vol.getNbRegions(),8);
    MeshDoctor doc(&m_vol);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    cad::FACManager manager;
    cad::GeomMeshLinker linker;

    manager.initAndLinkFrom3DMesh(&m_vol,&linker);

    manager.write_surfaces("fac_model.vtk");
    // Check that all the cuvrve are strictly convex
    std::vector<cad::GeomCurve*> curves;
    manager.getCurves(curves);
    for(auto c:curves){
        ASSERT_EQ(c->getCurvatureInfo(),
                  cad::GeomCurve::Convex);
    }
}
/*----------------------------------------------------------------------------*/
TEST(GeomTopologyTestSuite, cube_convex_B45)
{
    Mesh m_vol(gmds::MeshModel(DIM3|R|F|E|N|
                               R2N|R2F|R2E|
                               F2N|F2R|F2E|
                               E2F|E2N|N2E|N2R|N2F));
    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/B45.vtk";
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

    manager.write_surfaces("fac_model.vtk");
    // Check that all the cuvrve are strictly convex
    std::vector<cad::GeomCurve*> curves;
    manager.getCurves(curves);
    int nb_convex =0;
    int nb_concave=0;
    for(auto c:curves){
        cad::GeomCurve::CurvatureInfo ci = c->getCurvatureInfo();
        if(ci==cad::GeomCurve::Convex)
            nb_convex++;
        else if(ci==cad::GeomCurve::Concave)
            nb_concave++;
    }
    ASSERT_EQ(17, nb_convex);
    ASSERT_EQ(1 , nb_concave);
}/*----------------------------------------------------------------------------*/
TEST(GeomTopologyTestSuite, case2D)
{
    gmds::Mesh m_surf(gmds::MeshModel(DIM3|F|E|N|F2N|N2F|E2N));
    gmds::Mesh m_bnd (gmds::MeshModel(DIM3|E|N|N2E|E2N));

    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/HolesInSquare0.vtk";
    gmds::IGMeshIOService ioService(&m_surf);

    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::F);
    vtkReader.read(vtk_file);

    gmds::MeshDoctor doc(&m_surf);
    doc.updateUpwardConnectivity();
    doc.buildBoundaryCells();
    cad::FACManager manager;
    cad::GeomMeshLinker linker;

    manager.initAndLinkFrom2DMesh(&m_surf,&linker);

    // Check curve --> point connectivity
    std::vector<cad::GeomCurve*> curves;
    manager.getCurves(curves);
    int nb_pts_2 =0;
    int nb_pts_0=0;
    for(auto c:curves){
        std::vector<cad::GeomPoint*> pts = c->points();
        if(pts.size()==0)
            nb_pts_0++;
        else if(pts.size()==2)
            nb_pts_2++;
    }
    ASSERT_EQ(1, nb_pts_0);
    ASSERT_EQ(5, nb_pts_2);
    // Check point --> curve connectivity
    std::vector<cad::GeomPoint*> points;
    manager.getPoints(points);
    for(auto p:points){
        ASSERT_EQ(p->curves().size(),2);
    }
}