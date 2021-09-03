//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/igalgo/BoundaryExtractor3D.h>
#include <gmds/cad/FACManager.h>
#include <gmds/cad/FACSurface.h>
#include <gmds/cad/FACVolume.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/IGMeshIOService.h>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(FacManagerTestSuite, fromSurfMesh)
{
    // WE WRITE
    gmds::Mesh m_vol(gmds::MeshModel(DIM3|R|F|E|N|
                                     R2N|R2F|R2E|
                                     F2N|F2R|F2E|
                                     E2F|E2N|N2E));
    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/tet_in_box.vtk";
    
    gmds::IGMeshIOService ioService(&m_vol);

    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::R);
    vtkReader.read(vtk_file);

    gmds::MeshDoctor doc(&m_vol);
    doc.buildFacesAndR2F();
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();


    cad::FACManager manager;

    manager.initFrom3DMesh(&m_vol);


    ASSERT_EQ(manager.getNbPoints(),8);
    ASSERT_EQ(manager.getNbSurfaces(),6);
    ASSERT_EQ(manager.getNbCurves(),12);

    std::vector<cad::GeomSurface*> surfs;
    manager.getSurfaces(surfs);

    for(auto s:surfs){
        ASSERT_NEAR(s->computeArea(),100, 0.1);
        math::Point c(0,0,0);
        ASSERT_EQ(c.distance(s->closestPoint(c)),5);
    }

    std::vector<cad::GeomCurve*> curvs;
    manager.getCurves(curvs);
    for(auto c:curvs){
        ASSERT_EQ(c->length(),10);
        math::Point p(0,0,0);
        ASSERT_NEAR(p.distance(c->closestPoint(p)),7.07107,0.01);
    }

    std::vector<cad::GeomPoint*> pnts;
    manager.getPoints(pnts);
    for(auto p:pnts){
        math::Point c(0,0,0);
        ASSERT_NEAR(c.distance(p->point()),8.66,0.01);
    }
}
