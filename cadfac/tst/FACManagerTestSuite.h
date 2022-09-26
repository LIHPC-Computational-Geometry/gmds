//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/igalgo/BoundaryExtractor3D.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/cadfac/FACSurface.h>
#include <gmds/cadfac/FACVolume.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
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

		  math::Point p_far(100000000,0,0);
		  math::Point p_on = c->closestPoint(p_far);
		  ASSERT_LT(p_on.distance(math::Point (0,0,0)), 10.);

		  c->project(p_far);
		  ASSERT_LT(p_far.distance(math::Point (0,0,0)), 10.);
		  ASSERT_DOUBLE_EQ(p_on.distance(p_far), 0.);
    }

    std::vector<cad::GeomPoint*> pnts;
    manager.getPoints(pnts);
    for(auto p:pnts){
        math::Point c(0,0,0);
        ASSERT_NEAR(c.distance(p->point()),8.66,0.01);
    }
}
/*----------------------------------------------------------------------------*/
TEST(FacManagerTestSuite, surf_projection)
{
    // WE WRITE
    gmds::Mesh m_vol(gmds::MeshModel(DIM3|R|F|E|N|
                                     R2N|R2F|R2E|
                                     F2N|F2R|F2E|
                                     E2F|E2N|N2E));
    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/S24_CAD_test.vtk";

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

    cad::GeomSurface* s1 = manager.getSurface(1);

    TCoord min_s1[3];
    TCoord max_s1[3];
    s1->computeBoundingBox(min_s1,max_s1);

    math::Point p(-10.6674, 24.1967, 109.342);

    s1->project(p);


}
/*----------------------------------------------------------------------------*/
TEST(FacManagerTestSuite, project)
{
	gmds::Mesh m_vol(gmds::MeshModel(DIM3|R|F|E|N|
	                                 R2N|R2F|R2E|
	                                 F2N|F2R|F2E|
	                                 E2F|E2N|N2E));
	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/simpleCube.vtk";

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

	cad::GeomSurface* s = manager.getSurface(2);

	TCoord min_s[3];
	TCoord max_s[3];
	s->computeBoundingBox(min_s,max_s);

	//math::Point p0(min_s1[0], min_s1[1], min_s1[2]);
	//math::Point p1(max_s1[0], max_s1[1], max_s1[2]);
	//std::cout <<p0<<" "<<p1<<std::endl;

	math::Point ps(0.1, -0.1, 0.1);
	math::Point ps_on(0.1, 0., 0.1);

	s->project(ps);
	//std::cout<<ps<<std::endl;
	ASSERT_NEAR(ps.distance(ps_on), 0., 10e-15);

	cad::GeomCurve* c = manager.getCurve(1);

	TCoord min_c[3];
	TCoord max_c[3];
	c->computeBoundingBox(min_c,max_c);

//	math::Point p0(min_c[0], min_c[1], min_c[2]);
//	math::Point p1(max_c[0], max_c[1], max_c[2]);
//	std::cout <<p0<<" "<<p1<<std::endl;

	math::Point pc(-0.1, 0.2, 3);
	math::Point pc_on(0., 0.2, 0.);

	c->project(pc);
	//std::cout<<ps<<std::endl;
	ASSERT_NEAR(pc.distance(pc_on), 0., 10e-15);
}
/*----------------------------------------------------------------------------*/
TEST(FacManagerTestSuite, automatic_blocks_classification)
{
	gmds::Mesh m(gmds::MeshModel(gmds::DIM2|gmds::F|gmds::E|gmds::N|gmds::N2E|gmds::E2N|gmds::E2F|gmds::N2F|gmds::F2N));

	std::string vtk_file = "/home/calderans/dev/v3.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::E|gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.updateUpwardConnectivity();

	cad::FACManager manager;
	cad::GeomMeshLinker linkerTri;
	manager.initAndLinkFrom2DMesh(&m, &linkerTri);

	Mesh blocks(gmds::MeshModel(gmds::DIM2|gmds::F|gmds::E|gmds::N|gmds::N2E|gmds::E2N|gmds::E2F|gmds::N2F|gmds::F2N));

	

	cad::GeomMeshLinker linker;
	linker.setGeometry(&manager);
	linker.setMesh(&blocks);
}