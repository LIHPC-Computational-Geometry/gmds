/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/IGMeshIOService.h>
#include <unit_test_config.h>

#include<iostream>
#include <random>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(FacManagerFastRetrievalTestSuite, fromSurfMesh)
{
    // WE WRITE
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
	 manager.buildGTSTree(&m_vol);

	 std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-0.1, 1.1);

	 for(int i=0; i<999; i++) {
		 double x = dis(rd);
		 double y = dis(rd);
		 double z = dis(rd);
		 math::Point p(x, y, z);
		 if(x<0. || x>1. || y<0. || y>1. || z<0. || z>1.) {
			 ASSERT_FALSE(manager.is_in(p));
		 } else {
			 ASSERT_TRUE(manager.is_in(p));
		 }
	 }
}