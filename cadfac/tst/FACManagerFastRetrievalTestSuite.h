/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/IGMeshIOService.h>
#include <unit_test_config.h>
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
	 manager.buildGTSTree();

    math::Point p0(10.01, -1.2, 80.);
	 math::Point p1( 9.99, -1.2, 80.);
	 ASSERT_FALSE(manager.is_in(p0));
	 ASSERT_TRUE(manager.is_in(p1));
}