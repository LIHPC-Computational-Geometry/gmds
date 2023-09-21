/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <filesystem>
/*----------------------------------------------------------------------------*/
#include "gmds/io/IGMeshIOService.h"
#include "gmds/io/VTKWriter.h"
#include <gmds/blocking/Blocking.h>
#include <gmds/blocking/SheetCollapse.h>
#include <gmds/blocking/WriterDartsVTK.h>
/*----------------------------------------------------------------------------*/
void
setUpFrom(gmds::cad::FACManager &AGeomManager, const std::string& AFileName)
{
	gmds::Mesh m_vol(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::F | gmds::E | gmds::N |
	                                 gmds::R2N | gmds::R2F | gmds::R2E | gmds::F2N |
	                                 gmds::F2R | gmds::F2E
	                                 | gmds::E2F | gmds::E2N | gmds::N2E));
	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir + "/" + AFileName;
	gmds::IGMeshIOService ioService(&m_vol);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::R);
	vtkReader.read(vtk_file);
	gmds::MeshDoctor doc(&m_vol);
	doc.buildFacesAndR2F();
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	AGeomManager.initFrom3DMesh(&m_vol);
}
/*----------------------------------------------------------------------------*/
TEST(PolyCubePipelineTestSuite, simple_read)
{
	gmds::Mesh m_vol(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::N | gmds::R2N));
	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir + "/B0_blocking.vtk";
	gmds::IGMeshIOService ioService(&m_vol);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::R);
	vtkReader.read(vtk_file);
	std::cout<<"Nodes: "<<m_vol.getNbNodes()<<std::endl;
	std::cout<<"Blocks: "<<m_vol.getNbRegions()<<std::endl;
	//	ASSERT_EQ(bl3d.nbBlocks(), 24);
}
/*----------------------------------------------------------------------------*/
