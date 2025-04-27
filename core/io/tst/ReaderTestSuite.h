/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/MeditReader.h>
#include <gmds/io/IGMeshIOService.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
TEST(MeshClass, testIOVTK_NF)
{
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N));

    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/triangulated_quad.vtk";

    gmds::IGMeshIOService ioService(&m);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::F);
    vtkReader.read(vtk_file);

    ASSERT_EQ(m.getNbNodes(),29);
    ASSERT_EQ(m.getNbFaces(),40);

}
/*----------------------------------------------------------------------------*/
TEST(MeshClass, testIOVTK_NR)
{
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::R|gmds::N));

    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/tet_in_box.vtk";

    gmds::IGMeshIOService ioService(&m);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::R);
    vtkReader.read(vtk_file);

    ASSERT_EQ(m.getNbNodes(),29);
    ASSERT_EQ(m.getNbRegions(),61);
}
/*----------------------------------------------------------------------------*/
TEST(MeshClass, testIOVTK_NEF)
{
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::E|gmds::N));

    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/triangulated_quad.vtk";

    gmds::IGMeshIOService ioService(&m);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::E|gmds::F);
    vtkReader.read(vtk_file);

    ASSERT_EQ(m.getNbNodes(),29);
    ASSERT_EQ(m.getNbEdges(),16);
    ASSERT_EQ(m.getNbFaces(),40);

}
/*----------------------------------------------------------------------------*/
TEST(MeshClass, testIOMedit)
{
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::F2N));

    std::string dir(TEST_SAMPLES_DIR);
    std::string med_file = dir+"/triangulated_quad.mesh";

    gmds::IGMeshIOService ioService(&m);

    gmds::MeditReader meditReader(&ioService);

    meditReader.setCellOptions(gmds::N|gmds::F);

    meditReader.read(med_file);
    ASSERT_EQ(m.getNbNodes(),29);
    ASSERT_EQ(m.getNbFaces(),40);

}
/*----------------------------------------------------------------------------*/
TEST(MeshClass, testIOVTK_several_types)
{
	 gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::R | gmds::N | gmds::R2N));

	 std::string dir(TEST_SAMPLES_DIR);
	 std::string vtk_file = dir+"/several_types.vtk";

	 gmds::IGMeshIOService ioService(&m);
	 gmds::VTKReader vtkReader(&ioService);
	 vtkReader.setCellOptions(gmds::N|gmds::R);
	 vtkReader.read(vtk_file);

	 ASSERT_EQ(m.getNbNodes(),14);
	 ASSERT_EQ(m.getNbRegions(),4);
	 ASSERT_EQ(m.getNbTetrahedra(),1);
	 ASSERT_EQ(m.getNbHexahedra(),2);
}
/*----------------------------------------------------------------------------*/