//
// Created by ledouxf on 1/22/19.
//

/*----------------------------------------------------------------------------*/
#include <cstdlib>
#include <dirent.h>
#include <gtest/gtest.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/MeditReader.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/frame/CrossFieldGeneration2D.h>
#include <iostream>

#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
TEST(Frame2dTestClass, test1)
{
    // WE READ
    gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
				 gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/half_disk.vtk";

        gmds::IGMeshIOService ioService(&m);
        gmds::VTKReader vtkReader(&ioService);
        vtkReader.setCellOptions(gmds::N|gmds::F);
        vtkReader.read(vtk_file);
    
        gmds::MeshDoctor doc(&m);
    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    CrossFieldGeneration2D field_generator(&m);
    field_generator.setDebugPrefix("frame2d_debug_unit_test");
    field_generator.execute(CrossFieldGeneration2D::laplace_solve);

    gmds::VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::F);
    vtkWriter.setDataOptions(gmds::N|gmds::F);
    vtkWriter.write("frame2d_unit_test.vtk");

    ASSERT_TRUE(true);
}
