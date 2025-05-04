/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gmds/igalgo/SurfaceReorient.h>
#include <iostream>
#include <gmds/io/IGMeshIOService.h>
#include <unit_test_config.h>
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(SurfaceReorientTestClass, testReorient)
{
    // WE WRITE
    Mesh m(MeshModel(DIM3|F|N|F2N|N2F));

    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/reorient3D_test.vtk";
    gmds::IGMeshIOService ioService(&m);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::N|gmds::F);
    vtkReader.read(vtk_file);

    gmds::MeshDoctor doc(&m);
    doc.updateUpwardConnectivity();
    gmds::SurfaceReorient reorient(&m);

    ASSERT_TRUE(reorient.isValid());

	 auto nb_reorient = reorient.execute();
    ASSERT_TRUE(nb_reorient!=0);

	 gmds::VTKWriter w(&ioService);
	 w.setCellOptions(gmds::N|gmds::F);
	 w.write("toto.vtk");


}