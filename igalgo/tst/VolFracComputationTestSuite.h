/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/VTKReader.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/igalgo/VolFracComputation.h>
#include <iostream>
#include <gmds/io/IGMeshIOService.h>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
TEST(VolFracComputation, DISABLED_test2D)
{
	gmds::Mesh mImprint(MeshModel(DIM2|F|N|F2N));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/triangulated_quad.vtk";
	gmds::IGMeshIOService ioService(&mImprint);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(vtk_file);

	gmds::Mesh m(MeshModel(DIM2|F|N|F2N));
	gmds::GridBuilder gb(&m,2);
	gb.execute(3,.5, 4, .5);

	gmds::Variable<double>* vf = m.newVariable<double, gmds::GMDS_FACE>( "vf");

	gmds::volfraccomputation_2d(&m, &mImprint, vf);

	gmds::IGMeshIOService ioService_write(&m);
	gmds::VTKWriter vtkWriter(&ioService_write);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("VolFracComputation_test2D.vtk");
}
