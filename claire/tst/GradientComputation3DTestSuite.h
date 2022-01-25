//
// Created by rochec on 24/01/2022.
//

#include <gmds/claire/LevelSet.h>
#include <gmds/claire/LevelSetFromIntToOut.h>
#include <gmds/claire/GradientComputation3D.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gtest/gtest.h>
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

TEST(GradientComputation3DTestClass, GradientComputation3D_Test1)
{
	Mesh m(MeshModel(DIM3 | R | F | E | N |
	                 R2N | F2N | E2N | R2F | F2R |
	                 F2E | E2F | R2E | N2R | N2F | N2E));
	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/Cube.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::R);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doctor(&m);
	doctor.buildFacesAndR2F();
	doctor.buildEdgesAndX2E();
	doctor.updateUpwardConnectivity();

	int markFrontNodes = m.newMark<gmds::Node>();

	// Initialisation de la marque pour noter quels fronts sont à avancer
	for(auto id:m.nodes()){
		Node n = m.get<Node>(id);
		double coord_x = n.X() ;
		if ( abs(coord_x) <= pow(10,-6) ) {
			// For this test case, the front to advance is the boundary where x=0
			m.mark<Node>(id,markFrontNodes);
		}
	}

	std::cout << "Fin de l'initialisation des marques" << std::endl ;

	LevelSetNaif ls(&m, markFrontNodes);
	LevelSetNaif::STATUS result_ls = ls.execute();

	ASSERT_EQ(LevelSetNaif::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	GradientComputation3D grad3D(&m, m.getVariable<double,GMDS_NODE>("distance"));
	GradientComputation3D::STATUS result = grad3D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("GradientComputation3D_Test1_Result.vtk");

	ASSERT_EQ(GradientComputation3D::SUCCESS, result);
}


TEST(GradientComputation3DTestClass, GradientComputation3D_Test2)
{
	Mesh m(MeshModel(DIM3 | R | F | E | N |
	                 R2N | F2N | E2N | R2F | F2R |
	                 F2E | E2F | R2E | N2R | N2F | N2E));
	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/B0.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::R);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doctor(&m);
	doctor.buildFacesAndR2F();
	doctor.buildEdgesAndX2E();
	doctor.updateUpwardConnectivity();

	int markFrontNodes = m.newMark<gmds::Node>();

	// Initialisation des marques sur les bords internes et externes pour noter quels fronts sont à avancer
	for(auto id:m.nodes()){
		Node n = m.get<Node>(id);
		double coord_y = n.Y() ;
		double coord_x = n.X() ;
		double rayon;
		rayon = sqrt( (pow(coord_x, 2) + pow(coord_y + 2.5, 2)) ) ;
		if ( (rayon - 2.5) < pow(10,-3)) {
			// For this test case, the front to advance is the boundary where x²+y²=2.5
			m.mark<Node>(id,markFrontNodes);
		}
	}

	std::cout << "Fin de l'initialisation des marques" << std::endl ;

	LevelSetNaif ls(&m, markFrontNodes);
	LevelSetNaif::STATUS result_ls = ls.execute();

	ASSERT_EQ(LevelSetNaif::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	GradientComputation3D grad3D(&m, m.getVariable<double,GMDS_NODE>("distance"));
	GradientComputation3D::STATUS result = grad3D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("GradientComputation3D_Test2_Result.vtk");

	ASSERT_EQ(GradientComputation3D::SUCCESS, result);
}
