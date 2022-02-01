//
// Created by rochec on 26/01/2022.
//

#include <gmds/claire/LevelSet.h>
#include <gmds/claire/LevelSetFromIntToOut.h>
#include <gmds/claire/GradientComputation2D.h>
#include <gmds/claire/PointFollowingVectorField2D.h>
#include <gmds/claire/PointFollowingVectorField3D.h>
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

/* CAS TEST CLASSE PointFollowingVectorField2D */

TEST(PointFollowingVectorField2DTestClass, PointFollowingVectorField2D_Test1)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
	                             gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	//std::string vtk_file = dir+"/Carre.vtk";
	std::string vtk_file = dir+"/Carre_maxsize_0.01.vtk";		// Same geo but refined mesh

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
	doc.orient2DFaces();

	//Get the boundary node ids
	BoundaryOperator2D bnd_op(&m);
	std::vector<TCellID> bnd_node_ids;
	bnd_op.getBoundaryNodes(bnd_node_ids);

	// Initialisation des marques sur le front à avancer
	int markFrontNodes = m.newMark<gmds::Node>();
	for(auto id:bnd_node_ids){
		Node n = m.get<Node>(id);
		double coord_y = n.Y() ;
		if (coord_y == 0) {
			// For the square test case, the front to advance is the boundary where y=0
			m.mark<Node>(id,markFrontNodes);
		}
	}

	// Calcul des Level Set
	LevelSetNaif ls(&m, markFrontNodes);
	LevelSetNaif::STATUS result_ls = ls.execute();
	ASSERT_EQ(LevelSetNaif::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	// Calcul du gradient du champ de Level Set
	GradientComputation2D grad2D(&m, m.getVariable<double,GMDS_NODE>("distance"));
	GradientComputation2D::STATUS result_grad = grad2D.execute();
	ASSERT_EQ(GradientComputation2D::SUCCESS, result_grad);

	// Placement du point P à la distance souhaitée suivant le champ de gradient
	math::Point M(0.5, 0.5, 0.0);
	double distance = 0.3;
	PointFollowingVectorField2D pfvf2D(&m, M, distance, m.getVariable<math::Vector3d ,GMDS_FACE>("gradient_2D"));
	PointFollowingVectorField2D::STATUS result = pfvf2D.execute();

	// Placement du point P à la distance souhaitée suivant le champ de gradient
	M.setXYZ(0.3, 0.3, 0.0);
	distance = 0.5;
	pfvf2D = PointFollowingVectorField2D(&m, M, distance, m.getVariable<math::Vector3d ,GMDS_FACE>("gradient_2D"));
	result = pfvf2D.execute();

	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("PointFollowingVectorField2D_Test1_Result.vtk");

	ASSERT_EQ(PointFollowingVectorField2D::SUCCESS, result);
}



/* CAS TEST CLASSE PointFollowingVectorField3D */

TEST(PointFollowingVectorField3DTestClass, PointFollowingVectorField3D_Test1)
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
		double coord_y = n.Y() ;
		if ( abs(coord_y) <= pow(10,-6) ) {
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
	GradientComputation3D::STATUS result_grad = grad3D.execute();
	ASSERT_EQ(GradientComputation3D::SUCCESS, result_grad);

	// Placement du point P à la distance souhaitée suivant le champ de gradient
	math::Point M(0.1, 0.0, 0.3);
	double distance = 0.8;
	PointFollowingVectorField3D pfvf3D(&m, M, distance, m.getVariable<math::Vector3d ,GMDS_REGION>("gradient_3D"));
	PointFollowingVectorField3D::STATUS result = pfvf3D.execute();

	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("PointFollowingVectorField3D_Test1_Result.vtk");

	ASSERT_EQ(PointFollowingVectorField3D::SUCCESS, result);
}