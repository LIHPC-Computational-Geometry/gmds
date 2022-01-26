//
// Created by rochec on 26/01/2022.
//

#include <gmds/claire/LevelSet.h>
#include <gmds/claire/LevelSetFromIntToOut.h>
#include <gmds/claire/GradientComputation2D.h>
#include <gmds/claire/PointFollowingVectorField2D.h>
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


TEST(PointFollowingVectorField2DTestClass, PointFollowingVectorField2D_Test1)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
	                             gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/Carre_maxsize_0.01.vtk";

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

	int markFrontNodes = m.newMark<gmds::Node>();
	for(auto id:bnd_node_ids){
		Node n = m.get<Node>(id);
		double coord_y = n.Y() ;
		if (coord_y == 0) {
			// For the square test case, the front to advance is the boundary where y=0
			m.mark<Node>(id,markFrontNodes);
		}
	}

	LevelSetNaif ls(&m, markFrontNodes);
	LevelSetNaif::STATUS result_ls = ls.execute();

	ASSERT_EQ(LevelSetNaif::SUCCESS, result_ls);

	m.unmarkAll<Node>(markFrontNodes);
	m.freeMark<Node>(markFrontNodes);

	GradientComputation2D grad2D(&m, m.getVariable<double,GMDS_NODE>("distance"));
	GradientComputation2D::STATUS result_grad = grad2D.execute();
	ASSERT_EQ(GradientComputation2D::SUCCESS, result_grad);

	math::Point M(0.5, 0.5, 0.0);
	double distance = 0.3;
	PointFollowingVectorField2D pfvf2D(&m, M, distance, m.getVariable<math::Vector3d ,GMDS_FACE>("gradient_2D"));
	PointFollowingVectorField2D::STATUS result = pfvf2D.execute();

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