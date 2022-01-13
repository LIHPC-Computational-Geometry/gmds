//
// Created by rochec on 13/01/2022.
//

#ifndef GMDS_LEVELSET2DTESTSUITE_H
#define GMDS_LEVELSET2DTESTSUITE_H

#include <gmds/claire/LevelSet2D.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gtest/gtest.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(LevelSet2DTestClass, LevelSet2D_Test1)
{
	// Test
	Mesh m(MeshModel(DIM3 | R | F | N | F2N | N2F));

	GridBuilder gb(&m,2);
	gb.execute(3,1.0, 4, 1.0);

	ASSERT_EQ(m.getNbNodes(),12);
	ASSERT_EQ(m.getNbFaces(),6);

	//==================================================================
	// MESH PREPARATION
	//==================================================================
	MeshDoctor doctor(&m);
	doctor.updateUpwardConnectivity();

	//Get the boundary node ids
	BoundaryOperator2D bnd_op(&m);
	std::vector<TCellID> bnd_node_ids;
	bnd_op.getBoundaryNodes(bnd_node_ids);
	ASSERT_EQ(10, bnd_node_ids.size());

	Variable<int>* var_bnd = m.newVariable<int,GMDS_NODE>("constraint");
	for(auto id:bnd_node_ids){
		var_bnd->value(id)=1;
	}

	Node n5 = m.get<Node>(5);
	math::Point p5 = n5.point();
	math::Point p5_new(p5.X()+0.2,p5.Y()+0.3,p5.Z());
	n5.setPoint(p5_new);

	LevelSet2D ls(&m);
	LevelSet2D::STATUS result = ls.execute();

	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("LevelSet2D_Test1_Result.vtk");

	ASSERT_EQ(LevelSet2D::SUCCESS, result);
}

#endif     // GMDS_LEVELSET2DTESTSUITE_H
