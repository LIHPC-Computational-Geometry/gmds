//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#include <gmds/claire/Smooth2D.h>
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
TEST(ClaireTestClass, testGrid2D_1)
{
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

	Smooth2D smoother(&m,var_bnd);
	Smooth2D::STATUS result = smoother.execute();


	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("smooth2D_result.vtk");

	ASSERT_EQ(Smooth2D::SUCCESS, result);
}
