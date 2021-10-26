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

	Node n5 = m.get<Node>(5);
	math::Point p5 = n5.point();
	math::Point p5_new(p5.X()+0.2,p5.Y()+0.3,p5.Z());
	n5.setPoint(p5_new);

	Smooth2D smoother(&m,var_bnd);
	Smooth2D::STATUS result = smoother.execute();


	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("smooth2D_result.vtk");

	ASSERT_EQ(Smooth2D::SUCCESS, result);
}



TEST(ClaireTestClass, testGrid2D_Perturbation_1)
{
	Mesh m(MeshModel(DIM3 | R | F | N | F2N | N2F));

	int Nx = 6;								// Nombre de noeuds dans la direction x
	int Ny = 6;								// Nombre de noeuds dans la direction y
	double dx = 1.0/Nx;					// Pas d'espace suivant l'axe x
	double dy = 1.0/Ny;					// Pas d'espace suivant l'axe y

	GridBuilder gb(&m,2);
	gb.execute(Nx, 1.0, Ny, 1.0);

	ASSERT_EQ(m.getNbNodes(),Nx*Ny);
	ASSERT_EQ(m.getNbFaces(),(Nx-1)*(Ny-1));

	//==================================================================
	// MESH PREPARATION
	//==================================================================
	MeshDoctor doctor(&m);
	doctor.updateUpwardConnectivity();

	//Get the boundary node ids
	BoundaryOperator2D bnd_op(&m);
	std::vector<TCellID> bnd_node_ids;
	bnd_op.getBoundaryNodes(bnd_node_ids);
	int Nbr_BoundaryNodes = 2*Nx + 2*(Ny-2);								// Nombre de noeuds sur les bords
	ASSERT_EQ(Nbr_BoundaryNodes, bnd_node_ids.size());

	Variable<int>* var_bnd = m.newVariable<int,GMDS_NODE>("constraint");
	for(auto id:bnd_node_ids){
		var_bnd->value(id)=1;								// Les noeuds des bords sont marqués par l'entier 1
	}

	Node n = m.get<Node>(1);
	math::Point n_coords = n.point();
	// Pb : sur quoi on itère ?
	for(auto id:var_bnd->value()){
		if(var_bnd->value(id)==0){
			n = m.get<Node>(id);
			std::cout << "Numéro node :" << n << std::endl ;
			n.setX( n_coords.X() + 0.5);
			n.setY( n_coords.Y() + 0.5);
		}
	}

	Smooth2D smoother(&m,var_bnd);
	Smooth2D::STATUS result = smoother.execute();


	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("smooth2D_test_2_result.vtk");

	ASSERT_EQ(Smooth2D::SUCCESS, result);
}