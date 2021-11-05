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
	double Lx = 1.0;						// Longueur du domaine dans la direction x
	double Ly = 1.0;						// Longueur du domaine dans la direction y
	double dx = Lx/(Nx-1);				// Pas d'espace suivant l'axe x
	double dy = Ly/(Ny-1);				// Pas d'espace suivant l'axe y

	GridBuilder gb(&m,2);
	gb.execute(Nx, dx, Ny, dy);

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
	ASSERT_EQ(Nbr_BoundaryNodes, bnd_node_ids.size());			// Vérification : est-ce qu'on a tous les noeuds des bords

	Variable<int>* var_bnd = m.newVariable<int,GMDS_NODE>("constraint");
	for(auto id:bnd_node_ids){
		var_bnd->value(id)=1;								// Les noeuds des bords sont marqués par l'entier 1
	}

	// ------------ PERTURBATION ALEATOIRE DU MAILLAGE ------------
	// Initialisation de la variable n qui est un noeud et n_coords qui récupère les coordonnées du noeud
	Node n = m.get<Node>(1);
	math::Point n_coords = n.point();
	double c1 = 0.0;
	double c2 = 0.0;
	for(auto id:m.nodes()){
		if(var_bnd->value(id)==0){
			n = m.get<Node>(id);
		   n_coords = n.point();
			c1 = rand()/ double(RAND_MAX) ;
			c2 = rand()/ double(RAND_MAX) ;
			//std::cout << "Numéro node :" << n << std::endl ;
			//std::cout << "Nombre aléatoire c1 =" << c1 << std::endl ;
			n.setX( n_coords.X() + (2.0*c1-1.0) * dx);
			n.setY( n_coords.Y() + (2.0*c2-1.0) * dy);
		}
	}
	// ---------------------------------------------------------------

	IGMeshIOService ioService_geom_init(&m);
	VTKWriter writer_geom_init(&ioService_geom_init);
	writer_geom_init.setCellOptions(N|F);
	writer_geom_init.setDataOptions(N|F);
	writer_geom_init.write("smooth2D_Perturbation_1_init.vtk");

	Smooth2D smoother(&m,var_bnd, 20);
	Smooth2D::STATUS result = smoother.execute();

	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("smooth2D_Perturbation_1_result.vtk");

	ASSERT_EQ(Smooth2D::SUCCESS, result);
}



TEST(ClaireTestClass, testGrid2D_Quart_Cylindre)
{
	Mesh m(MeshModel(DIM3 | R | F | N | F2N | N2F));

	int Nx = 30;								// Nombre de noeuds dans la direction x
	int Ny = 30;								// Nombre de noeuds dans la direction y
	double Lx = 1.0;						// Longueur du domaine dans la direction x
	double Ly = 1.0;						// Longueur du domaine dans la direction y
	double dx = Lx/(Nx-1);				// Pas d'espace suivant l'axe x
	double dy = Ly/(Ny-1);				// Pas d'espace suivant l'axe y

	double R_int = 1.0;					// Rayon du cylindre
	double R_ext = 5.0;					// Rayon de la frontière extérieure

	GridBuilder gb(&m,2);
	gb.execute(Nx, dx, Ny, dy);

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
	ASSERT_EQ(Nbr_BoundaryNodes, bnd_node_ids.size());			// Vérification : est-ce qu'on a tous les noeuds des bords

	Variable<int>* var_bnd = m.newVariable<int,GMDS_NODE>("constraint");
	for(auto id:bnd_node_ids){
		var_bnd->value(id)=1;								// Les noeuds des bords sont marqués par l'entier 1
	}


	// ------------ REDEFINITION DES BORDS DU MAILLAGE ------------
	Node n = m.get<Node>(1);
	math::Point n_coords = n.point();
	double theta(0.0);
	// On parcours les noeuds du bord pour les replacer et former un domaine
	// quart de cylindre creux.
	// ATTENTION : ne fonctionne que grace à l'observation de la manière dont
	// sont numérotés les noeuds lors de la génération du maillage avec
	// GridBuilder
	for(auto id:bnd_node_ids){
			n = m.get<Node>(id);
			n_coords = n.point();
			//std::cout << "Numéro node :" << n << std::endl ;
			//std::cout << "Nombre aléatoire c1 =" << c1 << std::endl ;
			n.setX( n_coords.X() );
			n.setY( n_coords.Y() );
		   // "Bord gauche" du domaine carré
		   if (id <= Ny-1) {
			   n.setX( - R_int - (double(id)/(Ny-1))*(R_ext-R_int) );
				n.setY( 0.0 );
		   }
		   // "Bord droit" du domaine carré
		   else if (id >= (Ny-1)*Nx) {
			   n.setX( 0.0 );
			   n.setY( R_int + ((double(id)-(Nx-1.0)*Ny)/(Ny-1.0))*(R_ext-R_int) );
		   }
		   // "Bord bas et haut" du domaine carré
		   else {
			   // "Bord bas" du domaine carré
			   if ( id%Ny == 0 ) {
				   int k = id / Ny;
				   theta = M_PI - (double(k) / (Nx - 1.0)) * M_PI / 2.0;
				   n.setX(R_int*cos(theta) );
				   n.setY(R_int*sin(theta) );
			   }
			   // "Bord haut" du domaine carré
			   else {
				   int k = ( (id+1) / Ny) - 1;
				   theta = M_PI - (double(k) / (Nx - 1.0)) * M_PI / 2.0;
				   n.setX(R_ext*cos(theta) );
				   n.setY(R_ext*sin(theta) );
			   }
		   }
	}
	// ---------------------------------------------------------------

	IGMeshIOService ioService_geom_init(&m);
	VTKWriter writer_geom_init(&ioService_geom_init);
	writer_geom_init.setCellOptions(N|F);
	writer_geom_init.setDataOptions(N|F);
	writer_geom_init.write("smooth2D_test_quart_cylindre_init.vtk");

	Smooth2D smoother(&m,var_bnd, 1000);
	Smooth2D::STATUS result = smoother.execute();

	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("smooth2D_test_quart_cylindre_result.vtk");

	ASSERT_EQ(Smooth2D::SUCCESS, result);
}