//
// Created by ledouxf on 1/22/19.
//
/*----------------------------------------------------------------------------*/
#include <gmds/claire/Smooth2D.h>
#include <gmds/claire/Grid_Smooth2D.h>
#include <gmds/claire/ResLU.h>
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

	Smooth2D smoother(&m,var_bnd,1);
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


TEST(ClaireTestClass, testGrid2D_Plaque_Plane)
{
	Mesh m(MeshModel(DIM3 | R | F | N | F2N | N2F));

	int Nx = 50;								// Nombre de noeuds dans la direction x
	int Ny = 50;								// Nombre de noeuds dans la direction y
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


	// ------------ REDEFINITION DES BORDS DU MAILLAGE ------------
	Node n = m.get<Node>(1);
	math::Point n_coords = n.point();
	// On parcours les noeuds du bord pour les replacer et former un domaine
	// plaque plane.
	// ATTENTION : ne fonctionne que grace à l'observation de la manière dont
	// sont numérotés les noeuds lors de la génération du maillage avec
	// GridBuilder
	for(auto id:bnd_node_ids){
		n = m.get<Node>(id);
		n_coords = n.point();
		// "Bord gauche" du domaine carré
		if (id <= Ny-1) {
			n.setX( 0.0 );
			n.setY( double(id)/(Ny-1.0)*200.0 );
		}
		// "Bord droit" du domaine carré
		else if (id >= (Ny-1)*Nx) {
			n.setX( 1000.0 );
			n.setY( ((double(id)-(Nx-1.0)*Ny)/(Ny-1.0))*400.0 );
		}
		// "Bord bas et haut" du domaine carré
		else {
			// "Bord bas" du domaine carré
			if ( id%Ny == 0 ) {
				int k = id / Ny;
				n.setX( (double(k) / (Nx - 1.0))*1000.0 );
				n.setY(0.0 );
			}
			// "Bord haut" du domaine carré
			else {
				int k = ( (id+1) / Ny) - 1;
				n.setX( (double(k) / (Nx - 1.0))*1000.0 );
				n_coords = n.point();
				n.setY(0.2*n_coords.X() + 200.0 );
			}
		}
	}
	// ---------------------------------------------------------------

	IGMeshIOService ioService_geom_init(&m);
	VTKWriter writer_geom_init(&ioService_geom_init);
	writer_geom_init.setCellOptions(N|F);
	writer_geom_init.setDataOptions(N|F);
	writer_geom_init.write("smooth2D_test_plaque_plane_init.vtk");

	Smooth2D smoother(&m,var_bnd, 2000);
	Smooth2D::STATUS result = smoother.execute();

	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("smooth2D_test_plaque_plane_result.vtk");

	ASSERT_EQ(Smooth2D::SUCCESS, result);
}


TEST(ClaireTestClass, testUnstructuredMesh)
{
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::F2N | gmds::N2F));

	gmds::Node n0 = m.newNode(0,0,0);
	gmds::Node n1 = m.newNode(1,0,0);
	gmds::Node n2 = m.newNode(0,1,0);
	gmds::Node n3 = m.newNode(0.2,0.2,0);

	m.newTriangle(n0,n1,n3);
	m.newTriangle(n0,n3,n2);
	m.newTriangle(n1,n2,n3);

	gmds::MeshDoctor doc(&m);
	doc.updateUpwardConnectivity();

	//Get the boundary node ids
	BoundaryOperator2D bnd_op(&m);
	std::vector<TCellID> bnd_node_ids;
	bnd_op.getBoundaryNodes(bnd_node_ids);

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
	writer_geom.write("smooth2D_UnstructuredMesh.vtk");

	ASSERT_EQ(Smooth2D::SUCCESS, result);
}


TEST(ClaireTestClass, testUnstructuredMesh_2)
{
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::F2N | gmds::N2F));

	gmds::Node n0 = m.newNode(0,0,0);
	gmds::Node n1 = m.newNode(1,0,0);
	gmds::Node n2 = m.newNode(0,1,0);
	gmds::Node n3 = m.newNode(0.2,0.2,0);
	gmds::Node n4 = m.newNode(1, 1, 0);

	m.newTriangle(n0,n1,n3);
	m.newTriangle(n0,n3,n2);
	m.newTriangle(n1,n3,n4);
	m.newTriangle(n2,n3,n4);

	gmds::MeshDoctor doc(&m);
	doc.updateUpwardConnectivity();

	//Get the boundary node ids
	BoundaryOperator2D bnd_op(&m);
	std::vector<TCellID> bnd_node_ids;
	bnd_op.getBoundaryNodes(bnd_node_ids);

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
	writer_geom.write("smooth2D_UnstructuredMesh_2.vtk");

	ASSERT_EQ(Smooth2D::SUCCESS, result);
}




TEST(ClaireTestClass, testGrid_Smooth2D_1)
{
	Blocking2D m;
	Node n1 = m.newBlockCorner(0,0);
	Node n2 = m.newBlockCorner(1,0);
	Node n3 = m.newBlockCorner(1,1);
	Node n4=  m.newBlockCorner(0,1);

	Blocking2D::Block b1 = m.newBlock(n1,n2,n3,n4);

	Node n5 = m.newBlockCorner(2,0,0);
	Node n6 = m.newBlockCorner(2,1.5,0);
	Blocking2D::Block b2 = m.newBlock(n2,n5,n6,n3);
	b1.seNbDiscretizationI(10);
	b1.seNbDiscretizationJ(10);
	b2.seNbDiscretizationI(10);
	b2.seNbDiscretizationJ(10);

	m.initializeGridPoints();

	b1 = m.block(0);
	int Nx = b1.getNbDiscretizationI();
	int Ny = b1.getNbDiscretizationJ();

	// Perturbation of the mesh
	// Boucle sur les noeuds internes du bloc b0
	for (int i=1; i<Nx-1; i++) {
		for (int j=1; j<Ny-1; j++) {
			b1(i,j).setX(0.0);
			b1(i,j).setY(0.0);
		}
	}

	b1 = m.block(1);
	Nx = b1.getNbDiscretizationI();
	Ny = b1.getNbDiscretizationJ();

	// Perturbation of the mesh
	// Boucle sur les noeuds internes du bloc b1
	for (int i=1; i<Nx-1; i++) {
		for (int j=1; j<Ny-1; j++) {
			b1(i,j).setX(0.0);
			b1(i,j).setY(0.0);
		}
	}

	IGMeshIOService ios(&m);
	VTKWriter writer(&ios);
	writer.setCellOptions(N|F);
	writer.setDataOptions(N|F);
	writer.write("testGrid_Smooth2D_1_init.vtk");

	Grid_Smooth2D smoother(&m);
	Grid_Smooth2D::STATUS result = smoother.execute();

	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("testGrid_Smooth2D_1_result.vtk");

}










TEST(ClaireTestClass, ResLU_Test1)
{
	Eigen::SparseMatrix<double> A(3,3);
	Eigen::VectorXd b(3);

	A.coeffRef(0,0) = 2 ;
	A.coeffRef(0,1) = -1 ;
	A.coeffRef(0,2) = 0 ;
	A.coeffRef(1,0) = -1 ;
	A.coeffRef(1,1) = 2 ;
	A.coeffRef(1,2) = -1 ;
	A.coeffRef(2,0) = 0 ;
	A.coeffRef(2,1) = -1 ;
	A.coeffRef(2,2) = 2 ;

	b[0] = 1;
	b[1] = 1;
	b[2] = 1;

	ResLU reslu(A, b);
	ResLU::STATUS resultat = reslu.execute();

	Eigen::SparseMatrix<double> LU = reslu.getLU();

	ASSERT_EQ(LU.coeffRef(0,0),2);
	ASSERT_EQ(LU.coeffRef(0,1),-1);
	ASSERT_EQ(LU.coeffRef(0,2),0);
	ASSERT_EQ(LU.coeffRef(1,0),-0.5);
	ASSERT_EQ(LU.coeffRef(1,1),3.0/2.0);
	ASSERT_EQ(LU.coeffRef(1,2),-1);
	ASSERT_EQ(LU.coeffRef(2,0),0);
	ASSERT_EQ(LU.coeffRef(2,1),-2.0/3.0);
	//ASSERT_EQ(LU.coeffRef(2,2),4.0/3.0);

}

