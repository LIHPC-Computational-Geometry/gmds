//
// Created by rochec on 22/03/2022.
//

#include <gmds/math/Line.h>
#include <gmds/claire/Utils.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/math/BezierSurface.h>
#include <gmds/math/BezierHex.h>
#include <gmds/claire/Blocking3D.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gtest/gtest.h>
#include <iostream>
#include <unit_test_config.h>
#include <Eigen/Sparse>
#include <Eigen/Eigen>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

TEST(UtilsTestClass, Utils_Test1)
{
	// Test de la méthode math::Utils::distFromNodeIds

	Mesh m(MeshModel(DIM3 | R | F | N | F2N | N2F));

	GridBuilder gb(&m,2);
	gb.execute(3,1.0, 4, 1.0);

	double eps(pow(10,-5));

	ASSERT_NEAR(math::Utils::distFromNodeIds(&m, 0, 8), 2.0, eps);
	ASSERT_NEAR(math::Utils::distFromNodeIds(&m, 0, 3), 3.0, eps);
	ASSERT_NEAR(math::Utils::distFromNodeIds(&m, 0, 11), 3.60555, eps);
	ASSERT_NEAR(math::Utils::distFromNodeIds(&m, 3, 1), 2.0, eps);
	ASSERT_NEAR(math::Utils::distFromNodeIds(&m, 7, 9), 2.23607, eps);
	ASSERT_NEAR(math::Utils::distFromNodeIds(&m, 1, 7), 2.23607, eps);
	ASSERT_NEAR(math::Utils::distFromNodeIds(&m, 10, 3), 2.23607, eps);

	//std::cout << math::Utils::distFromNodeIds(&m, 10, 3) << std::endl;

	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("Utils_Test1.vtk");

}
/*----------------------------------------------------------------------------*/
TEST(UtilsTestClass, Utils_Test2)
{
	// Test de la méthode math::Utils::CommonEdge

	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::N2E | gmds::N2F | gmds::F2N | gmds::E2N | gmds::F2E | gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir + "/Aero/Poubelle/Carre.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	ASSERT_EQ(math::Utils::CommonEdge(&m, 0, 4), 212);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 51, 54), 141);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 62, 85), 14);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 16, 17), 115);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 78, 67), 1);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 17, 40), 116);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 86, 51), 203);

	ASSERT_EQ(math::Utils::CommonEdge(&m, 44, 84), NullID);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 56, 97), NullID);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 45, 38), NullID);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 87, 11), NullID);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 20, 16), NullID);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 24, 25), NullID);
	ASSERT_EQ(math::Utils::CommonEdge(&m, 14, 12), NullID);

	/*
	TCellID id1 = 86;
	TCellID id2 = 51;
	std::cout << math::Utils::CommonEdge(&m, id1, id2) << std::endl;
	Node n1 = m.get<Node>(id1);
	Node n2 = m.get<Node>(id2);
	std::vector<Edge> adj_edges_1 = n1.get<Edge>() ;
	std::vector<Edge> adj_edges_2 = n2.get<Edge>();
	std::cout << "---------" << std::endl;
	for (auto e:adj_edges_1){
		std::cout << e.id() << std::endl;
	}
	std::cout << "---------" << std::endl;
	for (auto e:adj_edges_2){
		std::cout << e.id() << std::endl;
	}
	 */

	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("Utils_Test2.vtk");

}
/*----------------------------------------------------------------------------*/
TEST(UtilsTestClass, Utils_Test3)
{
	// Test de la méthode math::Utils::MeshCleaner

	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::N2E | gmds::N2F | gmds::F2N | gmds::E2N | gmds::F2E | gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir + "/Aero/Poubelle/Carre.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	math::Point p0(-1,-1,0);
	Node n0 = m.newNode(p0);

	math::Point p1(-60,-12,0);
	Node n1 = m.newNode(p1);

	math::Point p2(3,10,0);
	Node n2 = m.newNode(p2);

	math::Point p3(-13,4,0);
	Node n3 = m.newNode(p3);

	math::Utils::MeshCleaner(&m);

	for (auto n_id:m.nodes()){
		ASSERT_NE(n_id, n0.id());
		ASSERT_NE(n_id, n1.id());
		ASSERT_NE(n_id, n2.id());
		ASSERT_NE(n_id, n3.id());
	}

	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("Utils_Test3.vtk");

}
/*----------------------------------------------------------------------------*/
TEST(ClaireTestClass, Utils_AdjacentNodes)
{
	// Test
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::N2E | gmds::N2F | gmds::F2N | gmds::E2N | gmds::F2E | gmds::E2F));

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir + "/Aero/Poubelle/Carre.vtk";

	gmds::IGMeshIOService ioService(&m);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N | gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	TCellID n_id = 1;
	Node n = m.get<Node>(n_id);
	std::vector<Node> adj_nodes = math::Utils::AdjacentNodes(&m, n);
	ASSERT_EQ(adj_nodes.size(),3);
	ASSERT_EQ(adj_nodes[0].id(),11);
	ASSERT_EQ(adj_nodes[1].id(),91);
	ASSERT_EQ(adj_nodes[2].id(),10);

	n_id = 4;
	n = m.get<Node>(n_id);
	adj_nodes = math::Utils::AdjacentNodes(&m, n);
	ASSERT_EQ(adj_nodes.size(),4);
	ASSERT_EQ(adj_nodes[0].id(),5);
	ASSERT_EQ(adj_nodes[1].id(),79);
	ASSERT_EQ(adj_nodes[2].id(),0);
	ASSERT_EQ(adj_nodes[3].id(),90);

	n_id = 92;
	n = m.get<Node>(n_id);
	adj_nodes = math::Utils::AdjacentNodes(&m, n);
	ASSERT_EQ(adj_nodes.size(), 6);
	ASSERT_EQ(adj_nodes[0].id(), 69);
	ASSERT_EQ(adj_nodes[1].id(), 70);
	ASSERT_EQ(adj_nodes[2].id(), 71);
	ASSERT_EQ(adj_nodes[3].id(), 50);
	ASSERT_EQ(adj_nodes[4].id(), 86);
	ASSERT_EQ(adj_nodes[5].id(), 84);
	//std::cout << adj_nodes[5].id() << std::endl;

	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("Utils_Test4.vtk");

}
/*----------------------------------------------------------------------------*/
TEST(ClaireTestClass, Utils_BuildMesh2DFromBlocking2D)
{
	// Test
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::N2E | gmds::N2F | gmds::F2N | gmds::E2N | gmds::F2E | gmds::E2F));

	Blocking2D b;
	Node n1 = b.newBlockCorner(0,0);
	Node n2 = b.newBlockCorner(1,0);
	Node n3 = b.newBlockCorner(1,1);
	Node n4=  b.newBlockCorner(0,1);

	Blocking2D::Block b0 = b.newBlock(n1,n2,n3,n4);

	Node n5 = b.newBlockCorner(2,0,0);
	Node n6 = b.newBlockCorner(2,1.5,0);
	Blocking2D::Block b1 = b.newBlock(n2,n5,n6,n3);

	ASSERT_EQ(b0.id(), b.block(0).id());
	b0.setNbDiscretizationI(11);
	b0.setNbDiscretizationJ(11);
	b1.setNbDiscretizationI(11);
	b1.setNbDiscretizationJ(11);
	b.initializeGridPoints();

	TInt mark_block_nodes = m.newMark<Node>();
	TInt mark_first_layer = m.newMark<Node>();
	TInt mark_farfield_nodes = m.newMark<Node>();
	math::Utils::BuildMesh2DFromBlocking2D(&b, &m, mark_block_nodes, mark_first_layer, mark_farfield_nodes);
	m.unmarkAll<Node>(mark_block_nodes);
	m.freeMark<Node>(mark_block_nodes);
	m.unmarkAll<Node>(mark_first_layer);
	m.freeMark<Node>(mark_first_layer);
	m.unmarkAll<Node>(mark_farfield_nodes);
	m.freeMark<Node>(mark_farfield_nodes);


	ASSERT_EQ(m.getNbFaces(), 200);
	ASSERT_EQ(m.getNbNodes(), 231);

	IGMeshIOService ioService_geom(&b);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("Utils_BuildMesh2DFromBlocking2D_Blocking.vtk");

	IGMeshIOService ioService_geom_mesh(&m);
	VTKWriter writer_geom_mesh(&ioService_geom_mesh);
	writer_geom_mesh.setCellOptions(N|F);
	writer_geom_mesh.setDataOptions(N|F);
	writer_geom_mesh.write("Utils_BuildMesh2DFromBlocking2D_Mesh.vtk");

}
/*----------------------------------------------------------------------------*/
TEST(ClaireTestClass, Utils_WeightedPointOnBranch)
{
	{
		math::Point A({0.0, 0.0, 0.0});
		math::Point B({1.0, 0.0, 0.0});
		math::Point C({3.0, 0.0, 0.0});
		math::Point D = math::Utils::WeightedPointOnBranch(A, B, C, 0.5);
		ASSERT_FLOAT_EQ(D.X(), 1.5);
		ASSERT_FLOAT_EQ(D.Y(), 0.0);
	}

	{
		math::Point A({0.0, 0.0, 0.0});
		math::Point B({1.0, 0.0, 0.0});
		math::Point C({4.0, 0.0, 0.0});
		math::Point D = math::Utils::WeightedPointOnBranch(A, B, C, 0.75);
		ASSERT_FLOAT_EQ(D.X(), 3.0);
		ASSERT_FLOAT_EQ(D.Y(), 0.0);
	}

	{
		math::Point A({0.0, 0.0, 0.0});
		math::Point B({1.0, 0.0, 0.0});
		math::Point C({1.0, 3.0, 0.0});
		math::Point D = math::Utils::WeightedPointOnBranch(A, B, C, 0.75);
		ASSERT_FLOAT_EQ(D.X(), 1.0);
		ASSERT_FLOAT_EQ(D.Y(), 2.0);
	}

}
/*----------------------------------------------------------------------------*/
TEST(ClaireTestClass, Utils_isInTriangle)
{
	{
		math::Point T1({0.0, 0.0, 0.0});
		math::Point T2({1.0, 0.0, 0.0});
		math::Point T3({0.0, 1.0, 0.0});
		ASSERT_EQ(math::Utils::isInTriangle(T1, T2, T3, {0.25, 0.25, 0.0}), true);
		ASSERT_EQ(math::Utils::isInTriangle(T1, T2, T3, {1.0, 1.0, 0.0}), false);
		ASSERT_EQ(math::Utils::isInTriangle(T1, T2, T3, {0.5, 0.0, 0.0}), true);
	}

}
/*----------------------------------------------------------------------------*/
TEST(ClaireTestClass, Utils_isInTetra)
{
	{
		math::Point T1({0.0, 0.0, 0.0});
		math::Point T2({1.0, 0.0, 0.0});
		math::Point T3({0.0, 1.0, 0.0});
		math::Point T4({0.0, 0.0, 1.0});

		ASSERT_EQ(math::Utils::isInTetra(T1, T2, T3, T4, {0.25, 0.25, 0.0}), true);
		ASSERT_EQ(math::Utils::isInTetra(T1, T2, T3, T4, {0, 0.25, 0.0}), true);
		ASSERT_EQ(math::Utils::isInTetra(T1, T2, T3, T4, {1.0, 1.0, 1.0}), false);
	}

	{
		math::Point T1({-1.76213, -1.12673, -0.129545});
		math::Point T2({-1.28791, -1.47346, -0.211251});
		math::Point T3({-1, -1, -0.5});
		math::Point T4({-1.45698, -0.699361, -0.0343749});

		math::Point M({-1, -1, -0.333333});

		ASSERT_EQ(math::Utils::isInTetra(T1, T2, T3, T4, M), false);
	}

}
/*----------------------------------------------------------------------------*/
TEST(ClaireTestClass, Utils_CreateQuadAndConnectivities)
{
	// Test
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::N2E | gmds::N2F | gmds::F2N | gmds::E2N | gmds::F2E | gmds::E2F));

	Node n0 = m.newNode({0,0,0});
	Node n1 = m.newNode({0,1,0});
	Node n2 = m.newNode({1,1,0});
	Node n3 = m.newNode({1,0,0});

	Node n4 = m.newNode({2,0,0});
	Node n5 = m.newNode({2,1,0});

	ASSERT_EQ(m.getNbFaces(), 0);

	TCellID f_id = math::Utils::GetOrCreateQuadAndConnectivities(&m, n0.id(), n1.id(), n2.id(), n3.id());
	ASSERT_EQ(m.getNbFaces(), 1);
	ASSERT_EQ((n0.get<Face>()).size(), 1);
	ASSERT_EQ((n0.get<Edge>()).size(), 2);
	ASSERT_EQ((n3.get<Face>()).size(), 1);
	ASSERT_EQ((n3.get<Edge>()).size(), 2);

	TCellID f2_id = math::Utils::GetOrCreateQuadAndConnectivities(&m,n2.id(), n3.id(), n4.id(), n5.id());
	ASSERT_EQ(m.getNbFaces(), 2);
	ASSERT_EQ((n0.get<Face>()).size(), 1);
	ASSERT_EQ((n0.get<Edge>()).size(), 2);
	ASSERT_EQ((n3.get<Face>()).size(), 2);
	ASSERT_EQ((n3.get<Edge>()).size(), 3);

}
/*----------------------------------------------------------------------------*/
TEST(ClaireTestClass, Utils_CreateHexaNConnectivities)
{
	// Test
	gmds::Mesh m(gmds::MeshModel(gmds::MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
																F2E | E2F | R2E | E2R | N2R | N2F | N2E )));

	Node n0 = m.newNode({0,0,0});
	Node n1 = m.newNode({1,0,0});
	Node n2 = m.newNode({1,1,0});
	Node n3 = m.newNode({0,1,0});

	Node n4 = m.newNode({0,0,1});
	Node n5 = m.newNode({1,0,1});
	Node n6 = m.newNode({1,1,1});
	Node n7 = m.newNode({0,1,1});

	Node n8 = m.newNode({2,0,0});
	Node n9 = m.newNode({2,1,0});

	Node n10 = m.newNode({2,0,1});
	Node n11 = m.newNode({2,1,1});

	ASSERT_EQ(m.getNbFaces(), 0);
	ASSERT_EQ(m.getNbNodes(), 12);

	// Create a first hexa will all the connectivities
	TCellID r1_id = math::Utils::CreateHexaNConnectivities(&m, n0, n1, n2, n3, n4, n5, n6, n7);

	ASSERT_EQ(m.getNbEdges(), 12);
	ASSERT_EQ(m.getNbFaces(), 6);
	ASSERT_EQ(m.getNbRegions(), 1);
	ASSERT_EQ((n0.get<Edge>()).size(), 3);
	ASSERT_EQ((n0.get<Face>()).size(), 3);
	ASSERT_EQ((n0.get<Region>()).size(), 1);
	ASSERT_EQ((n1.get<Edge>()).size(), 3);
	ASSERT_EQ((n1.get<Face>()).size(), 3);
	ASSERT_EQ((n1.get<Region>()).size(), 1);
	ASSERT_EQ((n8.get<Edge>()).size(), 0);
	ASSERT_EQ((n8.get<Face>()).size(), 0);
	ASSERT_EQ((n8.get<Region>()).size(), 0);

	// Create a second hexa will all the connectivities
	TCellID r2_id = math::Utils::CreateHexaNConnectivities(&m, n1, n8, n9, n2, n5, n10, n11, n6);

	ASSERT_EQ(m.getNbEdges(), 20);
	ASSERT_EQ(m.getNbFaces(), 11);
	ASSERT_EQ(m.getNbRegions(), 2);
	ASSERT_EQ((n0.get<Edge>()).size(), 3);
	ASSERT_EQ((n0.get<Face>()).size(), 3);
	ASSERT_EQ((n0.get<Region>()).size(), 1);
	ASSERT_EQ((n1.get<Edge>()).size(), 4);
	ASSERT_EQ((n1.get<Face>()).size(), 5);
	ASSERT_EQ((n1.get<Region>()).size(), 2);
	ASSERT_EQ((n8.get<Edge>()).size(), 3);
	ASSERT_EQ((n8.get<Face>()).size(), 3);
	ASSERT_EQ((n8.get<Region>()).size(), 1);

}
/*----------------------------------------------------------------------------*/
TEST(ClaireTestClass, Utils_minEdgeLenght)
{
	// Test
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::N2E | gmds::N2F | gmds::F2N | gmds::E2N | gmds::F2E | gmds::E2F));

	Node n0 = m.newNode({0,0,0});
	Node n1 = m.newNode({1,0,0});
	Node n2 = m.newNode({1,1,0});
	Node n3 = m.newNode({0,1,0});

	m.newTriangle(n0,n1,n3);
	m.newTriangle(n3,n1,n2);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	ASSERT_FLOAT_EQ(1.0, math::Utils::minEdgeLenght(&m));

}
/*----------------------------------------------------------------------------*/
TEST(ClaireTestClass, Utils_getFacesAdjToEdgeInHexa)
{
	// Test
	gmds::Mesh m(gmds::MeshModel(gmds::MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                                             F2E | E2F | R2E | E2R | N2R | N2F | N2E )));

	Node n0 = m.newNode({0,0,0});
	Node n1 = m.newNode({1,0,0});
	Node n2 = m.newNode({1,1,0});
	Node n3 = m.newNode({0,1,0});

	Node n4 = m.newNode({0,0,1});
	Node n5 = m.newNode({1,0,1});
	Node n6 = m.newNode({1,1,1});
	Node n7 = m.newNode({0,1,1});

	Node n8 = m.newNode({2,0,0});
	Node n9 = m.newNode({2,1,0});

	Node n10 = m.newNode({2,0,1});
	Node n11 = m.newNode({2,1,1});

	ASSERT_EQ(m.getNbFaces(), 0);
	ASSERT_EQ(m.getNbNodes(), 12);

	// Create a first hexa will all the connectivities
	TCellID r1_id = math::Utils::CreateHexaNConnectivities(&m, n0, n1, n2, n3, n4, n5, n6, n7);

	// Create a second hexa will all the connectivities
	TCellID r2_id = math::Utils::CreateHexaNConnectivities(&m, n1, n8, n9, n2, n5, n10, n11, n6);

	ASSERT_EQ(m.getNbRegions(), 2);

	TCellID e_id = math::Utils::CommonEdge(&m, n0.id(), n1.id());
	std::vector<Face> adj_faces = math::Utils::getFacesAdjToEdgeInHexa(&m, e_id, r1_id);

	ASSERT_EQ(adj_faces.size(), 2);

}
/*----------------------------------------------------------------------------*/
TEST(ClaireTestClass, Utils_cutMeshUnderXAxis)
{
	gmds::Mesh m(gmds::MeshModel(gmds::MeshModel(DIM2 | F | E | N | F2N | E2N |
	                                             F2E | E2F | N2F | N2E )));

	Node n0 = m.newNode({0,0,0});
	Node n1 = m.newNode({1,0,0});
	Node n2 = m.newNode({1,1,0});
	Node n3 = m.newNode({0,1,0});

	Node n4 = m.newNode({1,-1,0});
	Node n5 = m.newNode({0,-1,0});

	Node n6 = m.newNode({2,0,0});
	Node n7 = m.newNode({2,1,0});

	m.newQuad(n0,n1,n2,n3);
	m.newQuad(n0,n1,n4,n5);
	m.newQuad(n1,n6,n7,n2);

	ASSERT_EQ(m.getNbNodes(), 8);
	ASSERT_EQ(m.getNbFaces(), 3);

	math::Utils::cutAxiBlocking2D(&m);

	ASSERT_EQ(m.getNbNodes(), 6);
	ASSERT_EQ(m.getNbFaces(), 2);
}
/*----------------------------------------------------------------------------*/
TEST(ClaireTestClass, Utils_oppositeEdgeInFace)
{
	gmds::Mesh m(gmds::MeshModel(gmds::MeshModel(DIM2 | F | E | N | F2N | E2N |
	                                             F2E | E2F | N2F | N2E )));

	Node n0 = m.newNode({0,0,0});
	Node n1 = m.newNode({1,0,0});
	Node n2 = m.newNode({1,1,0});
	Node n3 = m.newNode({0,1,0});

	m.newQuad(n0,n1,n2,n3);

	gmds::MeshDoctor doc(&m);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	Edge e_opp = math::Utils::oppositeEdgeInFace(&m,0,0);

	ASSERT_EQ(e_opp.id(), 2);
}
/*----------------------------------------------------------------------------*/
TEST(ClaireTestClass, Utils_buildEfromFandConnectivies)
{
	gmds::Mesh m(gmds::MeshModel(gmds::MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                                             F2E | E2F | R2E | E2R | N2R | N2F | N2E)));

	Node n0 = m.newNode({0,0,0});
	Node n1 = m.newNode({1,0,0});
	Node n2 = m.newNode({1,1,0});
	Node n3 = m.newNode({0,1,0});

	Node n4 = m.newNode({2,0,0});
	Node n5 = m.newNode({2,1,0});

	Face f0 = m.newQuad(n0,n1,n2,n3);
	Face f1 = m.newQuad(n1,n4,n5,n2);

	math::Utils::buildEfromFandConnectivies(&m);

	ASSERT_EQ(m.getNbEdges(), 7);
	ASSERT_EQ(n0.get<Face>().size(), 1);
	ASSERT_EQ(n1.get<Face>().size(), 2);
	ASSERT_EQ(n1.get<Edge>().size(), 3);

	TCellID e_id = math::Utils::CommonEdge(&m, n1.id(), n2.id());
	Edge e = m.get<Edge>(e_id);
	ASSERT_EQ(e.get<Face>().size(), 2);
	ASSERT_EQ(e.get<Node>().size(), 2);

	ASSERT_EQ(f0.get<Node>().size(), 4);
	ASSERT_EQ(f1.get<Edge>().size(), 4);

}
/*----------------------------------------------------------------------------*/
TEST(ClaireTestClass, Test_BlockingControlPoints)
{
	gmds::Mesh m = Mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                        F2E | E2F | R2E | E2R | N2R | N2F | N2E));

	Node n1 = m.newNode(0,0);
	Node n2 = m.newNode(1,0);
	Node n3 = m.newNode(1,1);
	Node n4 = m.newNode(0,1);

	Node n5 = m.newNode(2,0);
	Node n6 = m.newNode(2,1);

	m.newQuad(n1, n2, n3, n4);
	m.newQuad(n2, n5, n6, n3);

	Blocking2D b(m);
	Blocking2D b_controlpoints(m);

	int degree_Bezier(5);
	int bloc_discretization(20);

	for (auto bloc:b_controlpoints.allBlocks())
	{
		bloc.setNbDiscretizationI(degree_Bezier+1);
		bloc.setNbDiscretizationJ(degree_Bezier+1);
	}

	for (auto bloc:b.allBlocks())
	{
		bloc.setNbDiscretizationI(bloc_discretization);
		bloc.setNbDiscretizationJ(bloc_discretization);
	}

	b.initializeGridPoints();
	b_controlpoints.initializeGridPoints();
	srand(time(NULL));

	Variable<int>* var_couche_b = b.newVariable<int, GMDS_NODE>("GMDS_Couche");
	for (auto bloc:b.allBlocks())
	{
		for (int i=0;i<bloc.getNbDiscretizationI();i++)
		{
			for (int j=1;j<bloc.getNbDiscretizationJ();j++)
			{
				var_couche_b->set(bloc(i,j).id(), 1);
			}
		}
	}

	/*
	for (auto bloc:b_controlpoints.allBlocks())
	{
		// Control Points on EDGES
		for (int j=1;j<=degree_Bezier-1;j++)
		{
			double random_x = rand()/(5.0*RAND_MAX) -0.1 ;
			double random_y = rand()/(5.0*RAND_MAX) -0.1 ;
			bloc(0,j).setX(bloc(0,j).X()+random_x);
			bloc(0,j).setY(bloc(0,j).Y()+random_y);

			random_x = rand()/(5.0*RAND_MAX) -0.1 ;
			random_y = rand()/(5.0*RAND_MAX) -0.1 ;
			bloc(degree_Bezier,j).setX(bloc(degree_Bezier,j).X()+random_x);
			bloc(degree_Bezier,j).setY(bloc(degree_Bezier,j).Y()+random_y);
		}
		for (int i=1;i<=degree_Bezier-1;i++)
		{
			double random_x = rand()/(5.0*RAND_MAX) -0.1 ;
			double random_y = rand()/(5.0*RAND_MAX) -0.1 ;
			bloc(i,0).setX(bloc(i,0).X()+random_x);
			bloc(i,0).setY(bloc(i,0).Y()+random_y);

			random_x = rand()/(5.0*RAND_MAX) -0.1 ;
			random_y = rand()/(5.0*RAND_MAX) -0.1 ;
			bloc(i,degree_Bezier).setX(bloc(i,degree_Bezier).X()+random_x);
			bloc(i,degree_Bezier).setY(bloc(i,degree_Bezier).Y()+random_y);
		}
		// Control Points on SURFACE
		for (int i=1;i<=degree_Bezier-1;i++)
		{
		   for (int j=1;j<=degree_Bezier-1;j++)
		   {
		      double random_x = rand()/(5.0*RAND_MAX) -0.1 ;
		      double random_y = rand()/(5.0*RAND_MAX) -0.1 ;
		      bloc(i,j).setX(bloc(i,j).X()+random_x);
		      bloc(i,j).setY(bloc(i,j).Y()+random_y);
		   }
		}
	}
	*/

	// Try to compute the control points to interpolate
	//Eigen::Matrix4d mat_B;
	//Eigen::Vector4d ctrl_points_x;
	//Eigen::Vector4d ctrl_points_y;
	//Eigen::VectorXd interp_points_x;
	//Eigen::Vector4d interp_points_y;

	Eigen::MatrixXd mat_B(degree_Bezier+1, degree_Bezier+1);
	Eigen::VectorXd ctrl_points_x(degree_Bezier+1);
	Eigen::VectorXd ctrl_points_y(degree_Bezier+1);
	Eigen::VectorXd interp_points_x(degree_Bezier+1);
	Eigen::VectorXd interp_points_y(degree_Bezier+1);

	interp_points_x[0] = 0.0 ;
	interp_points_y[0] = 0.0 ;

	for (int i=1;i<degree_Bezier;i++)
	{
		interp_points_x[i] = 1.0*i/degree_Bezier;
		interp_points_y[i] = 0.05;
	}

	/*
	interp_points_x[1] = 0.25 ;
	interp_points_y[1] = 0.05 ;

	interp_points_x[2] = 0.5 ;
	interp_points_y[2] = 0.05 ;

	interp_points_x[3] = 0.75 ;
	interp_points_y[3] = 0.05 ;
	 */

	interp_points_x[degree_Bezier] = 1.0 ;
	interp_points_y[degree_Bezier] = 0.0 ;

	// Matrix Assembly
	for (int i=0;i<=degree_Bezier;i++)
	{
		for (int j=0;j<=degree_Bezier;j++)
		{
			//std::cout << "Degree " << degree_Bezier << ", " << j << ", u: " << 1.0*i/degree_Bezier << std::endl;
			double bij = math::Utils::BernsteinPolynomial(degree_Bezier, j, 1.0*i/degree_Bezier);
			mat_B(i,j) = bij;
		}
	}

	std::cout << mat_B << std::endl;

	//Eigen::Matrix4d mat_B_inv = mat_B.inverse();
	Eigen::MatrixXd mat_B_inv = mat_B.inverse();

	std::cout << mat_B_inv << std::endl;

	ctrl_points_x = mat_B_inv*interp_points_x;
	ctrl_points_y = mat_B_inv*interp_points_y;

	Blocking2D::Block bloc = b_controlpoints.block(0) ;
	for (int i=0;i<=degree_Bezier;i++)
	{
		bloc(i,0).setX(ctrl_points_x[i]);
		bloc(i,0).setY(ctrl_points_y[i]);
	}




	// Set the bloc discretization
	// Rely on the fact the blocks have the same IDs between b and b_controlpoints
	for (auto bloc:b.allBlocks())
	{
		Array2D<math::Point> Ctrl_Pts(degree_Bezier+1,degree_Bezier+1);
		Blocking2D::Block b_ctrl_pts = b_controlpoints.block(bloc.id()) ;
		for (int i=0;i<=degree_Bezier;i++)
		{
			for (int j=0;j<=degree_Bezier;j++)
			{
				Ctrl_Pts(i,j) = b_ctrl_pts(i,j).point();
			}
		}
		math::BezierSurface curved_bezier_surface(Ctrl_Pts);
		for (int i=0;i<bloc.getNbDiscretizationI();i++)
		{
			for (int j=0;j<bloc.getNbDiscretizationJ();j++)
			{
				double u = 1.0*i/(bloc.getNbDiscretizationI()-1) ;
				double v = 1.0*j/(bloc.getNbDiscretizationJ()-1.0) ;
				bloc(i,j).setPoint(curved_bezier_surface(u,v));
			}
		}
	}




	// Write the mesh results
	IGMeshIOService ioService(&b);
	VTKWriter writer(&ioService);
	writer.setCellOptions(N|F);
	writer.setDataOptions(N|F);
	writer.write("Utils_Test_BlockingControlPoints.vtk");

	IGMeshIOService ioService_ctrl_points(&b_controlpoints);
	VTKWriter writer_ctrl_points(&ioService_ctrl_points);
	writer_ctrl_points.setCellOptions(N|F);
	writer_ctrl_points.setDataOptions(N|F);
	writer_ctrl_points.write("Utils_Test_BlockingControlPoints_CtrlPoints.vtk");

	gmds::Mesh m_Blocking(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::E | gmds::N2E | gmds::N2F | gmds::F2N | gmds::E2N | gmds::F2E | gmds::E2F));
	TInt mark_block_nodes = m_Blocking.newMark<Node>();
	TInt mark_first_layer = m_Blocking.newMark<Node>();
	TInt mark_farfield_nodes = m_Blocking.newMark<Node>();
	math::Utils::BuildMesh2DFromBlocking2D(&b, &m_Blocking, mark_block_nodes, mark_first_layer, mark_farfield_nodes);
	m_Blocking.unmarkAll<Node>(mark_block_nodes);
	m_Blocking.freeMark<Node>(mark_block_nodes);
	m_Blocking.unmarkAll<Node>(mark_first_layer);
	m_Blocking.freeMark<Node>(mark_first_layer);
	m_Blocking.unmarkAll<Node>(mark_farfield_nodes);
	m_Blocking.freeMark<Node>(mark_farfield_nodes);

	IGMeshIOService ioService_m_blocking(&m_Blocking);
	VTKWriter writer_m_blocking(&ioService_m_blocking);
	writer_m_blocking.setCellOptions(N|F);
	writer_m_blocking.setDataOptions(N|F);
	writer_m_blocking.write("Utils_Test_BlockingControlPoints_Mesh.vtk");

	ASSERT_EQ(0, 0);

}
/*----------------------------------------------------------------------------*/
TEST(ClaireTestClass, Test_resizeMesh)
{
	gmds::Mesh m = Mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                              F2E | E2F | R2E | E2R | N2R | N2F | N2E));

	Node n1 = m.newNode(0,0);
	Node n2 = m.newNode(1,0);
	Node n3 = m.newNode(1,1);
	Node n4 = m.newNode(0,1);

	Face f = m.newQuad(n1, n2, n3, n4);

	math::Utils::resizeMesh(&m, 10.0);

	ASSERT_FLOAT_EQ(f.get<Node>()[1].X(), 10.0);
	ASSERT_FLOAT_EQ(f.get<Node>()[2].X(), 10.0);
	ASSERT_FLOAT_EQ(f.get<Node>()[2].Y(), 10.0);

}
/*----------------------------------------------------------------------------*/