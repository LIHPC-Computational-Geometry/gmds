//
// Created by rochec on 22/03/2022.
//

#include <gmds/math/Line.h>
#include <gmds/claire/Utils.h>
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