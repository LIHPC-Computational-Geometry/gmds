//
// Created by rochec on 12/05/23.
//

#include <gmds/aero/Front.h>
#include <gmds/aero/Utils.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gtest/gtest.h>
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/*                               TESTS UNITAIRES                              */
/*----------------------------------------------------------------------------*/

TEST(FrontTestClass, Front_2D)
{
	// 2D Mesh
	//          n7 o---------------------------o n6
	//             |\                         /|
	//					| \ n3                 n2 / |
	//					|	o---------------------o  |
	//             |  |                     |  |
	//             |  |                     |  |
	//             |  |        EMPTY        |  |
	//             |  |                     |  |
	//             |  |                     |  |
	//             |  |                     |  |
	//					|	o---------------------o  |
	//             | / n0                 n1 \ |
	//             |/                         \|
	//          n4 o---------------------------o n5




	// WE READ
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
	                             gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));

	Node n0 = m.newNode({-0.5, -0.5, 0.0});
	Node n1 = m.newNode({0.5, -0.5, 0.0});
	Node n2 = m.newNode({0.5, 0.5, 0.0});
	Node n3 = m.newNode({-0.5, 0.5, 0.0});

	Node n4 = m.newNode({-1.0, -1.0, 0.0});
	Node n5 = m.newNode({1.0, -1.0, 0.0});
	Node n6 = m.newNode({1.0, 1.0, 0.0});
	Node n7 = m.newNode({-1.0, 1.0, 0.0});

	std::map<TCellID, TCellID> map_nextnodes;
	map_nextnodes[n0.id()] = n4.id();
	map_nextnodes[n1.id()] = n5.id();
	map_nextnodes[n2.id()] = n6.id();
	map_nextnodes[n3.id()] = n7.id();

	Edge e0 = m.newEdge(n0, n1);
	Edge e1 = m.newEdge(n1, n2);
	Edge e2 = m.newEdge(n2, n3);
	Edge e3 = m.newEdge(n3, n0);

	// The front
	Front Fr = Front();

	Variable<int>* couche_id = m.newVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	couche_id->set(n0.id(), 0);
	couche_id->set(n1.id(), 0);
	couche_id->set(n2.id(), 0);
	couche_id->set(n3.id(), 0);

	couche_id->set(n4.id(), -1);
	couche_id->set(n5.id(), -1);
	couche_id->set(n6.id(), -1);
	couche_id->set(n7.id(), -1);

	Fr.initializeFromLayerId(&m, 0);

	ASSERT_EQ(Fr.getNodes().size(), 4);
	ASSERT_EQ(Fr.getEdges().size(), 4);

	Fr.initializeNodeType(&m, map_nextnodes);

	ASSERT_EQ(Fr.isFusionable(n0.id()), true);
	ASSERT_EQ(Fr.isFusionable(n1.id()), true);
	ASSERT_EQ(Fr.isMultiplicable(n1.id()), true);
	ASSERT_EQ(Fr.isMultiplicable(n2.id()), true);
	ASSERT_EQ(Fr.isMultiplicable(n3.id()), true);

	Fr.setMultorFusFromLimits(&m, 0.0, -1.0, -1.0);

	ASSERT_EQ(Fr.isFusionable(n0.id()), false);
	ASSERT_EQ(Fr.isFusionable(n1.id()), true);
	ASSERT_EQ(Fr.isFusionable(n2.id()), true);
	ASSERT_EQ(Fr.isFusionable(n3.id()), false);
	ASSERT_EQ(Fr.isMultiplicable(n0.id()), false);
	ASSERT_EQ(Fr.isMultiplicable(n1.id()), true);
	ASSERT_EQ(Fr.isMultiplicable(n2.id()), true);
	ASSERT_EQ(Fr.isMultiplicable(n3.id()), false);

	Fr.setFrontID(1);
	ASSERT_EQ(Fr.getFrontID(), 1);
	Fr.setFrontID(0);

	Fr.setNonFusionable(n1.id());
	ASSERT_EQ(Fr.isFusionable(n1.id()), false);
	Fr.setNonMultiplicable(n2.id());
	ASSERT_EQ(Fr.isMultiplicable(n2.id()), false);

	ASSERT_EQ(Fr.getNeighbors(n0.id()).size(), 2);
	ASSERT_EQ(Fr.getNeighbors(n3.id()).size(), 2);

	ASSERT_EQ(Fr.getNodeType(n0.id()), 0);
	Fr.setMultipleNode(n0.id());
	ASSERT_EQ(Fr.getNodeType(n0.id()), 1);
	Fr.setContractedNode(n0.id());
	ASSERT_EQ(Fr.getNodeType(n0.id()), 2);

	ASSERT_EQ(Fr.getIdealNode(n1.id()), n5.id());
	ASSERT_EQ(Fr.getIdealNode(n2.id()), n6.id());
	ASSERT_EQ(Fr.getNextNode(n0.id(), n1.id()), n4.id());

	// Creates the faces
	/*
	math::Utils::GetOrCreateQuadAndConnectivities(&m, n0.id(), n1.id(), n2.id(), n3.id());
	math::Utils::GetOrCreateQuadAndConnectivities(&m, n4.id(), n5.id(), n6.id(), n7.id());
	math::Utils::GetOrCreateQuadAndConnectivities(&m, n0.id(), n1.id(), n5.id(), n4.id());
	math::Utils::GetOrCreateQuadAndConnectivities(&m, n0.id(), n3.id(), n7.id(), n4.id());
	math::Utils::GetOrCreateQuadAndConnectivities(&m, n3.id(), n2.id(), n6.id(), n7.id());
	math::Utils::GetOrCreateQuadAndConnectivities(&m, n2.id(), n1.id(), n5.id(), n6.id());


	Variable<int>* var_LayerID = m.newVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	std::vector<TCellID> nodes_id;
	std::vector<TCellID> faces_id;
	for (auto n_id : m.nodes()) {
		var_LayerID->set(n_id, 0);
		nodes_id.push_back(n_id);
	}
	for (auto f_id:m.faces())
	{
		faces_id.push_back(f_id);
	}

	ASSERT_EQ(nodes_id.size(), 8);
	ASSERT_EQ(faces_id.size(), 6);

	Front_3D Front = Front_3D(1, nodes_id, faces_id);

	ASSERT_EQ(Front.getFrontID(), 1);
	Front.setFrontID(0);
	ASSERT_EQ(Front.getFrontID(), 0);

	std::vector<TCellID> tmp_nodes = Front.getNodes();
	std::vector<TCellID> tmp_faces = Front.getFaces();

	ASSERT_EQ(tmp_nodes.size(), 8);
	ASSERT_EQ(tmp_faces.size(), 6);

	std::vector<TCellID> ord_edges_around_n0 = Front.orderedFrontEdgesAroundNode(&m, n0.id());
	ASSERT_EQ(ord_edges_around_n0.size(), 3);

	TCellID e_id = math::Utils::CommonEdge(&m, n0.id(), n1.id());
	std::vector<TCellID> adjfaces = Front.edgeFacesOnFront(&m, e_id);
	ASSERT_EQ(adjfaces.size(), 2);

	 */
}