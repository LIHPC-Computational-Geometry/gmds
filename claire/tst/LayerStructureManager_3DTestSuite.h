//
// Created by rochec on 11/10/2023.
//
#include <gmds/claire/LayerStructureManager_3D.h>
#include <gmds/claire/Utils.h>
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
TEST(ClaireTestClass, test_LayerStructureManager_3D)
{
	gmds::Mesh m(gmds::MeshModel(gmds::MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R | F2E | E2F | R2E | E2R | N2R | N2F | N2E)));

	// Create nodes
	Node n0 = m.newNode({0, 0, 0});
	Node n1 = m.newNode({0, 0, 1});
	Node n2 = m.newNode({1, -1, 0});
	Node n3 = m.newNode({1, -1, 1});
	Node n4 = m.newNode({1, -0.2, 0});
	Node n5 = m.newNode({1, -0.2, 1});
	Node n6 = m.newNode({1, 0.2, 0});
	Node n7 = m.newNode({1, 0.2, 1});
	Node n8 = m.newNode({1, 1, 0});
	Node n9 = m.newNode({1, 1, 1});
	Node n10 = m.newNode({2, -1, 0});
	Node n11 = m.newNode({2, -1, 1});
	Node n12 = m.newNode({2, 0, 0});
	Node n13 = m.newNode({2, 0, 1});
	Node n14 = m.newNode({2, 1, 0});
	Node n15 = m.newNode({2, 1, 1});

	std::vector<TCellID> F0_nodes;
	for (auto n_id:m.nodes())
	{
		F0_nodes.push_back(n_id);
	}

	// Create faces
	TCellID f0_id = math::Utils::GetOrCreateQuadAndConnectivities(&m, n0.id(), n1.id(), n5.id(), n4.id());
	TCellID f1_id = math::Utils::GetOrCreateQuadAndConnectivities(&m, n0.id(), n1.id(), n7.id(), n6.id());
	TCellID f2_id = math::Utils::GetOrCreateQuadAndConnectivities(&m, n1.id(), n5.id(), n13.id(), n7.id());
	TCellID f3_id = math::Utils::GetOrCreateQuadAndConnectivities(&m, n0.id(), n4.id(), n12.id(), n6.id());
	TCellID f4_id = math::Utils::GetOrCreateQuadAndConnectivities(&m, n2.id(), n3.id(), n5.id(), n4.id());
	TCellID f5_id = math::Utils::GetOrCreateQuadAndConnectivities(&m, n6.id(), n7.id(), n9.id(), n8.id());
	TCellID f6_id = math::Utils::GetOrCreateQuadAndConnectivities(&m, n10.id(), n11.id(), n13.id(), n12.id());
	TCellID f7_id = math::Utils::GetOrCreateQuadAndConnectivities(&m, n12.id(), n13.id(), n15.id(), n14.id());
	TCellID f8_id = math::Utils::GetOrCreateQuadAndConnectivities(&m, n6.id(), n12.id(), n14.id(), n8.id());
	TCellID f9_id = math::Utils::GetOrCreateQuadAndConnectivities(&m, n4.id(), n12.id(), n10.id(), n2.id());
	TCellID f10_id = math::Utils::GetOrCreateQuadAndConnectivities(&m, n7.id(), n9.id(), n15.id(), n13.id());
	TCellID f11_id = math::Utils::GetOrCreateQuadAndConnectivities(&m, n11.id(), n3.id(), n5.id(), n13.id());
	TCellID f12_id = math::Utils::GetOrCreateQuadAndConnectivities(&m, n8.id(), n9.id(), n15.id(), n14.id());
	TCellID f13_id = math::Utils::GetOrCreateQuadAndConnectivities(&m, n2.id(), n10.id(), n11.id(), n3.id());

	std::vector<TCellID> F0_faces;
	for (auto f_id:m.faces())
	{
		F0_faces.push_back(f_id);
	}

	// Create the Front
	Front_3D F0 = Front_3D(0, F0_nodes, F0_faces);

	Variable<int>* var_front_edges_classification = m.getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
	Variable<int>* var_node_couche_id = m.getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");

	for (auto n_id:F0_nodes)
	{
		var_node_couche_id->set(n_id, 0);
	}

	// Write the front F0
	IGMeshIOService ioService_geom(&m);
	VTKWriter writer_geom(&ioService_geom);
	writer_geom.setCellOptions(N|F);
	writer_geom.setDataOptions(N|F);
	writer_geom.write("TestSuite_LayerStructureManager_3D.vtk");

	// Create fake nodes on each node of the front F0
	std::map<TCellID, TCellID> map_new_nodes;
	for (auto n_id:F0_nodes)
	{
		Node n_new = m.newNode(m.get<Node>(n_id).point());
		map_new_nodes[n_id] = n_new.id() ;
		var_node_couche_id->set(n_new.id(), 1);
	}

	// Classify some edges
	TCellID e_corner_id = math::Utils::CommonEdge(&m, n8.id(), n14.id());
	TCellID e_end_id = math::Utils::CommonEdge(&m, n4.id(), n5.id());
	TCellID e_reversal_id = math::Utils::CommonEdge(&m, n0.id(), n1.id());

	var_front_edges_classification->set(e_corner_id, 1);
	var_front_edges_classification->set(e_end_id, 2);
	var_front_edges_classification->set(e_reversal_id, 3);

	// Create Structure Manager
	LayerStructureManager_3D lsm = LayerStructureManager_3D(&m, &F0, map_new_nodes);

	ASSERT_EQ(lsm.isFaceTreated(f0_id), false);
	lsm.setFaceTreated(f0_id);
	ASSERT_EQ(lsm.isFaceTreated(f0_id), true);

	TCellID e_id = math::Utils::CommonEdge(&m, n0.id(), n1.id());
	ASSERT_EQ(lsm.isEdgeTreated(e_id), false);
	lsm.setEdgeTreated(e_id);
	ASSERT_EQ(lsm.isEdgeTreated(e_id), true);

	ASSERT_EQ(lsm.getFaceIdealNextNode(f0_id, n0.id()), map_new_nodes[n0.id()]);
	lsm.setFaceNextNode(f0_id, n0.id(), n1.id());
	ASSERT_EQ(lsm.getFaceNextNode(f0_id, n0.id()), n1.id());
	lsm.setFaceNextNode(f0_id, n0.id(), map_new_nodes[n0.id()]);

	// Create fake nodes to test structure
	Node n_8_5 = m.newNode(n8.point());
	Node n_8_8 = m.newNode(n8.point());
	Node n_14_5 = m.newNode(n14.point());
	Node n_14_8 = m.newNode(n14.point());

	// Test END structure
	ASSERT_EQ(lsm.isEndFaceCreated(e_end_id, n6.id()), false);
	lsm.setEndFaceCreated(e_end_id, n6.id());
	ASSERT_EQ(lsm.isEndFaceCreated(e_end_id, n6.id()), true);

	// Test CORNER structure
	ASSERT_EQ(lsm.isCornerFaceCreated(e_corner_id, n8.id()), false);
	ASSERT_EQ(lsm.isCornerFaceCreated(e_corner_id, n14.id()), false);
	lsm.setCornerFaceCreated(e_corner_id, n8.id());
	lsm.setCornerFaceCreated(e_corner_id, n14.id());
	ASSERT_EQ(lsm.isCornerFaceCreated(e_corner_id, n8.id()), true);
	ASSERT_EQ(lsm.isCornerFaceCreated(e_corner_id, n14.id()), true);
	lsm.setCornerNextAdjNode(e_corner_id, f5_id, n8.id(), n_8_5.id());
	lsm.setCornerNextAdjNode(e_corner_id, f8_id, n8.id(), n_8_8.id());
	ASSERT_EQ(lsm.getCornerNextAdjNode(e_corner_id, f5_id, n8.id()), n_8_5.id());
	ASSERT_EQ(lsm.getCornerNextAdjNode(e_corner_id, f8_id, n8.id()), n_8_8.id());
	lsm.setNextDiagNode(e_corner_id, n8.id(), lsm.getFaceIdealNextNode(f8_id, n8.id()));
	lsm.setNextDiagNode(e_corner_id, n14.id(), lsm.getFaceIdealNextNode(f8_id, n14.id()));
	ASSERT_EQ(lsm.getNextDiagNode(e_corner_id, n8.id()), map_new_nodes[n8.id()]);
	ASSERT_EQ(lsm.getNextDiagNode(e_corner_id, n14.id()), map_new_nodes[n14.id()]);

	// Test REVERSAL structure
	// Create fake nodes to test data structure
	Node n_adj_0_1 = m.newNode(n0.point());	// Adjacent node of n0 in face f1 side
	Node n_adj_1_1 = m.newNode(n1.point());
	Node n_adj_0_0 = m.newNode(n0.point());
	Node n_adj_1_0 = m.newNode(n1.point());
	Node n_diag_0_1 = m.newNode(n0.point());	// Diag node of node n0 in face f1 side
	Node n_diag_1_1 = m.newNode(n1.point());
	Node n_diag_0_0 = m.newNode(n0.point());
	Node n_diag_1_0 = m.newNode(n1.point());

	ASSERT_EQ(lsm.areReversalFacesCreated(e_reversal_id, n0.id()), false);
	ASSERT_EQ(lsm.areReversalFacesCreated(e_reversal_id, n1.id()), false);
	lsm.setReversalFacesCreated(e_reversal_id, n0.id());
	lsm.setReversalFacesCreated(e_reversal_id, n1.id());
	ASSERT_EQ(lsm.areReversalFacesCreated(e_reversal_id, n0.id()), true);
	ASSERT_EQ(lsm.areReversalFacesCreated(e_reversal_id, n1.id()), true);
	lsm.setReversalNextAdjNode(e_reversal_id, f0_id, n0.id(), n_adj_0_0.id());
	lsm.setReversalNextAdjNode(e_reversal_id, f0_id, n1.id(), n_adj_1_0.id());
	lsm.setReversalNextAdjNode(e_reversal_id, f1_id, n0.id(), n_adj_0_1.id());
	lsm.setReversalNextAdjNode(e_reversal_id, f1_id, n1.id(), n_adj_1_1.id());
	ASSERT_EQ(lsm.getReversalNextAdjNode(e_reversal_id, f0_id, n0.id()), n_adj_0_0.id());
	ASSERT_EQ(lsm.getReversalNextAdjNode(e_reversal_id, f0_id, n1.id()), n_adj_1_0.id());
	ASSERT_EQ(lsm.getReversalNextAdjNode(e_reversal_id, f1_id, n0.id()), n_adj_0_1.id());
	ASSERT_EQ(lsm.getReversalNextAdjNode(e_reversal_id, f1_id, n1.id()), n_adj_1_1.id());
	lsm.setReversalNextMediumNode(e_reversal_id, n0.id(), map_new_nodes[n0.id()]);
	lsm.setReversalNextMediumNode(e_reversal_id, n1.id(), map_new_nodes[n1.id()]);
	ASSERT_EQ(lsm.getReversalNextMediumNode(e_reversal_id, n0.id()), map_new_nodes[n0.id()]);
	ASSERT_EQ(lsm.getReversalNextMediumNode(e_reversal_id, n1.id()), map_new_nodes[n1.id()]);

}
/*----------------------------------------------------------------------------*/