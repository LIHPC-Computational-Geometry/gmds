//
// Created by rochec on 11/10/2023.
//
#include <gmds/claire/LayerStructureManager_3D.h>
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

/*----------------------------------------------------------------------------*/
TEST(ClaireTestClass, LayerStructureManager_3D)
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
	}

	// Create Structure Manager
	LayerStructureManager_3D lsm = LayerStructureManager_3D(&m, &F0, map_new_nodes);

	ASSERT_EQ(lsm.isFaceTreated(f0_id), false);

	ASSERT_EQ(0, 0);

}
/*----------------------------------------------------------------------------*/