//
// Created by rochec on 12/05/23.
//

#include <gmds/aero/Front_3D.h>
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

TEST(Front_3DTestClass, Front_3D)
{
	// WE READ
	gmds::Mesh m(gmds::MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                             F2E | E2F | R2E | E2R | N2R | N2F | N2E));

	Node n0 = m.newNode({-0.5, -0.5, -0.5});
	Node n1 = m.newNode({-0.5, 0.5, -0.5});
	Node n2 = m.newNode({0.5, 0.5, -0.5});
	Node n3 = m.newNode({0.5, -0.5, -0.5});

	Node n4 = m.newNode({-0.5, -0.5, 0.5});
	Node n5 = m.newNode({-0.5, 0.5, 0.5});
	Node n6 = m.newNode({0.5, 0.5, 0.5});
	Node n7 = m.newNode({0.5, -0.5, 0.5});

	// Creates the faces
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
}
