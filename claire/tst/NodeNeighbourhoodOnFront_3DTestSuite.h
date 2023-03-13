//
// Created by rochec on 21/12/2022.
//

#include <gmds/ig/Mesh.h>
#include <gmds/claire/NodeNeighbourhoodOnFront_3D.h>
#include <gmds/claire/Utils.h>
#include <gmds/claire/Front_3D.h>
#include <gtest/gtest.h>
#include <iostream>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

TEST(NodeNeighbourhoodOnFront_3DTestClass, NodeNeighbourhoodOnFront_3D_Test1)
{
	Mesh* m_mesh = new Mesh(gmds::MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
	                                     F2E | E2F | R2E | E2R | N2R | N2F | N2E));
	   
	Node n0 = m_mesh->newNode({0, 0, 0});

	Node n1 = m_mesh->newNode({-0.5, -0.5, 0.5});
	Node n2 = m_mesh->newNode({0.5, -0.5, 0.5});
	Node n3 = m_mesh->newNode({0.0, 0.0, 0.5});
	Node n4 = m_mesh->newNode({0.5, 0.0, 0.5});
	Node n5 = m_mesh->newNode({-0.5, 0.5, 0.5});
	Node n6 = m_mesh->newNode({0.0, 0.5, 0.5});

	Node n7 = m_mesh->newNode({0.5, 0.0, 0.0});
	Node n8 = m_mesh->newNode({0.0, 0.5, 0.0});
	Node n9 = m_mesh->newNode({0.5, 0.5, 0.0});

	Node n10 = m_mesh->newNode({-0.5, -0.5, -0.5});
	Node n11 = m_mesh->newNode({0.5, -0.5, -0.5});
	Node n12 = m_mesh->newNode({-0.5, 0.5, -0.5});
	Node n13 = m_mesh->newNode({0.5, 0.5, -0.5});

	Node n20 = m_mesh->newNode({0, -0.5, 0.5});
	Node n21 = m_mesh->newNode({-0.5, -0.5, 0});
	Node n22 = m_mesh->newNode({0.0, -0.5, 0});
	Node n23 = m_mesh->newNode({0.5, -0.5, 0});
	Node n24 = m_mesh->newNode({0.0, -0.5, -0.5});
	Node n25 = m_mesh->newNode({0.25, 0.0, 0.5});
	Node n26 = m_mesh->newNode({0.0, 0.0, 0.25});
	Node n27 = m_mesh->newNode({0.25, 0.0, 0.25});
	Node n28 = m_mesh->newNode({0.5, 0.0, 0.25});
	Node n29 = m_mesh->newNode({0.25, 0.0, 0.0});
	Node n30 = m_mesh->newNode({-0.5, 0.0, 0.5});
	Node n31 = m_mesh->newNode({0.0, 0.25, 0.5});
	Node n32 = m_mesh->newNode({-0.5, 0.0, 0.0});
	Node n33 = m_mesh->newNode({0.0, 0.25, 0.25});
	Node n34 = m_mesh->newNode({0.0, 0.25, 0.0});
	Node n35 = m_mesh->newNode({0.25, 0.25, 0.0});
	Node n36 = m_mesh->newNode({0.5, 0.25, 0.0});
	Node n37 = m_mesh->newNode({0.5, 0.0, -0.5});
	Node n38 = m_mesh->newNode({-0.5, 0.0, -0.5});
	Node n39 = m_mesh->newNode({0.0, 0.0, -0.5});
	Node n40 = m_mesh->newNode({0.0, 0.5, 0.25});
	Node n41 = m_mesh->newNode({-0.5, 0.5, 0.0});
	Node n42 = m_mesh->newNode({0.25, 0.5, 0.0});
	Node n43 = m_mesh->newNode({0.0, 0.5, -0.5});

	// Create the faces
	TCellID f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n1.id(), n3.id(), n31.id(), n30.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n30.id(), n31.id(), n6.id(), n5.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n1.id(), n20.id(), n25.id(), n3.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n20.id(), n2.id(), n4.id(), n25.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n0.id(), n29.id(), n35.id(), n34.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n29.id(), n7.id(), n36.id(), n35.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n34.id(), n35.id(), n42.id(), n8.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n35.id(), n36.id(), n9.id(), n42.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n10.id(), n24.id(), n39.id(), n38.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n24.id(), n11.id(), n37.id(), n39.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n12.id(), n43.id(), n39.id(), n38.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n39.id(), n37.id(), n13.id(), n43.id());

	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n2.id(), n4.id(), n28.id(), n23.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n11.id(), n7.id(), n28.id(), n23.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n11.id(), n7.id(), n36.id(), n37.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n13.id(), n9.id(), n36.id(), n37.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n3.id(), n31.id(), n33.id(), n26.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n0.id(), n34.id(), n33.id(), n26.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n40.id(), n33.id(), n34.id(), n8.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n40.id(), n33.id(), n31.id(), n6.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n1.id(), n30.id(), n32.id(), n21.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n10.id(), n38.id(), n32.id(), n21.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n12.id(), n38.id(), n32.id(), n41.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n5.id(), n30.id(), n32.id(), n41.id());

	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n5.id(), n6.id(), n40.id(), n41.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n12.id(), n8.id(), n40.id(), n41.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n12.id(), n8.id(), n42.id(), n43.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n13.id(), n9.id(), n42.id(), n43.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n0.id(), n29.id(), n27.id(), n26.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n3.id(), n25.id(), n27.id(), n26.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n28.id(), n27.id(), n25.id(), n4.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n28.id(), n27.id(), n29.id(), n7.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n24.id(), n22.id(), n21.id(), n10.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n24.id(), n22.id(), n23.id(), n11.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n20.id(), n22.id(), n23.id(), n2.id());
	f_id = math::Utils::GetOrCreateQuadAndConnectivities(m_mesh, n20.id(), n22.id(), n21.id(), n1.id());

	Variable<int>* couche_id = m_mesh->newVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* couche_face_id = m_mesh->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");

	std::vector<TCellID> front_nodes;
	std::vector<TCellID> front_faces;

	for (auto n_id:m_mesh->nodes())
	{
		couche_id->set(n_id, 0);
		front_nodes.push_back(n_id);
	}
	for (auto f_loc_id:m_mesh->faces())
	{
		couche_face_id->set(f_loc_id, 0);
		front_faces.push_back(f_loc_id);
	}

	Front_3D Front = Front_3D(0, front_nodes, front_faces);

	NodeNeighbourhoodOnFront_3D Node_Neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, &Front, n12.id());
	NodeNeighbourhoodOnFront_3D::STATUS result = Node_Neighbourhood.execute();

	std::vector<TCellID> orderedEdges = Node_Neighbourhood.getOrderedEdges() ;
	ASSERT_EQ(orderedEdges[0], 32);
	ASSERT_EQ(orderedEdges[1], 62);
	ASSERT_EQ(orderedEdges[2], 59);
	ASSERT_EQ(orderedEdges[3], 34);

	std::vector<TCellID> orderedFaces = Node_Neighbourhood.getOrderedFaces() ;
	ASSERT_EQ(orderedFaces[0], 26);
	ASSERT_EQ(orderedFaces[1], 25);
	ASSERT_EQ(orderedFaces[2], 22);
	ASSERT_EQ(orderedFaces[3], 10);

	std::vector<TCellID> adj_faces = Node_Neighbourhood.adjFacesToEdge(32);
	ASSERT_EQ(adj_faces[0], 10);
	ASSERT_EQ(adj_faces[1], 26);

	adj_faces = Node_Neighbourhood.adjFacesToEdge(59);
	ASSERT_EQ(adj_faces[0], 25);
	ASSERT_EQ(adj_faces[1], 22);

	TCellID next_f_id = Node_Neighbourhood.nextEdgeOfFace(10, 32);
	ASSERT_EQ(next_f_id, 34);

	next_f_id = Node_Neighbourhood.nextEdgeOfFace(10, 34);
	ASSERT_EQ(next_f_id, 32);

	next_f_id = Node_Neighbourhood.nextEdgeOfFace(22, 59);
	ASSERT_EQ(next_f_id, 34);

	next_f_id = Node_Neighbourhood.nextEdgeOfFace(26, 32);
	ASSERT_EQ(next_f_id, 62);

	TCellID face_adj = Node_Neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(32, 59, 62);
	ASSERT_EQ(face_adj, 10);

	face_adj = Node_Neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(32, 34, 62);
	ASSERT_EQ(face_adj, 10);

	face_adj = Node_Neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(32, 62, 34);
	ASSERT_EQ(face_adj, 26);

	face_adj = Node_Neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(32, 62, 59);
	ASSERT_EQ(face_adj, 26);

	std::vector<TCellID> faces_btw_2edges = Node_Neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(32,59,62);
	ASSERT_EQ(faces_btw_2edges.size(), 2);
	ASSERT_EQ(faces_btw_2edges[0], 10);
	ASSERT_EQ(faces_btw_2edges[1], 22);

	faces_btw_2edges = Node_Neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(62,34,59);
	ASSERT_EQ(faces_btw_2edges.size(), 2);
	ASSERT_EQ(faces_btw_2edges[0], 26);
	ASSERT_EQ(faces_btw_2edges[1], 10);

	faces_btw_2edges = Node_Neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(62,59,32);
	ASSERT_EQ(faces_btw_2edges.size(), 1);
	ASSERT_EQ(faces_btw_2edges[0], 25);

	gmds::IGMeshIOService ioService(m_mesh);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	vtkWriter.write("TEST_SUITE_NodeNeighbourhoodOnFront_3D.vtk");

	ASSERT_EQ(NodeNeighbourhoodOnFront_3D::SUCCESS, result);

}
