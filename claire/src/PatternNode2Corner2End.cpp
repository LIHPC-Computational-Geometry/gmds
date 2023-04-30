//
// Created by rochec on 28/04/23.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/PatternNode2Corner2End.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/claire/NodeNeighbourhoodOnFront_3D.h>
#include <gmds/claire/Utils.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

PatternNode2Corner2End::PatternNode2Corner2End(Mesh *AMesh, Front_3D *AFront, TCellID An_id, LayerStructureManager_3D *AStructManager,
                                       Mesh *AMeshT, FastLocalize *Afl,
                                       double dc, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField) :
  AbstractPatternNode(AMesh, AFront, An_id, AStructManager, AMeshT, Afl, dc, A_DistanceField, A_VectorField)
{

}
/*------------------------------------------------------------------------*/
void
PatternNode2Corner2End::computeNewHex()
{
	Node n = m_mesh->get<Node>(m_n_id);

	NodeNeighbourhoodOnFront_3D n_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, m_n_id);
	n_neighbourhood.execute();
	std::vector<TCellID> n_ordered_edges = n_neighbourhood.getOrderedEdges();

	Variable<int>* var_front_edges_classification = m_mesh->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
	std::vector<Edge> ec;
	std::vector<Edge> ee;
	for (auto e_loc_id:n_ordered_edges)
	{
		if (var_front_edges_classification->value(e_loc_id)==1)
		{
			ec.push_back(m_mesh->get<Edge>(e_loc_id));
		}
	}

	// Store the edges in ee vector in the same side of the ones in ec
	// This way, ee[0] will be the closest end edge to ec[0], and ee[1] the closest to ec[1]
	// ---->
	std::vector<TCellID> ec_0_faces = n_neighbourhood.adjFacesToEdge(ec[0].id());
	if (var_front_edges_classification->value(n_neighbourhood.nextEdgeOfFace(ec_0_faces[0], ec[0].id()))==2)
	{
		ee.push_back(m_mesh->get<Edge>(n_neighbourhood.nextEdgeOfFace(ec_0_faces[0], ec[0].id())));
	}
	else
	{
		ee.push_back(m_mesh->get<Edge>(n_neighbourhood.nextEdgeOfFace(ec_0_faces[1], ec[0].id())));
	}

	std::vector<TCellID> ec_1_faces = n_neighbourhood.adjFacesToEdge(ec[1].id());
	if (var_front_edges_classification->value(n_neighbourhood.nextEdgeOfFace(ec_1_faces[0], ec[1].id()))==2)
	{
		ee.push_back(m_mesh->get<Edge>(n_neighbourhood.nextEdgeOfFace(ec_1_faces[0], ec[1].id())));
	}
	else
	{
		ee.push_back(m_mesh->get<Edge>(n_neighbourhood.nextEdgeOfFace(ec_1_faces[1], ec[1].id())));
	}
	// <----

	Node nc0 = ec[0].getOppositeNode(n);
	Node nc1 = ec[1].getOppositeNode(n);

	// Get the 5 faces around the node n
	TCellID f0_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec[0].id(), ee[0].id(), ee[1].id());
	TCellID f1_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec[0].id(), ee[1].id(), ee[0].id());
	TCellID f2_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ee[1].id(), ec[0].id(), ec[1].id());
	TCellID f3_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ee[1].id(), ec[1].id(), ec[0].id());
	TCellID f4_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec[1].id(), ee[0].id(), ee[1].id());
	TCellID f5_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ee[0].id(), ec[1].id(), ec[0].id());

	// Get the 2 other nodes of Face f1 (the first ones are n and nc0) ---->
	TCellID f1_next_edge = n_neighbourhood.nextEdgeOfFace(f1_id, ec[0].id());
	TCellID n_f1_1_id = (m_mesh->get<Edge>(f1_next_edge)).getOppositeNodeId(m_n_id) ;
	Node n_f1_1 = m_mesh->get<Node>(n_f1_1_id);

	Edge f1_opp_edge = math::Utils::oppositeEdgeInFace(m_mesh, ec[0].id(), f1_id);
	Node n_f1_2 = f1_opp_edge.getOppositeNode(n_f1_1);
	// <----

	// Get the 2 other nodes of Face f4 (the first ones are n and nc1) ---->
	TCellID f4_next_edge = n_neighbourhood.nextEdgeOfFace(f4_id, ec[1].id());
	TCellID n_f4_1_id = (m_mesh->get<Edge>(f4_next_edge)).getOppositeNodeId(m_n_id) ;
	Node n_f4_1 = m_mesh->get<Node>(n_f4_1_id);

	Edge f4_opp_edge = math::Utils::oppositeEdgeInFace(m_mesh, ec[1].id(), f4_id);
	Node n_f4_2 = f4_opp_edge.getOppositeNode(n_f4_1);
	// <----

	Node n1 = m_mesh->get<Node>(m_StructManager->getFaceIdealNextNode(f0_id, m_n_id));
	Node n2 = m_mesh->get<Node>(m_StructManager->getFaceIdealNextNode(f1_id, n_f1_2.id()));
	Node n3 = m_mesh->get<Node>(m_StructManager->getFaceIdealNextNode(f1_id, n_f1_1.id()));
	Node n4 = m_mesh->get<Node>(m_StructManager->getFaceIdealNextNode(f4_id, n_f4_2.id()));
	Node n5 = m_mesh->get<Node>(m_StructManager->getFaceIdealNextNode(f4_id, n_f4_1.id()));

	// Create the two new hexas
	TCellID r_id = math::Utils::CreateHexaNConnectivities(m_mesh, n, nc0, n1, nc1, n_f1_1, n_f1_2, n2, n3);
	m_hex.push_back(m_mesh->get<Region>(r_id));
	r_id = math::Utils::CreateHexaNConnectivities(m_mesh, n, nc0, n1, nc1, n_f4_1, n5, n4, n_f4_2);
	m_hex.push_back(m_mesh->get<Region>(r_id));


	//===================//
	//      UPDATES      //
	//===================//

	// Update the faces and edges that can't create hexas anymore
	m_StructManager->setEdgeTreated(ec[0].id());
	m_StructManager->setEdgeTreated(ec[1].id());
	m_StructManager->setFaceTreated(f1_id);
	m_StructManager->setFaceTreated(f4_id);

	// Update the edge struc info for the 2 END edges ---->
	m_StructManager->setEndFaceCreated(ee[0].id(), m_n_id);
	m_StructManager->setEndFaceCreated(ee[1].id(), m_n_id);
	m_StructManager->setNextDiagNode(ee[0].id(), m_n_id, n5.id());
	m_StructManager->setNextDiagNode(ee[1].id(), m_n_id, n3.id());
	// <----

	// Update the edge struc info for the first CORNER edge ---->
	NodeNeighbourhoodOnFront_3D n_c0_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, nc0.id()) ;
	n_c0_neighbourhood.execute();
	for (auto edge_id:n_c0_neighbourhood.getOrderedEdges()) {
		if (var_front_edges_classification->value(edge_id) == 1 && edge_id != ec[0].id()) {
			m_StructManager->setCornerFaceCreated(edge_id, nc0.id());
			m_StructManager->setNextDiagNode(edge_id, nc0.id(), n4.id());

			TCellID edge_f1_last = n_c0_neighbourhood.nextEdgeOfFace(f1_id, ec[0].id());
			for (auto faces_id : n_c0_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[0].id(), edge_id, edge_f1_last)) {
				m_StructManager->setFaceNextNode(faces_id, nc0.id(), n5.id());
			}
			TCellID edge_f0 = n_c0_neighbourhood.nextEdgeOfFace(f0_id, ec[0].id());
			for (auto faces_id : n_c0_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[0].id(), edge_id, edge_f0)) {
				m_StructManager->setFaceNextNode(faces_id, nc0.id(), n1.id());
			}

			std::pair<TCellID, TCellID> pair_1(n_c0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[0].id(), edge_f1_last), nc0.id());
			std::pair<TCellID, TCellID> pair_2(n_c0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[0].id(), edge_f0), nc0.id());

			m_StructManager->setCornerNextAdjNode(edge_id, n_c0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[0].id(), edge_f1_last), nc0.id(), n5.id());
			m_StructManager->setCornerNextAdjNode(edge_id, n_c0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[0].id(), edge_f0), nc0.id(), n1.id());
		}
	}
	// <----

	// Update the edge struc info for the second CORNER edge ---->
	NodeNeighbourhoodOnFront_3D n_c1_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, nc1.id()) ;
	n_c1_neighbourhood.execute();
	for (auto edge_id:n_c1_neighbourhood.getOrderedEdges()) {
		if (var_front_edges_classification->value(edge_id) == 1 && edge_id != ec[1].id()) {
			m_StructManager->setCornerFaceCreated(edge_id, nc1.id());
			m_StructManager->setNextDiagNode(edge_id, nc1.id(), n2.id());

			TCellID edge_f4_last = n_c1_neighbourhood.nextEdgeOfFace(f4_id, ec[1].id());
			for (auto faces_id : n_c1_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[1].id(), edge_id, edge_f4_last)) {
				m_StructManager->setFaceNextNode(faces_id, nc1.id(), n3.id());
			}
			TCellID edge_f3 = n_c1_neighbourhood.nextEdgeOfFace(f3_id, ec[1].id());
			for (auto faces_id : n_c1_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[1].id(), edge_id, edge_f3)) {
				m_StructManager->setFaceNextNode(faces_id, nc1.id(), n1.id());
			}

			std::pair<TCellID, TCellID> pair_1(n_c1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[1].id(), edge_f4_last), nc1.id());
			std::pair<TCellID, TCellID> pair_2(n_c1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[1].id(), edge_f3), nc1.id());

			m_StructManager->setCornerNextAdjNode(edge_id, n_c1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[1].id(), edge_f4_last), nc1.id(), n3.id());
			m_StructManager->setCornerNextAdjNode(edge_id, n_c1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[1].id(), edge_f3), nc1.id(), n1.id());
		}
	}
	// <----
}