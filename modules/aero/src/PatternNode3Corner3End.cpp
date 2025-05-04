//
// Created by rochec on 28/04/23.
//
/*------------------------------------------------------------------------*/
#include <gmds/aero/PatternNode3Corner3End.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/aero/NodeNeighbourhoodOnFront_3D.h>
#include <gmds/aero/Utils.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

PatternNode3Corner3End::PatternNode3Corner3End(Mesh *AMesh, Front_3D *AFront, TCellID An_id, LayerStructureManager_3D *AStructManager,
                                       Mesh *AMeshT, FastLocalize *Afl,
                                       double dc, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField) :
  AbstractPatternNode(AMesh, AFront, An_id, AStructManager, AMeshT, Afl, dc, A_DistanceField, A_VectorField)
{

}
/*------------------------------------------------------------------------*/
void
PatternNode3Corner3End::computeNewHex()
{
	Node n = m_mesh->get<Node>(m_n_id);

	NodeNeighbourhoodOnFront_3D n_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, m_n_id);
	n_neighbourhood.execute();
	std::vector<TCellID> n_ordered_edges = n_neighbourhood.getOrderedEdges();

	Variable<int>* var_front_edges_classification = m_mesh->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
	std::vector<Edge> ec;
	std::vector<Edge> ee;
	if (var_front_edges_classification->value(n_ordered_edges[0]) == 1)
	{
		ec.push_back(m_mesh->get<Edge>(n_ordered_edges[0]));
		ec.push_back(m_mesh->get<Edge>(n_ordered_edges[2]));
		ec.push_back(m_mesh->get<Edge>(n_ordered_edges[4]));
		ee.push_back(m_mesh->get<Edge>(n_ordered_edges[1]));
		ee.push_back(m_mesh->get<Edge>(n_ordered_edges[3]));
		ee.push_back(m_mesh->get<Edge>(n_ordered_edges[5]));
	}
	else if (var_front_edges_classification->value(n_ordered_edges[0])==2)
	{
		ee.push_back(m_mesh->get<Edge>(n_ordered_edges[0]));
		ee.push_back(m_mesh->get<Edge>(n_ordered_edges[2]));
		ee.push_back(m_mesh->get<Edge>(n_ordered_edges[4]));
		ec.push_back(m_mesh->get<Edge>(n_ordered_edges[5]));
		ec.push_back(m_mesh->get<Edge>(n_ordered_edges[1]));
		ec.push_back(m_mesh->get<Edge>(n_ordered_edges[3]));
	}

	Node nc0 = ec[0].getOppositeNode(n);
	Node nc1 = ec[1].getOppositeNode(n);
	Node nc2 = ec[2].getOppositeNode(n);

	Node n1 = m_mesh->get<Node>(m_StructManager->getFaceIdealNextNode(n.get<Face>()[0].id(), m_n_id));

	math::Vector3d v01 = ( (nc0.point() - n.point()).normalize() + (nc1.point() - n.point()).normalize() + (n1.point() - n.point()).normalize() ).normalize() ;
	math::Vector3d v12 = ( (nc1.point() - n.point()).normalize() + (nc2.point() - n.point()).normalize() + (n1.point() - n.point()).normalize() ).normalize() ;
	math::Vector3d v02 = ( (nc0.point() - n.point()).normalize() + (nc2.point() - n.point()).normalize() + (n1.point() - n.point()).normalize() ).normalize() ;

	math::Point pe0 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, v01);
	math::Point pe1 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, v12);
	math::Point pe2 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, v02);

	Node ne0 = m_mesh->newNode(pe0);
	Node ne1 = m_mesh->newNode(pe1);
	Node ne2 = m_mesh->newNode(pe2);

	// Create the new hexa
	m_hex.push_back(m_mesh->get<Region>( math::Utils::CreateHexaNConnectivities(m_mesh, n, nc0, ne0, nc1, nc2, ne2, n1, ne1)));

	// Update the layer ids ---->
	Variable<int>* var_node_couche_id = m_mesh->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	var_node_couche_id->set(ne0.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(ne1.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(ne2.id(), m_Front->getFrontID()+1);
	// <----

	// Update the edge struc info for the 3 end edges ---->
	m_StructManager->setEndFaceCreated(ee[0].id(), m_n_id);
	m_StructManager->setEndFaceCreated(ee[1].id(), m_n_id);
	m_StructManager->setEndFaceCreated(ee[2].id(), m_n_id);
	m_StructManager->setNextDiagNode(ee[0].id(), m_n_id, ne0.id());
	m_StructManager->setNextDiagNode(ee[1].id(), m_n_id, ne1.id());
	m_StructManager->setNextDiagNode(ee[2].id(), m_n_id, ne2.id());
	// <----

	TCellID f0_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec[0].id(), ee[0].id(), ee[2].id());
	TCellID f1_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ee[0].id(), ec[1].id(), ec[0].id());
	TCellID f2_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec[1].id(), ee[1].id(), ee[0].id());
	TCellID f3_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ee[1].id(), ec[2].id(), ec[1].id());
	TCellID f4_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec[2].id(), ee[2].id(), ee[1].id());
	TCellID f5_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ee[2].id(), ec[0].id(), ec[2].id());

	// Update the edge struc info for the first corner edge ---->
	NodeNeighbourhoodOnFront_3D n_c0_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, nc0.id()) ;
	n_c0_neighbourhood.execute();
	for (auto edge_id:n_c0_neighbourhood.getOrderedEdges()) {
		if (var_front_edges_classification->value(edge_id) == 1 && edge_id != ec[0].id()) {
			m_StructManager->setCornerFaceCreated(edge_id, nc0.id());
			m_StructManager->setNextDiagNode(edge_id, nc0.id(), n1.id());

			TCellID edge_opp_f0_e0 = n_c0_neighbourhood.nextEdgeOfFace(f0_id, ec[0].id());
			for (auto faces_id : n_c0_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[0].id(), edge_id, edge_opp_f0_e0)) {
				m_StructManager->setFaceNextNode(faces_id, nc0.id(), ne2.id());
			}
			TCellID edge_opp_f5_ee2 = n_c0_neighbourhood.nextEdgeOfFace(f5_id, ec[0].id());
			for (auto faces_id : n_c0_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[0].id(), edge_id, edge_opp_f5_ee2)) {
				m_StructManager->setFaceNextNode(faces_id, nc0.id(), ne0.id());
			}

			std::pair<TCellID, TCellID> pair_1(n_c0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[0].id(), edge_opp_f5_ee2), nc0.id());
			std::pair<TCellID, TCellID> pair_2(n_c0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[0].id(), edge_opp_f0_e0), nc0.id());

			m_StructManager->setCornerNextAdjNode(edge_id, n_c0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[0].id(), edge_opp_f5_ee2), nc0.id(), ne0.id());
			m_StructManager->setCornerNextAdjNode(edge_id, n_c0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[0].id(), edge_opp_f0_e0), nc0.id(), ne2.id());
		}
	}
	// <----

	// Update the edge struc info for the second corner edge ---->
	NodeNeighbourhoodOnFront_3D n_c1_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, nc1.id()) ;
	n_c1_neighbourhood.execute();
	for (auto edge_id:n_c1_neighbourhood.getOrderedEdges()) {
		if (var_front_edges_classification->value(edge_id) == 1 && edge_id != ec[1].id()) {
			m_StructManager->setCornerFaceCreated(edge_id, nc1.id());
			m_StructManager->setNextDiagNode(edge_id, nc1.id(), n1.id());

			TCellID edge_opp_f1_e0 = n_c1_neighbourhood.nextEdgeOfFace(f1_id, ec[1].id());
			for (auto faces_id : n_c1_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[1].id(), edge_id, edge_opp_f1_e0)) {
				m_StructManager->setFaceNextNode(faces_id, nc1.id(), ne0.id());
			}
			TCellID edge_opp_f2_ee1 = n_c1_neighbourhood.nextEdgeOfFace(f2_id, ec[1].id());
			for (auto faces_id : n_c1_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[1].id(), edge_id, edge_opp_f2_ee1)) {
				m_StructManager->setFaceNextNode(faces_id, nc0.id(), ne0.id());
			}

			std::pair<TCellID, TCellID> pair_1(n_c1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[1].id(), edge_opp_f2_ee1), nc1.id());
			std::pair<TCellID, TCellID> pair_2(n_c1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[1].id(), edge_opp_f1_e0), nc1.id());

			m_StructManager->setCornerNextAdjNode(edge_id, n_c1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[1].id(), edge_opp_f2_ee1), nc1.id(), ne0.id());
			m_StructManager->setCornerNextAdjNode(edge_id, n_c1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[1].id(), edge_opp_f1_e0), nc1.id(), ne1.id());
		}
	}
	// <----

	// Update the edge struc info for the thirst corner edge ---->
	NodeNeighbourhoodOnFront_3D n_c2_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, nc2.id()) ;
	n_c2_neighbourhood.execute();
	for (auto edge_id:n_c2_neighbourhood.getOrderedEdges()) {
		if (var_front_edges_classification->value(edge_id) == 1 && edge_id != ec[2].id()) {
			m_StructManager->setCornerFaceCreated(edge_id, nc2.id());
			m_StructManager->setNextDiagNode(edge_id, nc2.id(), n1.id());

			TCellID edge_opp_f3_ee1 = n_c2_neighbourhood.nextEdgeOfFace(f3_id, ec[2].id());
			for (auto faces_id : n_c2_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[2].id(), edge_id, edge_opp_f3_ee1)) {
				m_StructManager->setFaceNextNode(faces_id, nc2.id(), ne2.id());
			}
			TCellID edge_opp_f4_ee2 = n_c2_neighbourhood.nextEdgeOfFace(f4_id, ec[2].id());
			for (auto faces_id : n_c2_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[2].id(), edge_id, edge_opp_f4_ee2)) {
				m_StructManager->setFaceNextNode(faces_id, nc2.id(), ne1.id());
			}

			std::pair<TCellID, TCellID> pair_1(n_c2_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[2].id(), edge_opp_f4_ee2), nc2.id());
			std::pair<TCellID, TCellID> pair_2(n_c2_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[2].id(), edge_opp_f3_ee1), nc2.id());

			m_StructManager->setCornerNextAdjNode(edge_id, n_c2_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[2].id(), edge_opp_f4_ee2), nc2.id(), ne1.id());
			m_StructManager->setCornerNextAdjNode(edge_id, n_c2_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[2].id(), edge_opp_f3_ee1), nc2.id(), ne2.id());
		}
	}
	// <----

	// Update the faces and edges that can't do actions anymore
	m_StructManager->setEdgeTreated(ec[0].id());
	m_StructManager->setEdgeTreated(ec[1].id());
	m_StructManager->setEdgeTreated(ec[2].id());
}