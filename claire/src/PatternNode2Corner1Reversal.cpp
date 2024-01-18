//
// Created by rochec on 28/04/23.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/PatternNode2Corner1Reversal.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/claire/NodeNeighbourhoodOnFront_3D.h>
#include <gmds/claire/Utils.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

PatternNode2Corner1Reversal::PatternNode2Corner1Reversal(Mesh *AMesh, Front_3D *AFront, TCellID An_id, LayerStructureManager_3D *AStructManager,
                                       Mesh *AMeshT, FastLocalize *Afl,
                                       double dc, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField) :
  AbstractPatternNode(AMesh, AFront, An_id, AStructManager, AMeshT, Afl, dc, A_DistanceField, A_VectorField)
{

}
/*------------------------------------------------------------------------*/
void
PatternNode2Corner1Reversal::computeNewHex()
{
	Node n = m_mesh->get<Node>(m_n_id);
	NodeNeighbourhoodOnFront_3D n_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, m_n_id);
	n_neighbourhood.execute();
	std::vector<TCellID> n_ordered_edges = n_neighbourhood.getOrderedEdges();
	std::vector<TCellID> n_ordered_faces = n_neighbourhood.getOrderedFaces();

	Variable<int>* var_front_edges_classification = m_mesh->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
	std::vector<Edge> ec;
	Edge er;
	for (auto e_loc_id:n_ordered_edges)
	{
		if (var_front_edges_classification->value(e_loc_id)==1)
		{
			ec.push_back(m_mesh->get<Edge>(e_loc_id));
		}
		else if (var_front_edges_classification->value(e_loc_id)==3)
		{
			er = m_mesh->get<Edge>(e_loc_id);
		}
	}

	Node n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11;
	n2 = m_mesh->get<Node>(m_StructManager->getFaceIdealNextNode(n_ordered_faces[0], m_n_id));

	// We build the 10 missing nodes ---->
	Node nr = er.getOppositeNode(n);
	Node nc0 = ec[0].getOppositeNode(n);
	Node nc1 = ec[1].getOppositeNode(n);

	math::Vector3d vr = (n.point()-nr.point()).normalize() ;
	math::Vector3d vc0 = (n.point()-nc0.point()).normalize() ;
	math::Vector3d vc1 = (n.point()-nc1.point()).normalize() ;

	Face f_rc0 = m_mesh->get<Face>(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(er.id(), ec[0].id(), ec[1].id()));
	Face f_rc1 = m_mesh->get<Face>(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(er.id(), ec[1].id(), ec[0].id()));
	Face f_c0c1 = m_mesh->get<Face>(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec[0].id(), ec[1].id(), er.id()));

	/*
	math::Vector3d v_rc0 = f_rc0.normal();
	gmds::Cell::Data data_rc0 = m_fl->find(f_rc0.center());
	Node n_closest_f_rc0 = m_meshT->get<Node>(data_rc0.id);
	math::Vector3d v_test = m_VectorField->value(n_closest_f_rc0.id()).normalize() ;
	if (v_rc0.dot(v_test) <= 0)
	{
		v_rc0 = - v_rc0;
	}

	math::Vector3d v_rc1 = f_rc1.normal();
	gmds::Cell::Data data_rc1 = m_fl->find(f_rc1.center());
	Node n_closest_f_rc1 = m_meshT->get<Node>(data_rc1.id);
	v_test = m_VectorField->value(n_closest_f_rc1.id()).normalize() ;
	if (v_rc1.dot(v_test) <= 0)
	{
		v_rc1 = - v_rc1;
	}

	math::Vector3d v_c0c1 = f_c0c1.normal();
	gmds::Cell::Data data_c0c1 = m_fl->find(f_c0c1.center());
	Node n_closest_f_c0c1 = m_meshT->get<Node>(data_c0c1.id);
	v_test = m_VectorField->value(n_closest_f_c0c1.id()).normalize() ;
	if (v_c0c1.dot(v_test) <= 0)
	{
		v_c0c1 = - v_c0c1;
	}
	*/

	math::Vector3d v_rc0 = m_Front->outgoingNormal(m_mesh, m_meshT, m_fl, m_VectorField, f_rc0.id()).normalize() ;
	math::Vector3d v_rc1 = m_Front->outgoingNormal(m_mesh, m_meshT, m_fl, m_VectorField, f_rc1.id()).normalize() ;
	math::Vector3d v_c0c1 = m_Front->outgoingNormal(m_mesh, m_meshT, m_fl, m_VectorField, f_c0c1.id()).normalize() ;

	math::Point p2 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, (v_rc0+v_rc1).normalize()+v_c0c1);
	n2.setPoint(p2);

	// Old way to compute new nodes position
	//math::Point p1 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, vr);
	//math::Point p3 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, (vc0+vc1).normalize());
	//math::Point p6 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, (vc0+vc1+v_rc0).normalize());
	//math::Point p10 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, (vc0+vc1+v_rc1).normalize());
	//math::Point p5 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, (vr+vc0+vc1+v_rc0).normalize());
	//math::Point p4 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, (vr+v_rc0).normalize());
	//math::Point p8 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, (vr+v_rc1).normalize());

	math::Point p1 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, v_c0c1);
	math::Point p3 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, v_rc0+v_rc1);
	math::Point p7 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, v_rc0);
	math::Point p4 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, v_c0c1+v_rc0);
	math::Point p6 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, v_rc0+(v_rc0+v_rc1).normalize());
	math::Point p5 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, v_c0c1+(v_rc0+v_rc1).normalize()+v_rc0);

	math::Point p11 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, v_rc1);
	math::Point p8 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, v_c0c1+v_rc1);
	math::Point p10 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, v_rc1+(v_rc0+v_rc1).normalize());
	math::Point p9 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, v_c0c1+(v_rc0+v_rc1).normalize()+v_rc1);

	n1 = m_mesh->newNode(p1);
	n3 = m_mesh->newNode(p3);
	n7 = m_mesh->newNode(p7);
	n4 = m_mesh->newNode(p4);
	n6 = m_mesh->newNode(p6);
	n5 = m_mesh->newNode(p5);
	n11 = m_mesh->newNode(p11);
	n8 = m_mesh->newNode(p8);
	n10 = m_mesh->newNode(p10);
	n9 = m_mesh->newNode(p9);
	// <----

	// Create the two new hexas
	TCellID r_id = math::Utils::CreateHexaNConnectivities(m_mesh, n, n1, n2, n3, n7, n4, n5, n6);
	m_hex.push_back(m_mesh->get<Region>(r_id));
	r_id = math::Utils::CreateHexaNConnectivities(m_mesh, n, n1, n2, n3, n11,n8, n9, n10);
	m_hex.push_back(m_mesh->get<Region>(r_id));

	// Update the layer ID for the 10 new nodes on the next layer ---->
	Variable<int>* var_node_couche_id = m_mesh->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	var_node_couche_id->set(n1.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n2.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n3.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n4.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n5.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n6.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n7.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n8.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n9.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n10.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n11.id(), m_Front->getFrontID()+1);
	// <----

	// Update the layer ID for the 6 new faces on the next layer ---->
	Variable<int>* var_face_couche_id = m_mesh->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");
	TCellID f_1_id = math::Utils::CommonFace3Nodes(m_mesh, n4.id(), n5.id(), n6.id());
	TCellID f_2_id = math::Utils::CommonFace3Nodes(m_mesh, n1.id(),  n2.id(), n5.id());
	TCellID f_3_id = math::Utils::CommonFace3Nodes(m_mesh, n3.id(),  n2.id(), n5.id());
	TCellID f_4_id = math::Utils::CommonFace3Nodes(m_mesh, n3.id(),  n2.id(), n9.id());
	TCellID f_5_id = math::Utils::CommonFace3Nodes(m_mesh, n1.id(),  n2.id(), n9.id());
	TCellID f_6_id = math::Utils::CommonFace3Nodes(m_mesh, n8.id(),  n9.id(), n10.id());
	var_face_couche_id->set(f_1_id, m_Front->getFrontID()+1);
	var_face_couche_id->set(f_2_id, m_Front->getFrontID()+1);
	var_face_couche_id->set(f_3_id, m_Front->getFrontID()+1);
	var_face_couche_id->set(f_4_id, m_Front->getFrontID()+1);
	var_face_couche_id->set(f_5_id, m_Front->getFrontID()+1);
	var_face_couche_id->set(f_6_id, m_Front->getFrontID()+1);
	// <----

	// Update the edge structure info <----
	std::pair<TCellID, TCellID> pair_rc0(f_rc0.id(), n.id());
	std::pair<TCellID, TCellID> pair_rc1(f_rc1.id(), n.id());
	std::pair<TCellID, TCellID> pair_c0c1(f_c0c1.id(), n.id());
	m_StructManager->setReversalFacesCreated(er.id(), m_n_id);
	m_StructManager->setReversalNextDiagNode(er.id(), f_rc0.id(), m_n_id, n6.id());
	m_StructManager->setReversalNextAdjNode(er.id(), f_rc0.id(), m_n_id, n7.id());
	m_StructManager->setReversalNextDiagNode(er.id(), f_rc1.id(), m_n_id, n10.id());
	m_StructManager->setReversalNextAdjNode(er.id(), f_rc1.id(), m_n_id, n11.id());
	m_StructManager->setReversalNextMediumNode(er.id(), m_n_id, n3.id());

	m_StructManager->setCornerFaceCreated(ec[0].id(), m_n_id);
	m_StructManager->setNextDiagNode(ec[0].id(), m_n_id, n4.id());
	m_StructManager->setCornerNextAdjNode(ec[0].id(), f_rc0.id(), m_n_id, n7.id());
	m_StructManager->setCornerNextAdjNode(ec[0].id(), f_c0c1.id(), m_n_id, n1.id());

	m_StructManager->setCornerFaceCreated(ec[1].id(), m_n_id);
	m_StructManager->setNextDiagNode(ec[1].id(), m_n_id, n8.id());
	m_StructManager->setCornerNextAdjNode(ec[1].id(), f_rc1.id(), m_n_id, n11.id());
	m_StructManager->setCornerNextAdjNode(ec[1].id(), f_c0c1.id(), m_n_id, n1.id());
	// ---->

	// Update the faces info ---->
	m_StructManager->setFaceNextNode(f_rc0.id(), m_n_id, n7.id());
	m_StructManager->setFaceNextNode(f_rc1.id(), m_n_id, n11.id());
	m_StructManager->setFaceNextNode(f_c0c1.id(), m_n_id, n1.id());
	// <----
}