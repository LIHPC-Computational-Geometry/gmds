//
// Created by rochec on 27/04/23.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/PatternNode3Corner.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/claire/NodeNeighbourhoodOnFront_3D.h>
#include <gmds/claire/Utils.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

PatternNode3Corner::PatternNode3Corner(Mesh *AMesh, Front_3D *AFront, TCellID An_id, LayerStructureManager_3D *AStructManager,
                                       Mesh *AMeshT, FastLocalize *Afl,
                                       double dc, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField) :
  AbstractPatternNode(AMesh, AFront, An_id, AStructManager, AMeshT, Afl, dc, A_DistanceField, A_VectorField)
{

}
/*------------------------------------------------------------------------*/
void
PatternNode3Corner::computeNewHex()
{

	//std::cout << "Template Node 3 Corner au noeud " << n_id << std::endl;

	NodeNeighbourhoodOnFront_3D n_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, m_n_id);
	n_neighbourhood.execute();

	std::vector<TCellID> n_ordered_edges = n_neighbourhood.getOrderedEdges();
	Variable<int>* edge_classification = m_mesh->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
	std::vector<TCellID> n_edges_corner;
	for (auto edge:n_ordered_edges)
	{
		if (edge_classification->value(edge) == 1)
		{
			n_edges_corner.push_back(edge);
		}
	}

	Edge ec_0 = m_mesh->get<Edge>(n_edges_corner[0]);
	Edge ec_1 = m_mesh->get<Edge>(n_edges_corner[1]);
	Edge ec_2 = m_mesh->get<Edge>(n_edges_corner[2]);

	Node n = m_mesh->get<Node>(m_n_id);
	Node n6 = m_mesh->get<Node>(m_StructManager->getFaceNextNode(n_neighbourhood.getOrderedFaces()[0],m_n_id));

	Node n_c0 = ec_0.getOppositeNode(n);
	Node n_c1 = ec_1.getOppositeNode(n);
	Node n_c2 = ec_2.getOppositeNode(n);

	math::Vector3d v1 = m_Front->outgoingNormal(m_mesh, m_meshT, m_fl, m_VectorField, n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec_1.id(), ec_2.id(), ec_0.id())) ;
	math::Vector3d v2 = m_Front->outgoingNormal(m_mesh, m_meshT, m_fl, m_VectorField, n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec_0.id(), ec_2.id(), ec_1.id())) ;
	math::Vector3d v3 = m_Front->outgoingNormal(m_mesh, m_meshT, m_fl, m_VectorField, n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec_0.id(), ec_1.id(), ec_2.id())) ;
	v1.normalize();
	v2.normalize();
	v3.normalize();
	math::Point p1 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc,  m_DistanceField, v1);
	math::Point p3 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc,  m_DistanceField, v2);
	math::Point p4 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc,  m_DistanceField, v3);
	math::Point p2 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc,  m_DistanceField, v1+v2);
	math::Point p5 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc,  m_DistanceField, v1+v3);
	math::Point p7 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc,  m_DistanceField, v2+v3);

	Node n1 = m_mesh->newNode(p1);
	Node n2 = m_mesh->newNode(p2);
	Node n3 = m_mesh->newNode(p3);
	Node n4 = m_mesh->newNode(p4);
	Node n5 = m_mesh->newNode(p5);
	Node n7 = m_mesh->newNode(p7);

	// Update the layer ids
	Variable<int>* var_node_couche_id = m_mesh->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	var_node_couche_id->set(n1.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n2.id(),  m_Front->getFrontID()+1);
	var_node_couche_id->set(n3.id(),  m_Front->getFrontID()+1);
	var_node_couche_id->set(n4.id(),  m_Front->getFrontID()+1);
	var_node_couche_id->set(n5.id(),  m_Front->getFrontID()+1);
	var_node_couche_id->set(n7.id(),  m_Front->getFrontID()+1);

	m_hex.push_back(m_mesh->get<Region>(math::Utils::CreateHexaNConnectivities(m_mesh, n, n1, n2, n3, n4, n5, n6, n7)));

	// Update the faces information of the front
	std::vector<TCellID> adj_faces = n_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec_0.id(), ec_1.id(), ec_2.id());
	for (auto f_adj_id:adj_faces)
	{
		m_StructManager->setFaceNextNode(f_adj_id, m_n_id, n4.id());
	}
	adj_faces = n_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec_1.id(), ec_2.id(), ec_0.id());
	for (auto f_adj_id:adj_faces)
	{
		m_StructManager->setFaceNextNode(f_adj_id, m_n_id, n1.id());
	}
	adj_faces = n_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec_0.id(), ec_2.id(), ec_1.id());
	for (auto f_adj_id:adj_faces)
	{
		m_StructManager->setFaceNextNode(f_adj_id, m_n_id, n3.id());
	}

	TCellID f_1 = math::Utils::CommonFace3Nodes(m_mesh, m_n_id, n_c0.id(), n_c1.id());
	TCellID f_2 = math::Utils::CommonFace3Nodes(m_mesh, m_n_id, n_c1.id(), n_c2.id());
	TCellID f_3 = math::Utils::CommonFace3Nodes(m_mesh, m_n_id, n_c2.id(), n_c0.id());

	//---------------------------//
	// Update the corner edges   //
	// information					  //
	//---------------------------//
	TCellID e1_id = math::Utils::CommonEdge(m_mesh, n.id(), n_c0.id());
	TCellID e2_id = math::Utils::CommonEdge(m_mesh, n.id(), n_c1.id());
	TCellID e3_id = math::Utils::CommonEdge(m_mesh, n.id(), n_c2.id());

	m_StructManager->setCornerFaceCreated(e1_id, n.id());
	m_StructManager->setCornerFaceCreated(e2_id, n.id());
	m_StructManager->setCornerFaceCreated(e3_id, n.id());

	/*
	std::pair<TCellID, TCellID> pair_1(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e1_id, e2_id, e3_id), n.id());
	std::pair<TCellID, TCellID> pair_2(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e1_id, e3_id, e2_id), n.id());
	std::pair<TCellID, TCellID> pair_3(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e2_id, e1_id, e3_id), n.id());
	std::pair<TCellID, TCellID> pair_4(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e2_id, e3_id, e1_id), n.id());
	std::pair<TCellID, TCellID> pair_5(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e3_id, e2_id, e1_id), n.id());
	std::pair<TCellID, TCellID> pair_6(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e3_id, e1_id, e2_id), n.id());
	 */


	m_StructManager->setCornerNextAdjNode(e1_id, n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e1_id, e2_id, e3_id), n.id(), n4.id());
	m_StructManager->setCornerNextAdjNode(e1_id, n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e1_id, e3_id, e2_id), n.id(), n3.id());

	m_StructManager->setCornerNextAdjNode(e2_id, n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e2_id, e1_id, e3_id), n.id(), n4.id());
	m_StructManager->setCornerNextAdjNode(e2_id, n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e2_id, e3_id, e1_id), n.id(), n1.id());

	m_StructManager->setCornerNextAdjNode(e3_id, n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e3_id, e2_id, e1_id), n.id(), n1.id());
	m_StructManager->setCornerNextAdjNode(e3_id, n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e3_id, e1_id, e2_id), n.id(), n3.id());

	m_StructManager->setNextDiagNode(e1_id, n.id(), n7.id());
	m_StructManager->setNextDiagNode(e2_id, n.id(), n5.id());
	m_StructManager->setNextDiagNode(e3_id, n.id(), n2.id());

	//---------------------------//
	// Update the face layer id  //
	//---------------------------//
	Variable<int>* var_face_couche_id = m_mesh->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");
	TCellID f_1_new_id = math::Utils::CommonFace3Nodes(m_mesh, n4.id(), n5.id(), n6.id());
	TCellID f_2_new_id = math::Utils::CommonFace3Nodes(m_mesh, n1.id(),  n2.id(), n5.id());
	TCellID f_3_new_id = math::Utils::CommonFace3Nodes(m_mesh, n2.id(),  n3.id(), n6.id());
	var_face_couche_id->set(f_1_new_id,  m_Front->getFrontID()+1);
	var_face_couche_id->set(f_2_new_id,  m_Front->getFrontID()+1);
	var_face_couche_id->set(f_3_new_id,  m_Front->getFrontID()+1);

}