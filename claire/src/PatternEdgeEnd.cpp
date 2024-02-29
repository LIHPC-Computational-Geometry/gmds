//
// Created by rochec on 28/04/23.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/PatternEdgeEnd.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/claire/NodeNeighbourhoodOnFront_3D.h>
#include <gmds/claire/Utils.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

PatternEdgeEnd::PatternEdgeEnd(Mesh *AMesh, Front_3D *AFront, TCellID Ae_id, LayerStructureManager_3D *AStructManager,
                                     Mesh *AMeshT, FastLocalize *Afl,
                                     double dc, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField) :
  AbstractPatternEdge(AMesh, AFront, Ae_id, AStructManager, AMeshT, Afl, dc, A_DistanceField, A_VectorField)
{

}
/*------------------------------------------------------------------------*/
void
PatternEdgeEnd::computeNewHex()
{

	Edge e = m_mesh->get<Edge>(m_e_id);
	std::vector<Node> e_nodes = e.get<Node>();

	std::vector<TCellID> e_front_faces = m_Front->edgeFacesOnFront(m_mesh, m_e_id) ;

	Variable<int>* var_node_couche_id = m_mesh->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* var_front_edges_classification = m_mesh->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");

	TCellID n0_id = e_nodes[0].id();
	TCellID n1_id = e_nodes[1].id();

	NodeNeighbourhoodOnFront_3D n0_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n0_id);
	n0_neighbourhood.execute();
	NodeNeighbourhoodOnFront_3D n1_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n1_id);
	n1_neighbourhood.execute();


	// First side
	TCellID e_00_id = n0_neighbourhood.nextEdgeOfFace(e_front_faces[0], m_e_id);
	Edge e_00 = m_mesh->get<Edge>(e_00_id);
	Node n_00 = e_00.getOppositeNode(m_mesh->get<Node>(n0_id));

	TCellID e_01_id = n0_neighbourhood.nextEdgeOfFace(e_front_faces[1], m_e_id);
	Edge e_01 = m_mesh->get<Edge>(e_01_id);
	Node n_01 = e_01.getOppositeNode(m_mesh->get<Node>(n0_id));

	Node n_02;
	bool side_n0_needUpdates(false);
	if (!m_StructManager->isEndFaceCreated(m_e_id, n0_id))
	{
		n_02 = m_mesh->get<Node>(m_StructManager->getFaceIdealNextNode(e_front_faces[0], n0_id));
		m_StructManager->setEndFaceCreated(m_e_id, n0_id);
		m_StructManager->setNextDiagNode(m_e_id, n0_id, n_02.id());
		side_n0_needUpdates = true;
	}
	else
	{
		n_02 = m_mesh->get<Node>(m_StructManager->getNextDiagNode(m_e_id, n0_id));
	}

	// Second side
	TCellID e_10_id = n1_neighbourhood.nextEdgeOfFace(e_front_faces[0], m_e_id);
	Edge e_10 = m_mesh->get<Edge>(e_10_id);
	Node n_10 = e_10.getOppositeNode(m_mesh->get<Node>(n1_id));

	TCellID e_11_id = n1_neighbourhood.nextEdgeOfFace(e_front_faces[1], m_e_id);
	Edge e_11 = m_mesh->get<Edge>(e_11_id);
	Node n_11 = e_11.getOppositeNode(m_mesh->get<Node>(n1_id));

	Node n_12;
	bool side_n1_needUpdates(false);
	if (!m_StructManager->isEndFaceCreated(m_e_id, n1_id))
	{
		n_12 = m_mesh->get<Node>(m_StructManager->getFaceIdealNextNode(e_front_faces[0], n1_id));
		m_StructManager->setEndFaceCreated(m_e_id, n1_id);
		m_StructManager->setNextDiagNode(m_e_id, n1_id, n_12.id());
		side_n1_needUpdates = true;
	}
	else
	{
		n_12 = m_mesh->get<Node>(m_StructManager->getNextDiagNode(m_e_id, n1_id));
	}

	// Create the hexa
	m_hex.push_back( m_mesh->get<Region>(math::Utils::CreateHexaNConnectivities(m_mesh, e_nodes[0], n_00, n_02, n_01, e_nodes[1], n_10, n_12, n_11)));

	//---------------------------//
	// Update the two new nodes  //
	// layer id                  //
	//---------------------------//
	var_node_couche_id->set(n_02.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n_12.id(), m_Front->getFrontID()+1);

	// Update the FaceInfo structures
	if (side_n0_needUpdates)
	{
		NodeNeighbourhoodOnFront_3D n00_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n_00.id());
		n00_neighbourhood.execute();
		NodeNeighbourhoodOnFront_3D n01_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n_01.id());
		n01_neighbourhood.execute();
		for (auto f_id:n00_neighbourhood.getOrderedFaces())
		{
			m_StructManager->setFaceNextNode(f_id, n_00.id(), n_02.id());
		}
		for (auto f_id:n01_neighbourhood.getOrderedFaces())
		{
			m_StructManager->setFaceNextNode(f_id, n_01.id(), n_02.id());
		}
	}
	if (side_n1_needUpdates)
	{
		NodeNeighbourhoodOnFront_3D n10_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n_10.id());
		n10_neighbourhood.execute();
		NodeNeighbourhoodOnFront_3D n11_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n_11.id());
		n11_neighbourhood.execute();
		for (auto f_id:n10_neighbourhood.getOrderedFaces())
		{
			m_StructManager->setFaceNextNode(f_id, n_10.id(), n_12.id());
		}
		for (auto f_id:n11_neighbourhood.getOrderedFaces())
		{
			m_StructManager->setFaceNextNode(f_id, n_11.id(), n_12.id());
		}
	}

	/*
	Edge e_opp_0 = math::Utils::oppositeEdgeInFace(m_mesh, m_e_id, e_front_faces[0]);
	TCellID f_opp_0 = m_Front->adjacentFaceOnFront(m_mesh, e_front_faces[0], e_opp_0.id());
	m_StructManager->setFaceNextNode(f_opp_0, n_00.id(), n_02.id());
	m_StructManager->setFaceNextNode(f_opp_0, n_10.id(), n_12.id());

	Edge e_opp_1 = math::Utils::oppositeEdgeInFace(m_mesh, m_e_id, e_front_faces[1]);
	TCellID f_opp_1 = m_Front->adjacentFaceOnFront(m_mesh, e_front_faces[1], e_opp_1.id());
	m_StructManager->setFaceNextNode(f_opp_1, n_01.id(), n_02.id());
	m_StructManager->setFaceNextNode(f_opp_1, n_11.id(), n_12.id());
	 */

	// Mark the two faces of the front as treated
	m_StructManager->setFaceTreated(e_front_faces[0]);
	m_StructManager->setFaceTreated(e_front_faces[1]);

}
/*------------------------------------------------------------------------*/