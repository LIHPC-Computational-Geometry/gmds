//
// Created by rochec on 28/04/23.
//
/*------------------------------------------------------------------------*/
#include <gmds/aero/PatternNode3End.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/aero/NodeNeighbourhoodOnFront_3D.h>
#include <gmds/aero/Utils.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

PatternNode3End::PatternNode3End(Mesh *AMesh, Front_3D *AFront, TCellID An_id, LayerStructureManager_3D *AStructManager,
                                       Mesh *AMeshT, FastLocalize *Afl,
                                       double dc, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField) :
  AbstractPatternNode(AMesh, AFront, An_id, AStructManager, AMeshT, Afl, dc, A_DistanceField, A_VectorField)
{

}
/*------------------------------------------------------------------------*/
void
PatternNode3End::computeNewHex()
{
	Node n = m_mesh->get<Node>(m_n_id);

	NodeNeighbourhoodOnFront_3D n_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, m_n_id);
	n_neighbourhood.execute();
	std::vector<TCellID> n_ordered_edges = n_neighbourhood.getOrderedEdges();

	Edge e_0 = m_mesh->get<Edge>(n_ordered_edges[0]);
	Edge e_1 = m_mesh->get<Edge>(n_ordered_edges[1]);
	Edge e_2 = m_mesh->get<Edge>(n_ordered_edges[2]);

	Node n_e0_opp = e_0.getOppositeNode(n);
	Node n_e1_opp = e_1.getOppositeNode(n);
	Node n_e2_opp = e_2.getOppositeNode(n);

	TCellID f0_id = math::Utils::CommonFace3Nodes(m_mesh, m_n_id, n_e0_opp.id(), n_e1_opp.id());
	TCellID f1_id = math::Utils::CommonFace3Nodes(m_mesh, m_n_id, n_e1_opp.id(), n_e2_opp.id());
	TCellID f2_id = math::Utils::CommonFace3Nodes(m_mesh, m_n_id, n_e2_opp.id(), n_e0_opp.id());

	Face f0 = m_mesh->get<Face>(f0_id);
	Face f1 = m_mesh->get<Face>(f1_id);
	Face f2 = m_mesh->get<Face>(f2_id);

	std::vector<Node> f0_nodes = f0.get<Node>();
	std::vector<Node> f1_nodes = f1.get<Node>();
	std::vector<Node> f2_nodes = f2.get<Node>();

	Node n_f0_diag;
	Node n_f1_diag;
	Node n_f2_diag;
	for (auto const& n_loc:f0_nodes)
	{
		if (n_loc.id() != m_n_id && n_loc.id() != n_e0_opp.id() && n_loc.id() != n_e1_opp.id())
		{
			n_f0_diag = n_loc;
		}
	}
	for (auto const& n_loc:f1_nodes)
	{
		if (n_loc.id() != m_n_id && n_loc.id() != n_e1_opp.id() && n_loc.id() != n_e2_opp.id())
		{
			n_f1_diag = n_loc;
		}
	}
	for (auto const& n_loc:f2_nodes)
	{
		if (n_loc.id() != m_n_id && n_loc.id() != n_e2_opp.id() && n_loc.id() != n_e0_opp.id())
		{
			n_f2_diag = n_loc;
		}
	}

	Node n_new = m_mesh->get<Node>(m_StructManager->getFaceIdealNextNode(f0_id, m_n_id));

	// Create the new hexa
	m_hex.push_back(m_mesh->get<Region>( math::Utils::CreateHexaNConnectivities(m_mesh, n, n_e0_opp, n_f0_diag, n_e1_opp, n_e2_opp, n_f2_diag, n_new, n_f1_diag)));

	// Update the faces and edges that can't do actions anymore
	m_StructManager->setEdgeTreated(e_0.id());
	m_StructManager->setEdgeTreated(e_1.id());
	m_StructManager->setEdgeTreated(e_2.id());
	m_StructManager->setFaceTreated(f0.id());
	m_StructManager->setFaceTreated(f1.id());
	m_StructManager->setFaceTreated(f2.id());

	// Update the edges infos
	Variable<int>* var_front_edges_classification = m_mesh->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
	NodeNeighbourhoodOnFront_3D n_e0_opp_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n_e0_opp.id()) ;
	n_e0_opp_neighbourhood.execute();
	NodeNeighbourhoodOnFront_3D n_e1_opp_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n_e1_opp.id()) ;
	n_e1_opp_neighbourhood.execute();
	NodeNeighbourhoodOnFront_3D n_e2_opp_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n_e2_opp.id()) ;
	n_e2_opp_neighbourhood.execute();
	int compteur(0);
	for (auto e_loc_id:n_e0_opp_neighbourhood.getOrderedEdges())
	{
		if (var_front_edges_classification->value(e_loc_id) == 2 && e_loc_id != e_0.id())
		{
			m_StructManager->setNextDiagNode(e_loc_id, n_e0_opp.id(), n_new.id());
			m_StructManager->setEndFaceCreated(e_loc_id, n_e0_opp.id());
			compteur++;
		}
	}
	if (compteur != 1)
	{
		std::cout << "Attention." << std::endl;
	}

	compteur = 0;
	for (auto e_loc_id:n_e1_opp_neighbourhood.getOrderedEdges())
	{
		if (var_front_edges_classification->value(e_loc_id) == 2 && e_loc_id != e_1.id())
		{
			m_StructManager->setNextDiagNode(e_loc_id, n_e1_opp.id(), n_new.id());
			m_StructManager->setEndFaceCreated(e_loc_id, n_e1_opp.id());
			compteur++;
		}
	}
	if (compteur != 1)
	{
		std::cout << "Attention." << std::endl;
	}

	compteur = 0;
	for (auto e_loc_id:n_e2_opp_neighbourhood.getOrderedEdges())
	{
		if (var_front_edges_classification->value(e_loc_id) == 2 && e_loc_id != e_2.id())
		{
			m_StructManager->setNextDiagNode(e_loc_id, n_e2_opp.id(), n_new.id());
			m_StructManager->setEndFaceCreated(e_loc_id, n_e2_opp.id());
			compteur++;
		}
	}
	if (compteur != 1)
	{
		std::cout << "Attention." << std::endl;
	}


	// Update the FaceInfo structures
	NodeNeighbourhoodOnFront_3D node_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n_f0_diag.id());
	node_neighbourhood.execute();
	for (auto face_id:node_neighbourhood.getOrderedFaces())
	{
		m_StructManager->setFaceNextNode(face_id, n_f0_diag.id(), n_new.id());
	}
	node_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n_f1_diag.id());
	node_neighbourhood.execute();
	for (auto face_id:node_neighbourhood.getOrderedFaces())
	{
		m_StructManager->setFaceNextNode(face_id, n_f1_diag.id(), n_new.id());
	}
	node_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n_f2_diag.id());
	node_neighbourhood.execute();
	for (auto face_id:node_neighbourhood.getOrderedFaces())
	{
		m_StructManager->setFaceNextNode(face_id, n_f2_diag.id(), n_new.id());
	}
}