//
// Created by rochec on 28/04/23.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/PatternNode2End1Reversal.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/claire/NodeNeighbourhoodOnFront_3D.h>
#include <gmds/claire/Utils.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

PatternNode2End1Reversal::PatternNode2End1Reversal(Mesh *AMesh, Front_3D *AFront, TCellID An_id, LayerStructureManager_3D *AStructManager,
                                       Mesh *AMeshT, FastLocalize *Afl,
                                       double dc, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField) :
  AbstractPatternNode(AMesh, AFront, An_id, AStructManager, AMeshT, Afl, dc, A_DistanceField, A_VectorField)
{

}
/*------------------------------------------------------------------------*/
void
PatternNode2End1Reversal::computeNewHex()
{
	Node n = m_mesh->get<Node>(m_n_id);
	NodeNeighbourhoodOnFront_3D n_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, m_n_id);
	n_neighbourhood.execute();
	std::vector<TCellID> n_ordered_edges = n_neighbourhood.getOrderedEdges();
	std::vector<TCellID> n_ordered_faces = n_neighbourhood.getOrderedFaces();

	Variable<int>* var_front_edges_classification = m_mesh->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
	std::vector<Edge> ee;
	Edge er;
	for (auto e_loc_id:n_ordered_edges)
	{
		if (var_front_edges_classification->value(e_loc_id)==2)
		{
			ee.push_back(m_mesh->get<Edge>(e_loc_id));
		}
		else if (var_front_edges_classification->value(e_loc_id)==3)
		{
			er = m_mesh->get<Edge>(e_loc_id);
		}
	}

	// Get the nodes for the 2 new hexas
	//Node n2 = m_mesh->get<Node>(m_FaceInfo[n_ordered_faces[0]].next_ideal_nodes[m_n_id]);
	Node nr = er.getOppositeNode(n);

	Face f_re0 = m_mesh->get<Face>(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(er.id(), ee[0].id(), ee[1].id()));
	Face f_re1 = m_mesh->get<Face>(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(er.id(), ee[1].id(), ee[0].id()));

	Face f_0 = m_mesh->get<Face>( n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ee[0].id(), ee[1].id(), er.id()));
	Face f_1 = m_mesh->get<Face>( n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ee[1].id(), ee[0].id(), er.id()));

	Edge e_0 = m_mesh->get<Edge>( n_neighbourhood.nextEdgeOfFace(f_0.id(), ee[0].id()));
	Edge e_1 = m_mesh->get<Edge>( n_neighbourhood.nextEdgeOfFace(f_1.id(), ee[1].id()));

	Face f_2 = m_mesh->get<Face>( n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e_0.id(), e_1.id(), ee[0].id()) );
	Face f_3 = m_mesh->get<Face>( n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e_1.id(), e_0.id(), ee[1].id()) );

	Edge e_2 = m_mesh->get<Edge>(n_neighbourhood.nextEdgeOfFace(f_2.id(), e_0.id()));

	Node n1 = e_2.getOppositeNode(n);
	Node n4 = e_0.getOppositeNode(n);
	Node n8 = e_1.getOppositeNode(n);

	Edge e_2_opp_f2 = math::Utils::oppositeEdgeInFace(m_mesh, e_2.id(), f_2.id()) ;
	Edge e_2_opp_f3 = math::Utils::oppositeEdgeInFace(m_mesh, e_2.id(), f_3.id()) ;

	Node n5 = e_2_opp_f2.getOppositeNode(n4);
	Node n9 = e_2_opp_f3.getOppositeNode(n8);

	Node n2 = m_mesh->get<Node>(m_StructManager->getFaceIdealNextNode(f_2.id(), n1.id()));
	Node n3 = m_mesh->get<Node>(m_StructManager->getFaceIdealNextNode(f_2.id(), n4.id()));
	Node n6 = m_mesh->get<Node>(m_StructManager->getFaceIdealNextNode(f_2.id(), n5.id()));
	Node n7 = m_mesh->get<Node>(m_StructManager->getFaceIdealNextNode(f_3.id(), n8.id()));
	Node n10 = m_mesh->get<Node>(m_StructManager->getFaceIdealNextNode(f_3.id(), n9.id()));
	// <----

	// Create the two new hexas
	TCellID r_id = math::Utils::CreateHexaNConnectivities(m_mesh, nr, n, n1, n2, n3, n4, n5, n6);
	m_hex.push_back(m_mesh->get<Region>(r_id));
	r_id = math::Utils::CreateHexaNConnectivities(m_mesh, nr, n, n1, n2, n7, n8, n9, n10);
	m_hex.push_back(m_mesh->get<Region>(r_id));

	// Mark the faces and the edge of the front where the hexa insertion is not possible anymore
	m_StructManager->setFaceTreated(f_2.id());
	m_StructManager->setFaceTreated(f_3.id());
	m_StructManager->setEdgeTreated(er.id());

	// Update the EdgeInfo for the two END edges ---->
	m_StructManager->setEndFaceCreated(ee[0].id(), m_n_id);
	m_StructManager->setEndFaceCreated(ee[1].id(), m_n_id);
	m_StructManager->setNextDiagNode(ee[0].id(), m_n_id, n3.id());
	m_StructManager->setNextDiagNode(ee[1].id(), m_n_id, n7.id());
	// <----

	// Update the EdgeInfo for the next REVERSAL edge ---->
	NodeNeighbourhoodOnFront_3D nr_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, nr.id());
	nr_neighbourhood.execute();
	Edge er_next;
	for (auto e_loc_id:nr_neighbourhood.getOrderedEdges())
	{
		if (e_loc_id != er.id()
		    && var_front_edges_classification->value(e_loc_id) == 3)
		{
			er_next = m_mesh->get<Edge>(e_loc_id);
		}
	}
	Edge e_opp_e0_fre0 = math::Utils::oppositeEdgeInFace(m_mesh, ee[0].id(), f_re0.id());
	Edge e_opp_e1_fre1 = math::Utils::oppositeEdgeInFace(m_mesh, ee[1].id(), f_re1.id());

	Face f_0_er_next = m_mesh->get<Face>( nr_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(er_next.id(), e_opp_e0_fre0.id(), e_opp_e1_fre1.id()) ) ;
	Face f_1_er_next = m_mesh->get<Face>( nr_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(er_next.id(), e_opp_e1_fre1.id(), e_opp_e0_fre0.id()) ) ;

	std::pair<TCellID, TCellID> pair_0(f_0_er_next.id(), nr.id());
	std::pair<TCellID, TCellID> pair_1(f_1_er_next.id(), nr.id());

	m_StructManager->setReversalFacesCreated(er_next.id(), nr.id());
	m_StructManager->setReversalNextAdjNode(er_next.id(), f_0_er_next.id(), nr.id(), n3.id());
	m_StructManager->setReversalNextAdjNode(er_next.id(), f_1_er_next.id(), nr.id(), n7.id());
	m_StructManager->setReversalNextDiagNode(er_next.id(), f_0_er_next.id(), nr.id(), n6.id());
	m_StructManager->setReversalNextDiagNode(er_next.id(), f_1_er_next.id(), nr.id(), n10.id());
	m_StructManager->setReversalNextMediumNode(er_next.id(), nr.id(), n2.id());
	// <----

	// Update the FaceInfo for the faces around the next reversal edge ---->
	for (auto f_loc:nr_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(er_next.id(), er.id(), e_opp_e1_fre1.id()))
	{
		m_StructManager->setFaceNextNode(f_loc, nr.id(), n3.id());
	}
	for (auto f_loc:nr_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(er_next.id(), er.id(), e_opp_e0_fre0.id()))
	{
		m_StructManager->setFaceNextNode(f_loc, nr.id(), n7.id());
	}
	// <----
}