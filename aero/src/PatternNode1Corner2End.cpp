//
// Created by rochec on 28/04/23.
//
/*------------------------------------------------------------------------*/
#include <gmds/aero/PatternNode1Corner2End.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/aero/NodeNeighbourhoodOnFront_3D.h>
#include <gmds/aero/Utils.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

PatternNode1Corner2End::PatternNode1Corner2End(Mesh *AMesh, Front_3D *AFront, TCellID An_id, LayerStructureManager_3D *AStructManager,
                                       Mesh *AMeshT, FastLocalize *Afl,
                                       double dc, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField) :
  AbstractPatternNode(AMesh, AFront, An_id, AStructManager, AMeshT, Afl, dc, A_DistanceField, A_VectorField)
{

}
/*------------------------------------------------------------------------*/
void
PatternNode1Corner2End::computeNewHex()
{
	Node n = m_mesh->get<Node>(m_n_id);

	NodeNeighbourhoodOnFront_3D n_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, m_n_id);
	n_neighbourhood.execute();

	Variable<int>* var_front_edges_classification = m_mesh->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");

	std::vector<Edge> ee;	// For the two edge END
	Edge ec;						// For the edge CORNER

	// Get the local feature edges
	for (auto edge_id:n_neighbourhood.getOrderedEdges())
	{
		if (var_front_edges_classification->value(edge_id)==2)
		{
			ee.push_back(m_mesh->get<Edge>(edge_id));
		}
		else if (var_front_edges_classification->value(edge_id)==1)
		{
			ec = m_mesh->get<Edge>(edge_id);
		}
	}

	Node n_e0 = ee[0].getOppositeNode(n);
	Node n_e1 = ee[1].getOppositeNode(n);

	std::vector<TCellID> ee_faces = n_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ee[0].id(), ee[1].id(), ec.id()) ;
	if (ee_faces.size() != 3)
	{
		std::cout << "Attention AeroExtrusion_3D: Template Node 1 CORNER, 2 END, ne respecte pas la condition des 3 faces entre les deux Edges de type END." << std::endl;
	}

	Face f1 = m_mesh->get<Face>(ee_faces[1]);

	Edge e1 = m_mesh->get<Edge>(n_neighbourhood.nextEdgeOfFace(ee_faces[0], ee[0].id()));
	Edge e3 = m_mesh->get<Edge>(n_neighbourhood.nextEdgeOfFace(ee_faces[2], ee[1].id()));
	Node n1 = e1.getOppositeNode(n);
	Node n3 = e3.getOppositeNode(n);

	Node n2;
	for (auto const& n_loc:f1.get<Node>())
	{
		if (n_loc.id() != n.id()
		    && n_loc.id() != n1.id()
		    && n_loc.id() != n3.id())
		{
			n2 = n_loc;
		}
	}

	Node nc = ec.getOppositeNode(n) ;

	//Node n4 = m_mesh->get<Node>(m_FaceInfo[f1.id()].next_ideal_nodes[n1.id()]);
	//Node n5 = m_mesh->get<Node>(m_FaceInfo[f1.id()].next_ideal_nodes[n2.id()]);
	//Node n6 = m_mesh->get<Node>(m_FaceInfo[f1.id()].next_ideal_nodes[n3.id()]);
	Node n4 = m_mesh->get<Node>(m_StructManager->getFaceIdealNextNode(f1.id(), n1.id()));
	Node n5 = m_mesh->get<Node>(m_StructManager->getFaceIdealNextNode(f1.id(), n2.id()));
	Node n6 = m_mesh->get<Node>(m_StructManager->getFaceIdealNextNode(f1.id(), n3.id()));

	// Create the new hexa
	m_hex.push_back(m_mesh->get<Region>( math::Utils::CreateHexaNConnectivities(m_mesh, n, n1, n2, n3, nc, n4, n5, n6)));

	// Update the faces and edges that can't do actions anymore
	m_StructManager->setEdgeTreated(ec.id());
	m_StructManager->setFaceTreated(f1.id());

	// Update the EdgeInfo structure for the END edges ---->
	m_StructManager->setEndFaceCreated(ee[0].id(), m_n_id);
	m_StructManager->setNextDiagNode(ee[0].id(), m_n_id, n4.id());
	m_StructManager->setEndFaceCreated(ee[1].id(), m_n_id);
	m_StructManager->setNextDiagNode(ee[1].id(), m_n_id, n6.id());
	//<----

	// Update the EdgeInfo structure for the CORNER edge ---->
	NodeNeighbourhoodOnFront_3D nc_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, nc.id());
	nc_neighbourhood.execute();

	std::vector<TCellID> nc_front_edges = nc_neighbourhood.getOrderedEdges();
	Edge ec_next;
	for (auto e_loc_id:nc_front_edges)
	{
		if (var_front_edges_classification->value(e_loc_id) == 1
		    && e_loc_id != ec.id())
		{
			ec_next = m_mesh->get<Edge>(e_loc_id) ;
		}
	}
	m_StructManager->setCornerFaceCreated(ec_next.id(), nc.id());
	m_StructManager->setNextDiagNode(ec_next.id(), nc.id(), n5.id());

	Face f_ce_0 = m_mesh->get<Face>( n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec.id(), ee[0].id(), ee[1].id()) );
	Face f_ce_1 = m_mesh->get<Face>( n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec.id(), ee[1].id(), ee[0].id()) );

	Edge e_e0_opp = math::Utils::oppositeEdgeInFace(m_mesh, ee[0].id(), f_ce_0.id()) ;
	Edge e_e1_opp = math::Utils::oppositeEdgeInFace(m_mesh, ee[1].id(), f_ce_1.id()) ;

	Face f_ec_next_0 = m_mesh->get<Face>(nc_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec_next.id(), e_e0_opp.id(), e_e1_opp.id()) );
	Face f_ec_next_1 = m_mesh->get<Face>(nc_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec_next.id(), e_e1_opp.id(), e_e0_opp.id()) );

	m_StructManager->setCornerNextAdjNode(ec_next.id(), f_ec_next_0.id(), nc.id(), n4.id());
	m_StructManager->setCornerNextAdjNode(ec_next.id(), f_ec_next_1.id(), nc.id(), n6.id());
	//<----

	// Update the FaceInfo structures ---->
	for (auto f_loc_id:nc_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec_next.id(), ec.id(), e_e1_opp.id()))
	{
		m_StructManager->setFaceNextNode(f_loc_id, nc.id(), n4.id());
	}
	for (auto f_loc_id:nc_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec_next.id(), ec.id(), e_e0_opp.id()))
	{
		m_StructManager->setFaceNextNode(f_loc_id, nc.id(), n6.id());
	}
	//<----
}