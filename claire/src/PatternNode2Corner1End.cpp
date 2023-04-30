//
// Created by rochec on 28/04/23.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/PatternNode2Corner1End.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/claire/NodeNeighbourhoodOnFront_3D.h>
#include <gmds/claire/Utils.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

PatternNode2Corner1End::PatternNode2Corner1End(Mesh *AMesh, Front_3D *AFront, TCellID An_id, LayerStructureManager_3D *AStructManager,
                                       Mesh *AMeshT, FastLocalize *Afl,
                                       double dc, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField) :
  AbstractPatternNode(AMesh, AFront, An_id, AStructManager, AMeshT, Afl, dc, A_DistanceField, A_VectorField)
{

}
/*------------------------------------------------------------------------*/
void
PatternNode2Corner1End::computeNewHex()
{
	Node n = m_mesh->get<Node>(m_n_id);

	NodeNeighbourhoodOnFront_3D n_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, m_n_id);
	n_neighbourhood.execute();

	Variable<int>* var_front_edges_classification = m_mesh->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");

	std::vector<Edge> ec;	// For the two edge CORNER
	Edge ee;						// For the edge END

	// Get the local feature edges
	for (auto edge_id:n_neighbourhood.getOrderedEdges())
	{
		if (var_front_edges_classification->value(edge_id)==1)
		{
			ec.push_back(m_mesh->get<Edge>(edge_id));
		}
		else if (var_front_edges_classification->value(edge_id)==2)
		{
			ee = m_mesh->get<Edge>(edge_id);
		}
	}

	Node n_c0 = ec[0].getOppositeNode(n);
	Node n_c1 = ec[1].getOppositeNode(n);

	math::Vector3d ve = ( ( n_c0.point() - n.point() ).normalize() + ( n_c1.point() - n.point() ).normalize() ).normalize() ;
	math::Point pe = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, ve);
	Node n_e = m_mesh->newNode(pe);

	math::Vector3d v1 = (n.point()- (ee.getOppositeNode(n)).point() ).normalize() ;
	math::Point p1 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, v1);
	Node n_1 = m_mesh->newNode(p1);

	math::Point p2 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n_c0.point(), m_dc, m_DistanceField, v1);
	Node n_2 = m_mesh->newNode(p2);

	math::Point p4 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n_c1.point(), m_dc, m_DistanceField, v1);
	Node n_4 = m_mesh->newNode(p4);

	math::Vector3d v3 = ( ( n_c0.point() - n.point() ).normalize() + ( n_c1.point() - n.point() ).normalize() + v1 ).normalize() ;
	math::Point p3 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n.point(), m_dc, m_DistanceField, v3);
	Node n_3 = m_mesh->newNode(p3);

	// Create the new hexa
	m_hex.push_back(m_mesh->get<Region>(math::Utils::CreateHexaNConnectivities(m_mesh, n, n_c0, n_e, n_c1, n_1, n_2, n_3, n_4)));

	// Update the layer ids
	Variable<int>* var_node_couche_id = m_mesh->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	var_node_couche_id->set(n_1.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n_2.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n_3.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n_4.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n_e.id(), m_Front->getFrontID()+1);
	Variable<int>* var_face_couche_id = m_mesh->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");
	var_face_couche_id->set(math::Utils::CommonFace3Nodes(m_mesh, n_1.id(), n_2.id(), n_3.id()), m_Front->getFrontID()+1);

	// Update the faces and edges that can't do actions anymore
	m_StructManager->setEdgeTreated(ec[0].id());
	m_StructManager->setEdgeTreated(ec[1].id());

	// Update the FaceInfo structures
	for (auto face_id:n_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[0].id(), ec[1].id(), ee.id()))
	{
		m_StructManager->setFaceNextNode(face_id, m_n_id, n_1.id());
	}

	m_StructManager->setFaceNextNode(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec[0].id(), ec[1].id(), ee.id()), n_c0.id(), n_2.id());
	m_StructManager->setFaceNextNode(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec[1].id(), ec[0].id(), ee.id()), n_c1.id(), n_4.id());

	TCellID f_e0_id = math::Utils::CommonFace3Nodes(m_mesh, n_c0.id(), m_n_id, ee.getOppositeNode(n).id());
	m_StructManager->setFaceNextNode(f_e0_id, n_c0.id(), n_e.id());
	TCellID f_e1_id = math::Utils::CommonFace3Nodes(m_mesh, n_c1.id(), m_n_id, ee.getOppositeNode(n).id());
	m_StructManager->setFaceNextNode(f_e1_id, n_c1.id(), n_e.id());


	// Update the EdgeInfo structure for the END edge
	m_StructManager->setEndFaceCreated(ee.id(), m_n_id);
	m_StructManager->setNextDiagNode(ee.id(), m_n_id, n_e.id());

	// Update the EdgeInfo structure for the two CORNER edges
	NodeNeighbourhoodOnFront_3D n_c0_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n_c0.id()) ;
	n_c0_neighbourhood.execute();
	for (auto edge_id:n_c0_neighbourhood.getOrderedEdges())
	{
		if (var_front_edges_classification->value(edge_id) == 1
		    && edge_id != ec[0].id())
		{
			m_StructManager->setCornerFaceCreated(edge_id, n_c0.id());
			m_StructManager->setNextDiagNode(edge_id, n_c0.id(), n_3.id());

			TCellID edge_f_e0 = n_c0_neighbourhood.nextEdgeOfFace(f_e0_id, ec[0].id()) ;
			for (auto faces_id:n_c0_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[0].id(), edge_id, edge_f_e0))
			{
				m_StructManager->setFaceNextNode(faces_id, n_c0.id(), n_2.id());
			}

			m_StructManager->setCornerNextAdjNode(edge_id, n_c0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[0].id(), edge_f_e0), n_c0.id(), n_2.id());
			m_StructManager->setCornerNextAdjNode(edge_id, n_c0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, edge_f_e0, ec[0].id()), n_c0.id(), n_e.id());

		}
	}

	NodeNeighbourhoodOnFront_3D n_c1_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n_c1.id()) ;
	n_c1_neighbourhood.execute();
	for (auto edge_id:n_c1_neighbourhood.getOrderedEdges())
	{
		if (var_front_edges_classification->value(edge_id) == 1
		    && edge_id != ec[1].id())
		{
			m_StructManager->setCornerFaceCreated(edge_id, n_c1.id());
			m_StructManager->setNextDiagNode(edge_id, n_c1.id(), n_3.id());

			TCellID edge_f_e1 = n_c1_neighbourhood.nextEdgeOfFace(f_e1_id, ec[1].id()) ;
			for (auto faces_id:n_c1_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[1].id(), edge_id, edge_f_e1))
			{
				m_StructManager->setFaceNextNode(faces_id, n_c1.id(), n_4.id());
			}

			m_StructManager->setCornerNextAdjNode(edge_id, n_c1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[1].id(), edge_f_e1), n_c1.id(), n_4.id());
			m_StructManager->setCornerNextAdjNode(edge_id, n_c1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, edge_f_e1, ec[1].id()), n_c1.id(), n_e.id());

		}
	}
}