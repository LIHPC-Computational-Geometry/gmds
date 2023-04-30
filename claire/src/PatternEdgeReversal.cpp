//
// Created by rochec on 28/04/23.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/PatternEdgeReversal.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/claire/NodeNeighbourhoodOnFront_3D.h>
#include <gmds/claire/Utils.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

PatternEdgeReversal::PatternEdgeReversal(Mesh *AMesh, Front_3D *AFront, TCellID Ae_id, LayerStructureManager_3D *AStructManager,
                               Mesh *AMeshT, FastLocalize *Afl,
                               double dc, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField) :
  AbstractPatternEdge(AMesh, AFront, Ae_id, AStructManager, AMeshT, Afl, dc, A_DistanceField, A_VectorField)
{

}
/*------------------------------------------------------------------------*/
void
PatternEdgeReversal::computeNewHex()
{
	Node n0 = ((m_mesh->get<Edge>(m_e_id)).get<Node>())[0];
	Node n1 = ((m_mesh->get<Edge>(m_e_id)).get<Node>())[1];

	Node n2, n3, n4, n5, n6, n7, n8, n9, n10, n11;

	std::vector<TCellID> e_faces = m_Front->edgeFacesOnFront(m_mesh, m_e_id);

	// First side
	if (m_StructManager->areReversalFacesCreated(m_e_id, n0.id()))	// The faces are already created
	{
		n4 = m_mesh->get<Node>(m_StructManager->getReversalNextAdjNode(m_e_id, e_faces[0], n0.id()));
		n7 = m_mesh->get<Node>(m_StructManager->getReversalNextDiagNode(m_e_id, e_faces[0], n0.id()));
		n8 = m_mesh->get<Node>(m_StructManager->getReversalNextAdjNode(m_e_id, e_faces[1], n0.id()));
		n11 = m_mesh->get<Node>(m_StructManager->getReversalNextDiagNode(m_e_id, e_faces[1], n0.id()));
		n3 = m_mesh->get<Node>(m_StructManager->getReversalNextMediumNode(m_e_id, n0.id()));
	}
	else	// Create the faces
	{
		n3 = m_mesh->get<Node>(m_StructManager->getFaceIdealNextNode(e_faces[0], n0.id()));
		math::Vector3d v = (n3.point()-n0.point()).normalize() ;

		math::Point f0_barycentre = m_mesh->get<Face>(e_faces[0]).center() ;
		math::Vector3d f0_normale = m_mesh->get<Face>(e_faces[0]).normal();
		gmds::Cell::Data data = m_fl->find(f0_barycentre);
		Node n_closest = m_meshT->get<Node>(data.id);
		math::Vector3d v1 = m_VectorField->value(n_closest.id()).normalize() ;
		if (f0_normale.dot(v1) <= 0)
		{
			f0_normale = - f0_normale;
		}

		math::Point p4 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n0.point(), m_dc, m_DistanceField, f0_normale);
		math::Point p7 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n0.point(), m_dc, m_DistanceField, (f0_normale+v).normalize());

		n4 = m_mesh->newNode(p4);
		n7 = m_mesh->newNode(p7);

		math::Point f1_barycentre = m_mesh->get<Face>(e_faces[1]).center() ;
		math::Vector3d f1_normale = m_mesh->get<Face>(e_faces[1]).normal();
		gmds::Cell::Data data_f1 = m_fl->find(f1_barycentre);
		Node n_closest_f1 = m_meshT->get<Node>(data_f1.id);
		math::Vector3d v_f1 = m_VectorField->value(n_closest_f1.id()).normalize() ;
		if (f1_normale.dot(v_f1) <= 0)
		{
			f1_normale = - f1_normale;
		}

		math::Point p8 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n0.point(), m_dc, m_DistanceField, f1_normale);
		math::Point p11 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n0.point(), m_dc, m_DistanceField, (f1_normale+v).normalize());

		n8 = m_mesh->newNode(p8);
		n11 = m_mesh->newNode(p11);

		//=====================//
		// Update the REVERSAL //
		// edge info for the	  //
		// next edge			  //
		//=====================//
		Variable<int>* edge_classification = m_mesh->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
		NodeNeighbourhoodOnFront_3D n0_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n0.id());
		n0_neighbourhood.execute();
		Edge next_edge;
		for (auto e_loc_id:n0_neighbourhood.getOrderedEdges())
		{
			if ((e_loc_id != m_e_id)
			    && edge_classification->value(e_loc_id)==3 )
			{
				next_edge = m_mesh->get<Edge>(e_loc_id);
			}
		}

		Edge e_f0 = m_mesh->get<Edge>(n0_neighbourhood.nextEdgeOfFace(e_faces[0], m_e_id));
		Edge e_f1 = m_mesh->get<Edge>(n0_neighbourhood.nextEdgeOfFace(e_faces[1], m_e_id));

		Face f0_next = m_mesh->get<Face>(n0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(next_edge.id(), e_f0.id(), e_f1.id())) ;
		Face f1_next = m_mesh->get<Face>(n0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(next_edge.id(), e_f1.id(), e_f0.id())) ;

		m_StructManager->setReversalFacesCreated(next_edge.id(), n0.id());
		m_StructManager->setReversalNextMediumNode(next_edge.id(), n0.id(), n3.id());
		m_StructManager->setReversalNextAdjNode(next_edge.id(), f0_next.id(), n0.id(), n4.id());
		m_StructManager->setReversalNextAdjNode(next_edge.id(), f1_next.id(), n0.id(), n8.id());
		m_StructManager->setReversalNextDiagNode(next_edge.id(), f0_next.id(), n0.id(), n7.id());
		m_StructManager->setReversalNextDiagNode(next_edge.id(), f1_next.id(), n0.id(), n11.id());
		//==============//
		// END UPDATES  //
		//==============//

	}

	// Second side
	if (m_StructManager->areReversalFacesCreated(m_e_id, n1.id()))	// The faces are already created
	{
		n5 = m_mesh->get<Node>(m_StructManager->getReversalNextAdjNode(m_e_id, e_faces[0], n1.id()));
		n6 = m_mesh->get<Node>(m_StructManager->getReversalNextDiagNode(m_e_id, e_faces[0], n1.id()));
		n9 = m_mesh->get<Node>(m_StructManager->getReversalNextAdjNode(m_e_id, e_faces[1], n1.id()));
		n10 = m_mesh->get<Node>(m_StructManager->getReversalNextDiagNode(m_e_id, e_faces[1], n1.id()));
		n2 = m_mesh->get<Node>(m_StructManager->getReversalNextMediumNode(m_e_id, n1.id()));
	}
	else	// Create the faces
	{
		n2 = m_mesh->get<Node>(m_StructManager->getFaceIdealNextNode(e_faces[0], n1.id())) ;
		math::Vector3d v = (n2.point()-n1.point()).normalize() ;

		math::Point f0_barycentre = m_mesh->get<Face>(e_faces[0]).center() ;
		math::Vector3d f0_normale = m_mesh->get<Face>(e_faces[0]).normal();
		gmds::Cell::Data data = m_fl->find(f0_barycentre);
		Node n_closest = m_meshT->get<Node>(data.id);
		math::Vector3d v1 = m_VectorField->value(n_closest.id()).normalize() ;
		if (f0_normale.dot(v1) <= 0)
		{
			f0_normale = - f0_normale;
		}

		math::Point p5 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n1.point(), m_dc, m_DistanceField, f0_normale);
		math::Point p6 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n1.point(), m_dc, m_DistanceField, (f0_normale+v).normalize());

		n5 = m_mesh->newNode(p5);
		n6 = m_mesh->newNode(p6);

		math::Point f1_barycentre = m_mesh->get<Face>(e_faces[1]).center() ;
		math::Vector3d f1_normale = m_mesh->get<Face>(e_faces[1]).normal();
		gmds::Cell::Data data_f1 = m_fl->find(f1_barycentre);
		Node n_closest_f1 = m_meshT->get<Node>(data_f1.id);
		math::Vector3d v_f1 = m_VectorField->value(n_closest_f1.id()).normalize() ;
		if (f1_normale.dot(v_f1) <= 0)
		{
			f1_normale = - f1_normale;
		}

		math::Point p9 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n1.point(), m_dc, m_DistanceField, f1_normale);
		math::Point p10 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, n1.point(), m_dc, m_DistanceField, (f1_normale+v).normalize());

		n9 = m_mesh->newNode(p9);
		n10 = m_mesh->newNode(p10);

		//=====================//
		// Update the REVERSAL //
		// edge info for the	  //
		// next edge			  //
		//=====================//
		Variable<int>* edge_classification = m_mesh->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
		NodeNeighbourhoodOnFront_3D n1_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n1.id());
		n1_neighbourhood.execute();
		Edge next_edge;
		for (auto e_loc_id:n1_neighbourhood.getOrderedEdges())
		{
			if ((e_loc_id != m_e_id)
			    && edge_classification->value(e_loc_id)==3 )
			{
				next_edge = m_mesh->get<Edge>(e_loc_id);
			}
		}

		Edge e_f0 = m_mesh->get<Edge>(n1_neighbourhood.nextEdgeOfFace(e_faces[0], m_e_id));
		Edge e_f1 = m_mesh->get<Edge>(n1_neighbourhood.nextEdgeOfFace(e_faces[1], m_e_id));

		Face f0_next = m_mesh->get<Face>(n1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(next_edge.id(), e_f0.id(), e_f1.id())) ;
		Face f1_next = m_mesh->get<Face>(n1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(next_edge.id(), e_f1.id(), e_f0.id())) ;

		m_StructManager->setReversalFacesCreated(next_edge.id(), n1.id());
		m_StructManager->setReversalNextMediumNode(next_edge.id(), n1.id(), n2.id());
		m_StructManager->setReversalNextAdjNode(next_edge.id(), f0_next.id(), n1.id(), n5.id());
		m_StructManager->setReversalNextAdjNode(next_edge.id(), f1_next.id(), n1.id(), n9.id());
		m_StructManager->setReversalNextDiagNode(next_edge.id(), f0_next.id(), n1.id(), n6.id());
		m_StructManager->setReversalNextDiagNode(next_edge.id(), f1_next.id(), n1.id(), n10.id());
		//==============//
		// END UPDATES //
		//==============//
	}

	TCellID r_id = math::Utils::CreateHexaNConnectivities(m_mesh, n0, n1, n2, n3, n4, n5, n6, n7);
	m_hex.push_back(m_mesh->get<Region>(r_id));
	r_id = math::Utils::CreateHexaNConnectivities(m_mesh, n0, n1, n2, n3, n8, n9, n10, n11);
	m_hex.push_back(m_mesh->get<Region>(r_id));

	// Update the layer id of the new nodes ---->
	Variable<int>* var_node_couche_id = m_mesh->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	var_node_couche_id->set(n4.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n7.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n3.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n11.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n8.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n5.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n6.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n2.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n10.id(), m_Front->getFrontID()+1);
	var_node_couche_id->set(n9.id(), m_Front->getFrontID()+1);
	// <----

	// Update the layer ID for the 6 new faces on the next layer ---->
	Variable<int>* var_face_couche_id = m_mesh->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");
	TCellID f_1_id = math::Utils::CommonFace3Nodes(m_mesh, n4.id(), n5.id(), n6.id());
	TCellID f_2_id = math::Utils::CommonFace3Nodes(m_mesh, n2.id(), n3.id(), n7.id());
	TCellID f_3_id = math::Utils::CommonFace3Nodes(m_mesh, n2.id(), n3.id(), n11.id());
	TCellID f_4_id = math::Utils::CommonFace3Nodes(m_mesh, n9.id(), n10.id(), n11.id());
	var_face_couche_id->set(f_1_id, m_Front->getFrontID()+1);
	var_face_couche_id->set(f_2_id, m_Front->getFrontID()+1);
	var_face_couche_id->set(f_3_id, m_Front->getFrontID()+1);
	var_face_couche_id->set(f_4_id, m_Front->getFrontID()+1);
	// <----

	// Update the Face Info structures ---->
	m_StructManager->setFaceNextNode(e_faces[0], n0.id(), n4.id());
	m_StructManager->setFaceNextNode(e_faces[0], n1.id(), n5.id());
	m_StructManager->setFaceNextNode(e_faces[1], n0.id(), n8.id());
	m_StructManager->setFaceNextNode(e_faces[1], n1.id(), n9.id());
	// <----
}