//
// Created by rochec on 28/04/23.
//

/*------------------------------------------------------------------------*/
#include <gmds/aero/PatternEdgeCorner.h>
#include <gmds/ig/Mesh.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/aero/NodeNeighbourhoodOnFront_3D.h>
#include <gmds/aero/Utils.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

PatternEdgeCorner::PatternEdgeCorner(Mesh *AMesh, Front_3D *AFront, TCellID Ae_id, LayerStructureManager_3D *AStructManager,
                                       Mesh *AMeshT, FastLocalize *Afl,
                                       double dc, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField) :
  AbstractPatternEdge(AMesh, AFront, Ae_id, AStructManager, AMeshT, Afl, dc, A_DistanceField, A_VectorField)
{

}
/*------------------------------------------------------------------------*/
void
PatternEdgeCorner::computeNewHex()
{
	Edge e = m_mesh->get<Edge>(m_e_id);
	std::vector<Node> e_nodes = e.get<Node>();

	std::vector<TCellID> e_front_faces = m_Front->edgeFacesOnFront(m_mesh, m_e_id) ;


	Variable<int>* var_node_couche_id = m_mesh->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* var_front_edges_classification = m_mesh->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");

	// First side
	TCellID n0_id = e_nodes[0].id();
	TCellID n0_1_id(NullID);
	TCellID n0_2_id(NullID);
	TCellID n0_3_id(NullID);
	//if (m_EdgeInfo[e_id].CORNER_n_face_created[n0_id])
	if (m_StructManager->isCornerFaceCreated(m_e_id, n0_id))
	{
		// Then, a face has already been inserted on node e_nodes[0] with a node template
		n0_1_id = m_StructManager->getCornerNextAdjNode(m_e_id, e_front_faces[0], n0_id);
		n0_2_id = m_StructManager->getNextDiagNode(m_e_id, n0_id);
		n0_3_id = m_StructManager->getCornerNextAdjNode(m_e_id, e_front_faces[1], n0_id);
	}
	else
	{
		//============================//
		// We create the 3 new nodes	//
		//============================//

		/*
		math::Vector3d f0_normal = m_Front->outgoingNormal(m_mesh, m_meshT, m_fl, m_VectorField, e_front_faces[0]) ;
		f0_normal.normalize();
		f0_normal = e.length()*f0_normal ;
		 */
		math::Vector3d f0_normal = computeNormaltoFacesAroundNodeSideFace(n0_id, e_front_faces[0]);

		/*
		math::Vector3d f1_normal = m_Front->outgoingNormal(m_mesh, m_meshT, m_fl, m_VectorField, e_front_faces[1]) ;
		f1_normal.normalize();
		f1_normal = e.length()*f1_normal ;
		 */
		math::Vector3d f1_normal = computeNormaltoFacesAroundNodeSideFace(n0_id, e_front_faces[1]);

		math::Point p_n0 = e_nodes[0].point();

		n0_2_id = m_StructManager->getFaceIdealNextNode(e_front_faces[0], e_nodes[0].id());
		//Node n0_2_new = m_mesh->get<Node>(n0_2_id);
		//math::Point p0_2_new = n0_2_new.point() ;
		math::Point p0_2_new = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, e_nodes[0].point(), m_dc, m_DistanceField, f0_normal.normalize()+f1_normal.normalize());
		Node n0_2_new = m_mesh->newNode(p0_2_new);
		n0_2_id = n0_2_new.id();
		math::Point p0_1_new = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, e_nodes[0].point(), m_dc, m_DistanceField, f0_normal.normalize());
		math::Point p0_3_new = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, e_nodes[0].point(), m_dc, m_DistanceField, f1_normal.normalize());

		Node n0_1_new = m_mesh->newNode(p0_1_new);
		Node n0_3_new = m_mesh->newNode(p0_3_new);

		n0_1_id = n0_1_new.id();
		n0_2_id = n0_2_new.id();
		n0_3_id = n0_3_new.id();

		//---------------------------//
		// Update the layer ids of   //
		// the new nodes				  //
		//---------------------------//
		var_node_couche_id->set(n0_1_id, m_Front->getFrontID()+1);
		var_node_couche_id->set(n0_2_id, m_Front->getFrontID()+1);
		var_node_couche_id->set(n0_3_id, m_Front->getFrontID()+1);

		//---------------------------//
		// Update the edge corner	  //
		// information around the	  //
		//	node							  //
		//---------------------------//
		std::vector<TCellID> n0_adj_edges_on_front = m_Front->orderedFrontEdgesAroundNode(m_mesh, n0_id);
		for (auto e_adj_id:n0_adj_edges_on_front)
		{
			if ( e_adj_id != m_e_id
			    && var_front_edges_classification->value(e_adj_id) == 1)
			{
				m_StructManager->setCornerFaceCreated(e_adj_id, n0_id);
				m_StructManager->setNextDiagNode(e_adj_id, n0_id, n0_2_id);

				NodeNeighbourhoodOnFront_3D n0_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n0_id);
				n0_neighbourhood.execute();

				//TCellID next_edge_id_0side = n0_neighbourhood.nextEdgeOfFace(e_front_faces[0], m_e_id) ;
				//TCellID next_edge_id_1side = n0_neighbourhood.nextEdgeOfFace(e_front_faces[1], m_e_id) ;
				//TCellID f_adj_0 = n0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e_adj_id, next_edge_id_0side, next_edge_id_1side);
				//TCellID f_adj_1 = n0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e_adj_id, next_edge_id_1side, next_edge_id_0side);
				TCellID f_adj_0 = n0_neighbourhood.adjFaceToEdge1InFaceSideStartingFromEdge2(e_adj_id, e_front_faces[0], m_e_id);
				TCellID f_adj_1 = n0_neighbourhood.adjFaceToEdge1InFaceSideStartingFromEdge2(e_adj_id, e_front_faces[1], m_e_id);

				m_StructManager->setCornerNextAdjNode(e_adj_id, f_adj_0, n0_id, n0_1_id);
				m_StructManager->setCornerNextAdjNode(e_adj_id, f_adj_1, n0_id, n0_3_id);

				m_StructManager->setFaceNextNode(f_adj_0, n0_id, n0_1_id);
				m_StructManager->setFaceNextNode(f_adj_1, n0_id, n0_3_id);
			}
		}


		//---------------------------//
		// Update the face 			  //
		// information					  //
		//---------------------------//
		m_StructManager->setFaceNextNode(e_front_faces[0], n0_id, n0_1_id);
		m_StructManager->setFaceNextNode(e_front_faces[1], n0_id, n0_3_id);

	}


	// Second side
	TCellID n1_id = e_nodes[1].id();
	TCellID n1_1_id(NullID);
	TCellID n1_2_id(NullID);
	TCellID n1_3_id(NullID);
	if (m_StructManager->isCornerFaceCreated(m_e_id, n1_id))
	{
		// Then, a face has already been inserted on node e_nodes[1] with a node template
		n1_1_id = m_StructManager->getCornerNextAdjNode(m_e_id, e_front_faces[0], n1_id);
		n1_2_id = m_StructManager->getNextDiagNode(m_e_id, n1_id);
		n1_3_id = m_StructManager->getCornerNextAdjNode(m_e_id, e_front_faces[1], n1_id);
	}
	else
	{
		//============================//
		// We create the 3 new nodes	//
		//============================//

		/*
		math::Vector3d f0_normal = m_Front->outgoingNormal(m_mesh, m_meshT, m_fl, m_VectorField, e_front_faces[0]) ;
		f0_normal.normalize();
		f0_normal = e.length()*f0_normal ;
		 */
		math::Vector3d f0_normal = computeNormaltoFacesAroundNodeSideFace(n1_id, e_front_faces[0]);

		/*
		math::Vector3d f1_normal = m_Front->outgoingNormal(m_mesh, m_meshT, m_fl, m_VectorField, e_front_faces[1]) ;
		f1_normal.normalize();
		f1_normal = e.length()*f1_normal ;
		 */
		math::Vector3d f1_normal = computeNormaltoFacesAroundNodeSideFace(n1_id, e_front_faces[1]);

		math::Point p_n1 = e_nodes[1].point();

		/*
		n1_2_id = m_StructManager->getFaceIdealNextNode(e_front_faces[0], e_nodes[1].id());
		Node n1_2_new = m_mesh->get<Node>(n1_2_id);
		math::Point p1_2_new = n1_2_new.point() ;
		 */
		math::Point p1_2_new = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, e_nodes[1].point(), m_dc, m_DistanceField, f0_normal.normalize()+f1_normal.normalize());
		Node n1_2_new = m_mesh->newNode(p1_2_new);
		//n1_2_id = n1_2_new.id();
		math::Point p1_1_new = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, e_nodes[1].point(), m_dc, m_DistanceField, f0_normal.normalize());
		math::Point p1_3_new = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, m_fl, e_nodes[1].point(), m_dc, m_DistanceField, f1_normal.normalize());

		Node n1_1_new = m_mesh->newNode(p1_1_new);
		Node n1_3_new = m_mesh->newNode(p1_3_new);

		n1_1_id = n1_1_new.id();
		n1_2_id = n1_2_new.id();
		n1_3_id = n1_3_new.id();

		//---------------------------//
		// Update the layer ids of   //
		// the new nodes				  //
		//---------------------------//
		var_node_couche_id->set(n1_1_id, m_Front->getFrontID()+1);
		var_node_couche_id->set(n1_2_id, m_Front->getFrontID()+1);
		var_node_couche_id->set(n1_3_id, m_Front->getFrontID()+1);

		//---------------------------//
		// Update the edge corner	  //
		// information around the	  //
		//	node							  //
		//---------------------------//
		std::vector<TCellID> n1_adj_edges_on_front = m_Front->orderedFrontEdgesAroundNode(m_mesh, n1_id);
		for (auto e_adj_id:n1_adj_edges_on_front)
		{
			if ( e_adj_id != m_e_id
			    && var_front_edges_classification->value(e_adj_id) == 1)
			{
				m_StructManager->setCornerFaceCreated(e_adj_id, n1_id);
				m_StructManager->setNextDiagNode(e_adj_id, n1_id, n1_2_id);

				NodeNeighbourhoodOnFront_3D n1_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n1_id);
				n1_neighbourhood.execute();

				//TCellID next_edge_id_0side = n1_neighbourhood.nextEdgeOfFace(e_front_faces[0], m_e_id) ;
				//TCellID next_edge_id_1side = n1_neighbourhood.nextEdgeOfFace(e_front_faces[1], m_e_id) ;
				//TCellID f_adj_0 = n1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e_adj_id, next_edge_id_0side, next_edge_id_1side);
				//TCellID f_adj_1 = n1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e_adj_id, next_edge_id_1side, next_edge_id_0side);
				TCellID f_adj_0 = n1_neighbourhood.adjFaceToEdge1InFaceSideStartingFromEdge2(e_adj_id, e_front_faces[0], m_e_id);
				TCellID f_adj_1 = n1_neighbourhood.adjFaceToEdge1InFaceSideStartingFromEdge2(e_adj_id, e_front_faces[1], m_e_id);

				m_StructManager->setCornerNextAdjNode(e_adj_id, f_adj_0, n1_id, n1_1_id);
				m_StructManager->setCornerNextAdjNode(e_adj_id, f_adj_1, n1_id, n1_3_id);

				m_StructManager->setFaceNextNode(f_adj_0, n1_id, n1_1_id);
				m_StructManager->setFaceNextNode(f_adj_1, n1_id, n1_3_id);
			}
		}


		//---------------------------//
		// Update the face 			  //
		// information					  //
		//---------------------------//
		m_StructManager->setFaceNextNode(e_front_faces[0], n1_id, n1_1_id);
		m_StructManager->setFaceNextNode(e_front_faces[1], n1_id, n1_3_id);

	}

	Node n0_1 = m_mesh->get<Node>(n0_1_id);
	Node n0_2 = m_mesh->get<Node>(n0_2_id);
	Node n0_3 = m_mesh->get<Node>(n0_3_id);
	Node n1_1 = m_mesh->get<Node>(n1_1_id);
	Node n1_2 = m_mesh->get<Node>(n1_2_id);
	Node n1_3 = m_mesh->get<Node>(n1_3_id);

	// Create the hexa
	m_hex.push_back(m_mesh->get<Region>( math::Utils::CreateHexaNConnectivities(m_mesh, e_nodes[0], n0_1, n0_2, n0_3,
	                                                                           e_nodes[1], n1_1, n1_2, n1_3)));

	//---------------------------//
	// Update the two new faces  //
	// layer id                  //
	//---------------------------//
	Variable<int>* var_face_couche_id = m_mesh->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");
	TCellID f_1_new_id = math::Utils::CommonFace3Nodes(m_mesh, n0_1_id, n0_2_id, n1_2_id);
	TCellID f_2_new_id = math::Utils::CommonFace3Nodes(m_mesh, n1_2_id,  n0_2.id(), n0_3.id());
	var_face_couche_id->set(f_1_new_id, m_Front->getFrontID()+1);
	var_face_couche_id->set(f_2_new_id, m_Front->getFrontID()+1);

}
/*------------------------------------------------------------------------*/
math::Vector3d
PatternEdgeCorner::computeNormaltoFacesAroundNodeSideFace(TCellID n_id, TCellID f_id)
{
	Node n = m_mesh->get<Node>(n_id);
	NodeNeighbourhoodOnFront_3D n_neighbourhood = NodeNeighbourhoodOnFront_3D(m_mesh, m_Front, n_id);
	n_neighbourhood.execute();

	Variable<int>* var_front_edges_classification = m_mesh->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");

	// Compute the next corner edge, to know when to stop
	TCellID next_corner_edge_id(m_e_id);
	for (auto e_id:n_neighbourhood.getOrderedEdges())
	{
		if (e_id != m_e_id
		    && var_front_edges_classification->value(e_id)==1)
		{
			next_corner_edge_id = e_id;
		}
	}

	std::vector<TCellID> faces = n_neighbourhood.facesBtwEdge1nEdge2inFaceSide(m_e_id, next_corner_edge_id, f_id) ;

	math::Vector3d normal;
	for (auto face_id:faces)
	{
		normal += m_Front->outgoingNormal(m_mesh, m_meshT, m_fl, m_VectorField, face_id).normalize() ;
	}

	return normal.normalize();

}
/*------------------------------------------------------------------------*/