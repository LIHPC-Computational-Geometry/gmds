//
// Created by rochec on 18/11/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/claire/Front_3D.h>
#include <gmds/claire/Utils.h>
#include <gmds/claire/AeroException.h>
#include <gmds/claire/AeroExtrusion_3D.h>
#include <gmds/claire/AdvectedPointRK4_3D.h>
#include <gmds/claire/FrontEdgesNodesClassification_3D.h>
#include <gmds/claire/NodeNeighbourhoodOnFront_3D.h>

#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <iostream>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AeroExtrusion_3D::AeroExtrusion_3D(Mesh *AMeshT, Mesh *AMeshH, ParamsAero Aparams_aero, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField) :
  m_meshT(AMeshT),
  m_meshH(AMeshH),
  m_fl(m_meshT)
{
	m_meshT = AMeshT;
	m_meshH = AMeshH;
	m_params_aero = Aparams_aero;
	m_DistanceField = A_DistanceField;
	m_VectorField = A_VectorField;
	m_iteration = 0;
}


/*------------------------------------------------------------------------*/
AeroExtrusion_3D::STATUS
AeroExtrusion_3D::execute()
{
	double pas_couche = 1.0/m_params_aero.nbr_couches ;

	std::vector<TCellID> surface_block_corners_Id;
	std::vector<TCellID> surface_block_faces_Id;

	for (auto n_id:m_meshH->nodes())
	{
		surface_block_corners_Id.push_back(n_id);
	}
	for (auto f_id:m_meshH->faces())
	{
		surface_block_faces_Id.push_back(f_id);
	}

	Front_3D Front_Geom = Front_3D(0, surface_block_corners_Id, surface_block_faces_Id);

	Front_3D Current_Front = Compute1stLayer(Front_Geom, m_meshT->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"), m_VectorField);

	// Compute the successive layers
	for (int i=2; i <= m_params_aero.nbr_couches; i++) {
		Current_Front = ComputeLayer(Current_Front, m_meshT->getVariable<double,GMDS_NODE>("GMDS_Distance"), i*pas_couche,
		                             m_VectorField);
	}

	Variable<int>* var_node_layer = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* var_face_layer = m_meshH->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");
	for (auto f_id:m_meshH->faces())
	{
		Face f = m_meshH->get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		int max_layer_index = std::max(var_node_layer->value(nodes[0].id()), var_node_layer->value(nodes[1].id()));
		max_layer_index = std::max(max_layer_index, var_node_layer->value(nodes[2].id()));
		max_layer_index = std::max(max_layer_index, var_node_layer->value(nodes[3].id()));
		var_face_layer->set(f_id, max_layer_index);
	}

	return AeroExtrusion_3D::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
std::map<TCellID, TCellID>
AeroExtrusion_3D::ComputeIdealPositions(Front_3D AFront, double dist_cible, Variable<double>* A_distance, Variable<math::Vector3d>* A_vectors)
{
	std::map<TCellID, TCellID> map_idealNextNodes;
	std::vector<TCellID> front_nodes = AFront.getNodes();

	for (auto n_id:front_nodes){
		Node n = m_meshH->get<Node>(n_id);
		math::Point M = n.point();
		AdvectedPointRK4_3D advpoint(m_meshT, &m_fl, M, dist_cible, A_distance, A_vectors);
		advpoint.execute();
		Node n_new = m_meshH->newNode(advpoint.getPend());
		map_idealNextNodes[n_id] = n_new.id() ;
	}

	return map_idealNextNodes;
}
/*------------------------------------------------------------------------*/
Front_3D
AeroExtrusion_3D::ComputeLayer(Front_3D Front_IN, Variable<double>* A_distance, double dist_cible, Variable<math::Vector3d>* A_vectors){

	std::cout << "---------> build layer: " << Front_IN.getFrontID()+1 << std::endl;

	std::map<TCellID, TCellID> map_new_nodes = ComputeIdealPositions(Front_IN, dist_cible, A_distance, A_vectors);

	// Init the face_info type for the faces of the front
	InitFaceStructInfo(Front_IN, map_new_nodes);

	// Variables
	Variable<int>* var_NODE_couche_id = m_meshH->getVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* var_front_edges_classification = m_meshH->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
	Variable<int>* var_front_nodes_classification = m_meshH->getOrCreateVariable<int, GMDS_NODE>("Nodes_Classification");


	// Mise à jour de l'indice de couche
	for (auto n_id:Front_IN.getNodes()){
		var_NODE_couche_id->set(map_new_nodes[n_id], Front_IN.getFrontID()+1);
	}

	// Front edges and nodes classification
	FrontEdgesNodesClassification_3D Classification = FrontEdgesNodesClassification_3D(m_meshH,
	                                                                                   &Front_IN,
	                                                                                   var_front_edges_classification,
	                                                                                   var_front_nodes_classification);
	Classification.execute();

	// Init the edge_info type for the edges of the front
	InitEdgeStructInfo(Front_IN);

	int mark_edgesTreated = m_meshH->newMark<Edge>();
	int mark_facesTreated = m_meshH->newMark<Face>();

	int mark_EdgesTemplates = Classification.getMarkEdgesTemplates();
	int mark_NodesTemplates = Classification.getMarkNodesTemplates();

	// Create the HEXA on each NODE classified
	std::map<TCellID, int> singular_nodes = getSingularNodes(Front_IN, var_front_edges_classification);
	std::cout << "Noeuds singuliers: " << singular_nodes.size() << std::endl;
	for (auto singu_node:singular_nodes)
	{
		TCellID n_id = singu_node.first ;
		int singu_type = singu_node.second;
		if (singu_type==1 && m_meshH->isMarked(m_meshH->get<Node>(n_id),mark_NodesTemplates))
		{
			TCellID r_id = TemplateNode3Corner(Front_IN, n_id, map_new_nodes, dist_cible);
		}
		else if (singu_type==2 && m_meshH->isMarked(m_meshH->get<Node>(n_id),mark_NodesTemplates))
		{
			TCellID r_id = TemplateNode2Corner1End(Front_IN, n_id, dist_cible, mark_edgesTreated, mark_facesTreated);
		}
		else if (singu_type==4 && m_meshH->isMarked(m_meshH->get<Node>(n_id),mark_NodesTemplates))
		{
			TCellID r_id = TemplateNode3End(Front_IN, n_id, mark_edgesTreated, mark_facesTreated);
		}
	}


	// Create the HEXA on each EDGE classified
	std::map<TCellID, int> singular_edges = getSingularEdges(Front_IN, var_front_edges_classification, mark_edgesTreated);
	std::cout << "Arêtes singuliéres: " << singular_edges.size() << std::endl;
	for (auto singu_edge:singular_edges)
	{
		TCellID e_id = singu_edge.first ;
		int singu_type = singu_edge.second;
		if (singu_type==1
		    && !m_meshH->isMarked(m_meshH->get<Edge>(e_id), mark_edgesTreated)
		    && m_meshH->isMarked(m_meshH->get<Edge>(e_id), mark_EdgesTemplates))
		{
			TCellID r_id = TemplateEdgeCorner(Front_IN, e_id, dist_cible);
		}
		else if (singu_type==2
		         && !m_meshH->isMarked(m_meshH->get<Edge>(e_id), mark_edgesTreated)
		         && m_meshH->isMarked(m_meshH->get<Edge>(e_id), mark_EdgesTemplates))
		{
			TCellID r_id = TemplateEdgeEnd(Front_IN, e_id, dist_cible, mark_edgesTreated, mark_facesTreated);
		}
		else
		{
			//std::cout << "Edge singularity not implemented yet" << std::endl;
		}
	}


	// Create the HEXA on each FACE of the front
	for (auto f_id:Front_IN.getFaces()){
		Face f = m_meshH->get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		if (!m_meshH->isMarked(m_meshH->get<Face>(f_id), mark_facesTreated))
		{
			TemplateFace(f_id, Front_IN, map_new_nodes);
		}
	}

	m_meshH->unmarkAll<Edge>(mark_edgesTreated);
	m_meshH->freeMark<Edge>(mark_edgesTreated);
	m_meshH->unmarkAll<Face>(mark_facesTreated);
	m_meshH->freeMark<Face>(mark_facesTreated);

	// Erase the nodes connected to nothing
	math::Utils::MeshCleaner(m_meshH);

	// Init the Front OUT
	Front_3D Front_OUT = InitFrontOUT(Front_IN);

	std::cout << "Dist layer: " << dist_cible << std::endl;
	for (auto n_id:m_meshH->nodes())
	{
		if (var_NODE_couche_id->value(n_id) == Front_OUT.getFrontID())
		{
			gmds::Cell::Data data = m_fl.find(m_meshH->get<Node>(n_id).point());
			//std::cout << "Dist: " << m_DistanceField->value(data.id) << std::endl;
		}
	}

	//===================//
	// FAST ANALYSIS		//
	//===================//
	for (auto e_id:m_meshH->edges())
	{
		Edge e = m_meshH->get<Edge>(e_id);
		if (var_NODE_couche_id->value(e.get<Node>()[0].id()) == Front_OUT.getFrontID()
		    && var_NODE_couche_id->value(e.get<Node>()[1].id()) == Front_OUT.getFrontID()
		    && Front_OUT.edgeFacesOnFront(m_meshH, e_id).size() != 2)
		{
			std::cout << "ATTENTION: the layer is not valid. At least one edge of the front has not 2 faces." << std::endl;
		}
	}

	return Front_OUT;
}
/*------------------------------------------------------------------------*/
Front_3D
AeroExtrusion_3D::Compute1stLayer(Front_3D A_Front_IN, Variable<double>* A_distance, Variable<math::Vector3d>* A_vectors)
{
	std::cout << "---------> build layer: " << A_Front_IN.getFrontID()+1 << std::endl;

	std::map<TCellID, TCellID> map_new_nodes = ComputeIdealPositions(A_Front_IN, m_params_aero.delta_cl, A_distance, A_vectors);

	// Init the face_info type
	InitFaceStructInfo(A_Front_IN, map_new_nodes);

	// Mise à jour de l'indice de couche
	Variable<int>*var_NODE_couche_id = m_meshH->getVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	for (auto n_id:A_Front_IN.getNodes()){
		var_NODE_couche_id->set(map_new_nodes[n_id], A_Front_IN.getFrontID()+1);
	}

	// Create the HEXA on each FACE of the front
	for (auto f_id:A_Front_IN.getFaces()){
		Face f = m_meshH->get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		TemplateFace(f_id, A_Front_IN, map_new_nodes);
	}

	// Erase the nodes connected to nothing
	math::Utils::MeshCleaner(m_meshH);

	// Init the Front OUT
	Front_3D Front_OUT = InitFrontOUT(A_Front_IN);

	//===================//
	// FAST ANALYSIS		//
	//===================//
	bool isHexMeshValid = math::Utils::isThisHexMeshValid(m_meshH);
	if (!isHexMeshValid)
	{
		std::cout << "ATTENTION: the mesh is not valid. An element may be reversed during the advancing front." << std::endl;
	}


	return Front_OUT;

}
/*------------------------------------------------------------------------*/
void
AeroExtrusion_3D::InitFaceStructInfo(Front_3D &Front, std::map<TCellID, TCellID> map_new_nodes)
{
	for (auto f_id:Front.getFaces())
	{
		Face f = m_meshH->get<Face>(f_id);
		std::vector<Node> f_nodes = f.get<Node>();

		m_FaceInfo[f_id].next_ideal_nodes[f_nodes[0].id()] = map_new_nodes[f_nodes[0].id()];
		m_FaceInfo[f_id].next_ideal_nodes[f_nodes[1].id()] = map_new_nodes[f_nodes[1].id()];
		m_FaceInfo[f_id].next_ideal_nodes[f_nodes[2].id()] = map_new_nodes[f_nodes[2].id()];
		m_FaceInfo[f_id].next_ideal_nodes[f_nodes[3].id()] = map_new_nodes[f_nodes[3].id()];

		m_FaceInfo[f_id].next_nodes[f_nodes[0].id()] = map_new_nodes[f_nodes[0].id()];
		m_FaceInfo[f_id].next_nodes[f_nodes[1].id()] = map_new_nodes[f_nodes[1].id()];
		m_FaceInfo[f_id].next_nodes[f_nodes[2].id()] = map_new_nodes[f_nodes[2].id()];
		m_FaceInfo[f_id].next_nodes[f_nodes[3].id()] = map_new_nodes[f_nodes[3].id()];

	}
}
/*------------------------------------------------------------------------*/
void
AeroExtrusion_3D::InitEdgeStructInfo(Front_3D &Front)
{
	Variable<int>* var_front_edges_classification = m_meshH->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
	Variable<int>* var_node_couche_id = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");

	for (auto e_id:m_meshH->edges())
	{
		Edge e = m_meshH->get<Edge>(e_id);
		std::vector<Node> e_nodes = e.get<Node>() ;
		if (var_node_couche_id->value(e_nodes[0].id()) == Front.getFrontID()
		    && var_node_couche_id->value(e_nodes[1].id()) == Front.getFrontID())
		{
			// So, the edge is on the front.
			if (var_front_edges_classification->value(e_id)==1)
			{
				// Init Corner Edge information
				m_EdgeInfo[e_id].singularity_type = 1;
				m_EdgeInfo[e_id].CORNER_n_face_created[e_nodes[0].id()] = false;
				m_EdgeInfo[e_id].CORNER_n_face_created[e_nodes[1].id()] = false;
			}
			else if (var_front_edges_classification->value(e_id)==2)
			{
				// Init Corner Edge information
				m_EdgeInfo[e_id].singularity_type = 2;
				m_EdgeInfo[e_id].END_n_face_created[e_nodes[0].id()] = false;
				m_EdgeInfo[e_id].END_n_face_created[e_nodes[1].id()] = false;
			}
		}
	}

}
/*------------------------------------------------------------------------*/
Front_3D
AeroExtrusion_3D::InitFrontOUT(Front_3D &Front_IN)
{
	std::vector<TCellID> new_front_nodes_id;
	std::vector<TCellID> new_front_faces_id;

	Variable<int>* var_NODE_couche_id = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>*var_FACE_couche_id = m_meshH->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");

	for (auto n_id:m_meshH->nodes())
	{
		if (var_NODE_couche_id->value(n_id) == Front_IN.getFrontID()+1)
		{
			new_front_nodes_id.push_back(n_id);
		}
	}
	for (auto f_id:m_meshH->faces())
	{
		if (var_FACE_couche_id->value(f_id) == Front_IN.getFrontID()+1)
		{
			new_front_faces_id.push_back(f_id);
		}
	}

	return Front_3D(Front_IN.getFrontID()+1, new_front_nodes_id, new_front_faces_id);

}
/*------------------------------------------------------------------------*/
std::map<TCellID, int>
AeroExtrusion_3D::getSingularNodes(Front_3D &AFront, Variable<int>* front_edges_classification)
{
	std::map<TCellID, int> sing_nodes;
	for (auto n_id:AFront.getNodes())
	{
		std::vector<TCellID> n_ordered_edges = AFront.orderedFrontEdgesAroundNode(m_meshH, n_id);
		int compteur_corner(0);
		int compteur_end(0);
		int compteur_reversal(0);
		for (auto e_id:n_ordered_edges)
		{
			Edge e = m_meshH->get<Edge>(e_id);
			if (front_edges_classification->value(e_id) == 1)
			{
				compteur_corner +=1;
			}
			else if (front_edges_classification->value(e_id) == 2)
			{
				compteur_end +=1;
			}
			else if (front_edges_classification->value(e_id) == 3)
			{
				compteur_reversal +=1;
			}
		}

		if ( compteur_corner == 0 && compteur_end == 0 && compteur_reversal == 0)
		{
			// The node is regular
		}
		else if (compteur_corner == 3 && compteur_end == 0 && compteur_reversal == 0)
		{
			sing_nodes[n_id] = 1;
		}
		else if (compteur_corner==2 && compteur_end==1)
		{
			sing_nodes[n_id] = 2;
		}
		else if (compteur_end == 3 && n_ordered_edges.size()==3)
		{
			sing_nodes[n_id] = 4;
		}
		else
		{
			/*
			std::cout << "-------" << std::endl;
			std::cout << "TemplateNode pas encore implémenté pour gérer cette configuration." << std::endl;
			std::cout << "NOEUD: " << n_id << std::endl;
			std::cout << "ARÊTES: " << n_ordered_edges.size() << std::endl;
			std::cout << "Corner: " << compteur_corner << std::endl;
			std::cout << "End: " << compteur_end << std::endl;
			std::cout << "Reversal: " << compteur_reversal << std::endl;
			*/
		}

	}
	return sing_nodes;
}
/*------------------------------------------------------------------------*/
std::map<TCellID, int>
AeroExtrusion_3D::getSingularEdges(Front_3D &AFront, Variable<int>* front_edges_classification, int mark_singEdgesTreated)
{
	std::map<TCellID, int> sing_edges;
	Variable<int>* var_node_couche_id = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	for (auto e_id:m_meshH->edges())
	{
		Edge e = m_meshH->get<Edge>(e_id);
		std::vector<Node> e_nodes = e.get<Node>();
		if (var_node_couche_id->value(e_nodes[0].id()) == AFront.getFrontID()
		    && var_node_couche_id->value(e_nodes[1].id()) == AFront.getFrontID())
		{
			if ( !m_meshH->isMarked(e, mark_singEdgesTreated)
			    && front_edges_classification->value(e_id) > 0)
			{
				sing_edges[e_id] = front_edges_classification->value(e_id);
			}
			if (front_edges_classification->value(e_id) > 1)
			{
				//std::cout << "TemplateEdge pas encore implémenté pour gérer cette configuration." << std::endl;
			}
		}

	}
	return sing_edges;
}
/*------------------------------------------------------------------------*/
TCellID
AeroExtrusion_3D::TemplateNode3Corner(Front_3D &AFront, TCellID n_id, std::map<TCellID, TCellID> map_new_nodes, double dc)
{
	//std::cout << "Template Node 3 Corner au noeud " << n_id << std::endl;
	TCellID r_id;

	NodeNeighbourhoodOnFront_3D n_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n_id);
	n_neighbourhood.execute();
	std::vector<TCellID> n_ordered_edges = n_neighbourhood.getOrderedEdges();
	Variable<int>* edge_classification = m_meshH->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
	std::vector<TCellID> n_edges_corner;
	for (auto edge:n_ordered_edges)
	{
		if (edge_classification->value(edge) == 1)
		{
			n_edges_corner.push_back(edge);
		}
	}

	Edge e_adj_1 = m_meshH->get<Edge>(n_edges_corner[0]);
	Edge e_adj_2 = m_meshH->get<Edge>(n_edges_corner[1]);
	Edge e_adj_3 = m_meshH->get<Edge>(n_edges_corner[2]);

	Node n = m_meshH->get<Node>(n_id);
	Node n6 = m_meshH->get<Node>(map_new_nodes[n_id]);

	Node n_adj_1 = e_adj_1.getOppositeNode(n);
	Node n_adj_2 = e_adj_2.getOppositeNode(n);
	Node n_adj_3 = e_adj_3.getOppositeNode(n);

	math::Vector3d v1 = n_adj_1.point()-n.point() ;
	math::Vector3d v2 = n_adj_2.point()-n.point() ;
	math::Vector3d v3 = n_adj_3.point()-n.point() ;

	// Compute positions of new nodes
	v1.normalize();
	v2.normalize();
	v3.normalize();
	math::Point p1 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, -v1);
	math::Point p3 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, -v2);
	math::Point p4 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, -v3);
	math::Point p2 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, -v1-v2); // (p3-n.point()) + (p1-n.point()) ;
	math::Point p5 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, -v1-v3); // p4 + (p1-n.point()) ;
	math::Point p7 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, -v2-v3); //p4 + (p3-n.point()) ;

	Node n1 = m_meshH->newNode(p1);
	Node n2 = m_meshH->newNode(p2);
	Node n3 = m_meshH->newNode(p3);
	Node n4 = m_meshH->newNode(p4);
	Node n5 = m_meshH->newNode(p5);
	Node n7 = m_meshH->newNode(p7);

	// Update the layer ids
	Variable<int>* var_node_couche_id = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	var_node_couche_id->set(n1.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n2.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n3.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n4.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n5.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n7.id(), AFront.getFrontID()+1);

	r_id = math::Utils::CreateHexaNConnectivities(m_meshH, n, n1, n2, n3, n4, n5, n6, n7);

	// Update the faces information of the front
	std::vector<TCellID> adj_faces = n_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(e_adj_1.id(), e_adj_2.id(), e_adj_3.id());
	for (auto f_adj_id:adj_faces)
	{
		m_FaceInfo[f_adj_id].next_nodes[n_id] = n4.id() ;
	}
	adj_faces = n_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(e_adj_2.id(), e_adj_3.id(), e_adj_1.id());
	for (auto f_adj_id:adj_faces)
	{
		m_FaceInfo[f_adj_id].next_nodes[n_id] = n1.id() ;
	}
	adj_faces = n_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(e_adj_1.id(), e_adj_3.id(), e_adj_2.id());
	for (auto f_adj_id:adj_faces)
	{
		m_FaceInfo[f_adj_id].next_nodes[n_id] = n3.id() ;
	}

	TCellID f_1 = math::Utils::CommonFace3Nodes(m_meshH, n_id, n_adj_1.id(), n_adj_2.id());
	TCellID f_2 = math::Utils::CommonFace3Nodes(m_meshH, n_id,  n_adj_2.id(), n_adj_3.id());
	TCellID f_3 = math::Utils::CommonFace3Nodes(m_meshH, n_id,  n_adj_3.id(), n_adj_1.id());
	/*
	m_FaceInfo[f_1].next_nodes[n_id] = n4.id() ;
	m_FaceInfo[f_2].next_nodes[n_id] = n1.id() ;
	m_FaceInfo[f_3].next_nodes[n_id] = n3.id() ;
	 */

	//---------------------------//
	// Update the corner edges   //
	// information					  //
	//---------------------------//
	TCellID e1_id = math::Utils::CommonEdge(m_meshH, n.id(), n_adj_1.id());
	TCellID e2_id = math::Utils::CommonEdge(m_meshH, n.id(), n_adj_2.id());
	TCellID e3_id = math::Utils::CommonEdge(m_meshH, n.id(), n_adj_3.id());

	m_EdgeInfo[e1_id].CORNER_n_face_created[n.id()] = true ;
	m_EdgeInfo[e2_id].CORNER_n_face_created[n.id()] = true ;
	m_EdgeInfo[e3_id].CORNER_n_face_created[n.id()] = true ;

	std::pair<TCellID, TCellID> pair_1(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e1_id, e2_id, e3_id), n.id());
	std::pair<TCellID, TCellID> pair_2(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e1_id, e3_id, e2_id), n.id());
	std::pair<TCellID, TCellID> pair_3(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e2_id, e1_id, e3_id), n.id());
	std::pair<TCellID, TCellID> pair_4(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e2_id, e3_id, e1_id), n.id());
	std::pair<TCellID, TCellID> pair_5(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e3_id, e2_id, e1_id), n.id());
	std::pair<TCellID, TCellID> pair_6(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e3_id, e1_id, e2_id), n.id());


	m_EdgeInfo[e1_id].CORNER_next_nodes[pair_1] = n4.id();
	m_EdgeInfo[e1_id].CORNER_next_nodes[pair_2] = n3.id();

	m_EdgeInfo[e2_id].CORNER_next_nodes[pair_3] = n4.id();
	m_EdgeInfo[e2_id].CORNER_next_nodes[pair_4] = n1.id();

	m_EdgeInfo[e3_id].CORNER_next_nodes[pair_5] = n1.id();
	m_EdgeInfo[e3_id].CORNER_next_nodes[pair_6] = n3.id();

	m_EdgeInfo[e1_id].diag_next_node[n.id()] = n7.id() ;
	m_EdgeInfo[e2_id].diag_next_node[n.id()] = n5.id() ;
	m_EdgeInfo[e3_id].diag_next_node[n.id()] = n2.id() ;

	//---------------------------//
	// Update the face layer id  //
	//---------------------------//
	Variable<int>* var_face_couche_id = m_meshH->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");
	TCellID f_1_new_id = math::Utils::CommonFace3Nodes(m_meshH, n4.id(), n5.id(), n6.id());
	TCellID f_2_new_id = math::Utils::CommonFace3Nodes(m_meshH, n1.id(),  n2.id(), n5.id());
	TCellID f_3_new_id = math::Utils::CommonFace3Nodes(m_meshH, n2.id(),  n3.id(), n6.id());
	var_face_couche_id->set(f_1_new_id, AFront.getFrontID()+1);
	var_face_couche_id->set(f_2_new_id, AFront.getFrontID()+1);
	var_face_couche_id->set(f_3_new_id, AFront.getFrontID()+1);

	gmds::IGMeshIOService ioService(m_meshH);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("AeroExtrusion_3D_"+std::to_string(m_iteration)+".vtk");
	m_iteration++;

	return r_id;
}
/*------------------------------------------------------------------------*/
TCellID
AeroExtrusion_3D::TemplateNode2Corner1End(Front_3D &AFront, TCellID n_id, double dc, int mark_edgesTreated, int mark_facesTreated)
{
	//std::cout << "Template Node 2 Corner 1 End au noeud " << n_id << std::endl;
	TCellID r_id;

	Node n = m_meshH->get<Node>(n_id);

	NodeNeighbourhoodOnFront_3D n_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n_id);
	n_neighbourhood.execute();

	Variable<int>* var_front_edges_classification = m_meshH->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");

	std::vector<Edge> ec;	// For the two edge CORNER
	Edge ee;						// For the edge END

	// Get the local feature edges
	for (auto edge_id:n_neighbourhood.getOrderedEdges())
	{
		if (var_front_edges_classification->value(edge_id)==1)
		{
			ec.push_back(m_meshH->get<Edge>(edge_id));
		}
		else if (var_front_edges_classification->value(edge_id)==2)
		{
			ee = m_meshH->get<Edge>(edge_id);
		}
	}

	Node n_c0 = ec[0].getOppositeNode(n);
	Node n_c1 = ec[1].getOppositeNode(n);

	TCellID f_end0_id = math::Utils::CommonFace3Nodes(m_meshH, n_c0.id(), n_id, ee.getOppositeNode(n).id()) ;
	TCellID f_end1_id = math::Utils::CommonFace3Nodes(m_meshH, n_c1.id(), n_id, ee.getOppositeNode(n).id()) ;

	math::Vector3d ve = ( ( n_c0.point() - n.point() ).normalize() + ( n_c1.point() - n.point() ).normalize() ).normalize() ;
	math::Point pe = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, ve);
	Node n_e = m_meshH->newNode(pe);

	math::Vector3d v1 = (n.point()- (ee.getOppositeNode(n)).point() ).normalize() ;
	math::Point p1 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, v1);
	Node n_1 = m_meshH->newNode(p1);

	math::Point p2 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n_c0.point(), dc, m_DistanceField, v1);
	Node n_2 = m_meshH->newNode(p2);

	math::Point p4 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n_c1.point(), dc, m_DistanceField, v1);
	Node n_4 = m_meshH->newNode(p4);

	math::Vector3d v3 = ( ( n_c0.point() - n.point() ).normalize() + ( n_c1.point() - n.point() ).normalize() + v1 ).normalize() ;
	math::Point p3 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, v3);
	Node n_3 = m_meshH->newNode(p3);

	// Create the new hexa
	r_id = math::Utils::CreateHexaNConnectivities(m_meshH, n, n_c0, n_e, n_c1, n_1, n_2, n_3, n_4);

	// Update the layer ids
	Variable<int>* var_node_couche_id = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	var_node_couche_id->set(n_1.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n_2.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n_3.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n_4.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n_e.id(), AFront.getFrontID()+1);
	Variable<int>* var_face_couche_id = m_meshH->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");
	var_face_couche_id->set(math::Utils::CommonFace3Nodes(m_meshH, n_1.id(), n_2.id(), n_3.id()), AFront.getFrontID()+1);

	// Update the faces and edges that can't do actions anymore
	m_meshH->mark(ec[0], mark_edgesTreated);
	m_meshH->mark(ec[1], mark_edgesTreated);

	// Update the FaceInfo structures
	for (auto face_id:n_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[0].id(), ec[1].id(), ee.id()))
	{
		m_FaceInfo[face_id].next_nodes[n_id] = n_1.id() ;
	}

	m_FaceInfo[n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec[0].id(), ec[1].id(), ee.id())].next_nodes[n_c0.id()] = n_2.id();
	m_FaceInfo[n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec[1].id(), ec[0].id(), ee.id())].next_nodes[n_c1.id()] = n_4.id();

	TCellID f_e0_id = math::Utils::CommonFace3Nodes(m_meshH, n_c0.id(), n_id, ee.getOppositeNode(n).id());
	m_FaceInfo[f_e0_id].next_nodes[n_c0.id()] = n_e.id();
	TCellID f_e1_id = math::Utils::CommonFace3Nodes(m_meshH, n_c1.id(), n_id, ee.getOppositeNode(n).id());
	m_FaceInfo[f_e1_id].next_nodes[n_c1.id()] = n_e.id();


	// Update the EdgeInfo structure for the END edge
	m_EdgeInfo[ee.id()].END_n_face_created[n_id] = true;
	m_EdgeInfo[ee.id()].diag_next_node[n_id] = n_e.id();

	// Update the EdgeInfo structure for the two CORNER edges
	NodeNeighbourhoodOnFront_3D n_c0_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n_c0.id()) ;
	n_c0_neighbourhood.execute();
	for (auto edge_id:n_c0_neighbourhood.getOrderedEdges())
	{
		if (var_front_edges_classification->value(edge_id) == 1
		    && edge_id != ec[0].id())
		{
			m_EdgeInfo[edge_id].CORNER_n_face_created[n_c0.id()] = true;
			m_EdgeInfo[edge_id].diag_next_node[n_c0.id()] = n_3.id();

			TCellID edge_f_e0 = n_c0_neighbourhood.nextEdgeOfFace(f_e0_id, ec[0].id()) ;
			for (auto faces_id:n_c0_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[0].id(), edge_id, edge_f_e0))
			{
				m_FaceInfo[faces_id].next_nodes[n_c0.id()] = n_2.id();
			}

			std::pair<TCellID, TCellID> pair_1(n_c0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[0].id(), edge_f_e0), n_c0.id());
			std::pair<TCellID, TCellID> pair_2(n_c0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, edge_f_e0, ec[0].id()), n_c0.id());

			m_EdgeInfo[edge_id].CORNER_next_nodes[pair_1] = n_2.id();
			m_EdgeInfo[edge_id].CORNER_next_nodes[pair_2] = n_e.id();

		}
	}

	NodeNeighbourhoodOnFront_3D n_c1_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n_c1.id()) ;
	n_c1_neighbourhood.execute();
	for (auto edge_id:n_c1_neighbourhood.getOrderedEdges())
	{
		if (var_front_edges_classification->value(edge_id) == 1
		    && edge_id != ec[1].id())
		{
			m_EdgeInfo[edge_id].CORNER_n_face_created[n_c1.id()] = true;
			m_EdgeInfo[edge_id].diag_next_node[n_c1.id()] = n_3.id();

			TCellID edge_f_e1 = n_c1_neighbourhood.nextEdgeOfFace(f_e1_id, ec[1].id()) ;
			for (auto faces_id:n_c1_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[1].id(), edge_id, edge_f_e1))
			{
				m_FaceInfo[faces_id].next_nodes[n_c1.id()] = n_4.id();
			}

			std::pair<TCellID, TCellID> pair_1(n_c1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[1].id(), edge_f_e1), n_c1.id());
			std::pair<TCellID, TCellID> pair_2(n_c1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, edge_f_e1, ec[1].id()), n_c1.id());

			m_EdgeInfo[edge_id].CORNER_next_nodes[pair_1] = n_4.id();
			m_EdgeInfo[edge_id].CORNER_next_nodes[pair_2] = n_e.id();

		}
	}

	gmds::IGMeshIOService ioService(m_meshH);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("AeroExtrusion_3D_"+std::to_string(m_iteration)+".vtk");
	m_iteration++;

	return r_id;
}
/*------------------------------------------------------------------------*/
TCellID
AeroExtrusion_3D::TemplateNode3End(Front_3D &AFront, TCellID n_id, int mark_edgesTreated, int mark_facesTreated)
{
	//std::cout << "Template Node 3 End au noeud " << n_id << std::endl;
	TCellID r_id;

	Node n = m_meshH->get<Node>(n_id);

	NodeNeighbourhoodOnFront_3D n_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n_id);
	n_neighbourhood.execute();
	std::vector<TCellID> n_ordered_edges = n_neighbourhood.getOrderedEdges();

	Edge e_0 = m_meshH->get<Edge>(n_ordered_edges[0]);
	Edge e_1 = m_meshH->get<Edge>(n_ordered_edges[1]);
	Edge e_2 = m_meshH->get<Edge>(n_ordered_edges[2]);

	Node n_e0_opp = e_0.getOppositeNode(n);
	Node n_e1_opp = e_1.getOppositeNode(n);
	Node n_e2_opp = e_2.getOppositeNode(n);

	TCellID f0_id = math::Utils::CommonFace3Nodes(m_meshH, n_id, n_e0_opp.id(), n_e1_opp.id());
	TCellID f1_id = math::Utils::CommonFace3Nodes(m_meshH, n_id, n_e1_opp.id(), n_e2_opp.id());
	TCellID f2_id = math::Utils::CommonFace3Nodes(m_meshH, n_id, n_e2_opp.id(), n_e0_opp.id());

	Face f0 = m_meshH->get<Face>(f0_id);
	Face f1 = m_meshH->get<Face>(f1_id);
	Face f2 = m_meshH->get<Face>(f2_id);

	std::vector<Node> f0_nodes = f0.get<Node>();
	std::vector<Node> f1_nodes = f1.get<Node>();
	std::vector<Node> f2_nodes = f2.get<Node>();

	Node n_f0_diag;
	Node n_f1_diag;
	Node n_f2_diag;
	for (auto n_loc:f0_nodes)
	{
		if (n_loc.id() != n_id && n_loc.id() != n_e0_opp.id() && n_loc.id() != n_e1_opp.id())
		{
			n_f0_diag = n_loc;
		}
	}
	for (auto n_loc:f1_nodes)
	{
		if (n_loc.id() != n_id && n_loc.id() != n_e1_opp.id() && n_loc.id() != n_e2_opp.id())
		{
			n_f1_diag = n_loc;
		}
	}
	for (auto n_loc:f2_nodes)
	{
		if (n_loc.id() != n_id && n_loc.id() != n_e2_opp.id() && n_loc.id() != n_e0_opp.id())
		{
			n_f2_diag = n_loc;
		}
	}

	Node n_new = m_meshH->get<Node>(m_FaceInfo[f0_id].next_ideal_nodes[n_id]);

	// Create the new hexa
	r_id = math::Utils::CreateHexaNConnectivities(m_meshH, n, n_e0_opp, n_f0_diag, n_e1_opp, n_e2_opp, n_f2_diag, n_new, n_f1_diag);

	// Update the faces and edges that can't do actions anymore
	m_meshH->mark(e_0, mark_edgesTreated);
	m_meshH->mark(e_1, mark_edgesTreated);
	m_meshH->mark(e_2, mark_edgesTreated);
	m_meshH->mark(f0, mark_facesTreated);
	m_meshH->mark(f1, mark_facesTreated);
	m_meshH->mark(f2, mark_facesTreated);

	// Update the edges infos
	Variable<int>* var_front_edges_classification = m_meshH->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
	NodeNeighbourhoodOnFront_3D n_e0_opp_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n_e0_opp.id()) ;
	n_e0_opp_neighbourhood.execute();
	NodeNeighbourhoodOnFront_3D n_e1_opp_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n_e1_opp.id()) ;
	n_e1_opp_neighbourhood.execute();
	NodeNeighbourhoodOnFront_3D n_e2_opp_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n_e2_opp.id()) ;
	n_e2_opp_neighbourhood.execute();
	int compteur(0);
	for (auto e_loc_id:n_e0_opp_neighbourhood.getOrderedEdges())
	{
		if (var_front_edges_classification->value(e_loc_id) == 2 && e_loc_id != e_0.id())
		{
			m_EdgeInfo[e_loc_id].diag_next_node[n_e0_opp.id()] = n_new.id() ;
			m_EdgeInfo[e_loc_id].END_n_face_created[n_e0_opp.id()] = true ;
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
			m_EdgeInfo[e_loc_id].diag_next_node[n_e1_opp.id()] = n_new.id() ;
			m_EdgeInfo[e_loc_id].END_n_face_created[n_e1_opp.id()] = true ;
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
			m_EdgeInfo[e_loc_id].diag_next_node[n_e2_opp.id()] = n_new.id() ;
			m_EdgeInfo[e_loc_id].END_n_face_created[n_e2_opp.id()] = true ;
			compteur++;
		}
	}
	if (compteur != 1)
	{
		std::cout << "Attention." << std::endl;
	}


	// Update the FaceInfo structures
	NodeNeighbourhoodOnFront_3D node_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n_f0_diag.id());
	node_neighbourhood.execute();
	for (auto face_id:node_neighbourhood.getOrderedFaces())
	{
		m_FaceInfo[face_id].next_nodes[n_f0_diag.id()] = n_new.id();
	}
	node_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n_f1_diag.id());
	node_neighbourhood.execute();
	for (auto face_id:node_neighbourhood.getOrderedFaces())
	{
		m_FaceInfo[face_id].next_nodes[n_f1_diag.id()] = n_new.id();
	}
	node_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n_f2_diag.id());
	node_neighbourhood.execute();
	for (auto face_id:node_neighbourhood.getOrderedFaces())
	{
		m_FaceInfo[face_id].next_nodes[n_f2_diag.id()] = n_new.id();
	}

	gmds::IGMeshIOService ioService(m_meshH);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("AeroExtrusion_3D_"+std::to_string(m_iteration)+".vtk");
	m_iteration++;

	return r_id;
}
/*------------------------------------------------------------------------*/
TCellID
AeroExtrusion_3D::TemplateEdgeCorner(Front_3D &AFront, TCellID e_id, double dc)
{
	//std::cout << "Template Edge Corner" << std::endl;
	TCellID r_id(NullID);

	Edge e = m_meshH->get<Edge>(e_id);
	std::vector<Node> e_nodes = e.get<Node>();

	std::vector<TCellID> e_front_faces = AFront.edgeFacesOnFront(m_meshH, e_id) ;

	Variable<int>* var_node_couche_id = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* var_front_edges_classification = m_meshH->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");

	// First side
	TCellID n0_id = e_nodes[0].id();
	TCellID n0_1_id(NullID);
	TCellID n0_2_id(NullID);
	TCellID n0_3_id(NullID);
	if (m_EdgeInfo[e_id].CORNER_n_face_created[n0_id])
	{
		// Then, a face has already been inserted on node e_nodes[0] with a node template
		std::pair<TCellID, TCellID> pair_1(e_front_faces[0], n0_id);
		n0_1_id = m_EdgeInfo[e_id].CORNER_next_nodes[pair_1] ;
		n0_2_id = m_EdgeInfo[e_id].diag_next_node[n0_id];
		std::pair<TCellID, TCellID> pair_2(e_front_faces[1], n0_id);
		n0_3_id = m_EdgeInfo[e_id].CORNER_next_nodes[pair_2] ;
	}
	else
	{
		//============================//
		// We create the 3 new nodes	//
		//============================//

		math::Vector3d f0_normal = AFront.outgoingNormal(m_meshH, e_front_faces[0]) ;
		f0_normal.normalize();
		f0_normal = e.length()*f0_normal ;

		math::Vector3d f1_normal = AFront.outgoingNormal(m_meshH, e_front_faces[1]) ;
		f1_normal.normalize();
		f1_normal = e.length()*f1_normal ;

		math::Point p_n0 = e_nodes[0].point();

		n0_2_id = m_FaceInfo[e_front_faces[0]].next_ideal_nodes[e_nodes[0].id()];
		Node n0_2_new = m_meshH->get<Node>(n0_2_id);
		math::Point p0_2_new = n0_2_new.point() ;
		math::Point p0_1_new = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, e_nodes[0].point(), dc, m_DistanceField, f0_normal.normalize());
		math::Point p0_3_new = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, e_nodes[0].point(), dc, m_DistanceField, f1_normal.normalize());

		Node n0_1_new = m_meshH->newNode(p0_1_new);
		Node n0_3_new = m_meshH->newNode(p0_3_new);

		n0_1_id = n0_1_new.id();
		n0_2_id = n0_2_new.id();
		n0_3_id = n0_3_new.id();

		//---------------------------//
		// Update the layer ids of   //
		// the new nodes				  //
		//---------------------------//
		var_node_couche_id->set(n0_1_id, AFront.getFrontID()+1);
		var_node_couche_id->set(n0_2_id, AFront.getFrontID()+1);
		var_node_couche_id->set(n0_3_id, AFront.getFrontID()+1);

		//---------------------------//
		// Update the edge corner	  //
		// information around the	  //
		//	node							  //
		//---------------------------//
		std::vector<TCellID> n0_adj_edges_on_front = AFront.orderedFrontEdgesAroundNode(m_meshH, n0_id);
		for (auto e_adj_id:n0_adj_edges_on_front)
		{
			if ( e_adj_id != e_id
			    && var_front_edges_classification->value(e_adj_id) == 1)
			{
				m_EdgeInfo[e_adj_id].CORNER_n_face_created[n0_id] = true;
				m_EdgeInfo[e_adj_id].diag_next_node[n0_id] = n0_2_id;

				NodeNeighbourhoodOnFront_3D n0_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n0_id);
				n0_neighbourhood.execute();

				TCellID next_edge_id_0side = n0_neighbourhood.nextEdgeOfFace(e_front_faces[0], e_id) ;
				TCellID next_edge_id_1side = n0_neighbourhood.nextEdgeOfFace(e_front_faces[1], e_id) ;

				TCellID f_adj_0 = n0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e_adj_id, next_edge_id_0side, next_edge_id_1side);
				TCellID f_adj_1 = n0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e_adj_id, next_edge_id_1side, next_edge_id_0side);

				//TCellID f_adj_0 = AFront.adjacentFaceOnFront(m_meshH, e_front_faces[0], e_adj_id);
				//TCellID f_adj_1 = AFront.adjacentFaceOnFront(m_meshH, e_front_faces[1], e_adj_id);

				std::pair<TCellID, TCellID> pair_0(f_adj_0, n0_id) ;
				std::pair<TCellID, TCellID> pair_1(f_adj_1, n0_id) ;

				m_EdgeInfo[e_adj_id].CORNER_next_nodes[pair_0] = n0_1_id ;
				m_EdgeInfo[e_adj_id].CORNER_next_nodes[pair_1] = n0_3_id ;


				m_FaceInfo[f_adj_0].next_nodes[n0_id] = n0_1_id ;
				m_FaceInfo[f_adj_1].next_nodes[n0_id] = n0_3_id ;
			}
		}


		//---------------------------//
		// Update the face 			  //
		// information					  //
		//---------------------------//
		m_FaceInfo[e_front_faces[0]].next_nodes[n0_id] = n0_1_id ;
		m_FaceInfo[e_front_faces[1]].next_nodes[n0_id] = n0_3_id ;

	}


	// Second side
	TCellID n1_id = e_nodes[1].id();
	TCellID n1_1_id(NullID);
	TCellID n1_2_id(NullID);
	TCellID n1_3_id(NullID);
	if (m_EdgeInfo[e_id].CORNER_n_face_created[n1_id])
	{
		// Then, a face has already been inserted on node e_nodes[1] with a node template
		std::pair<TCellID, TCellID> pair_1(e_front_faces[0], n1_id);
		n1_1_id = m_EdgeInfo[e_id].CORNER_next_nodes[pair_1] ;
		n1_2_id = m_EdgeInfo[e_id].diag_next_node[n1_id];
		std::pair<TCellID, TCellID> pair_2(e_front_faces[1], n1_id);
		n1_3_id = m_EdgeInfo[e_id].CORNER_next_nodes[pair_2] ;
	}
	else
	{
		//============================//
		// We create the 3 new nodes	//
		//============================//

		math::Vector3d f0_normal = AFront.outgoingNormal(m_meshH, e_front_faces[0]) ;
		f0_normal.normalize();
		f0_normal = e.length()*f0_normal ;

		math::Vector3d f1_normal = AFront.outgoingNormal(m_meshH, e_front_faces[1]) ;
		f1_normal.normalize();
		f1_normal = e.length()*f1_normal ;

		math::Point p_n1 = e_nodes[1].point();

		n1_2_id = m_FaceInfo[e_front_faces[0]].next_ideal_nodes[e_nodes[1].id()];
		Node n1_2_new = m_meshH->get<Node>(n1_2_id);
		math::Point p1_2_new = n1_2_new.point() ;
		math::Point p1_1_new = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, e_nodes[1].point(), dc, m_DistanceField, f0_normal.normalize());
		math::Point p1_3_new = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, e_nodes[1].point(), dc, m_DistanceField, f1_normal.normalize());

		Node n1_1_new = m_meshH->newNode(p1_1_new);
		Node n1_3_new = m_meshH->newNode(p1_3_new);

		n1_1_id = n1_1_new.id();
		n1_2_id = n1_2_new.id();
		n1_3_id = n1_3_new.id();

		//---------------------------//
		// Update the layer ids of   //
		// the new nodes				  //
		//---------------------------//
		var_node_couche_id->set(n1_1_id, AFront.getFrontID()+1);
		var_node_couche_id->set(n1_2_id, AFront.getFrontID()+1);
		var_node_couche_id->set(n1_3_id, AFront.getFrontID()+1);

		//---------------------------//
		// Update the edge corner	  //
		// information around the	  //
		//	node							  //
		//---------------------------//
		std::vector<TCellID> n1_adj_edges_on_front = AFront.orderedFrontEdgesAroundNode(m_meshH, n1_id);
		for (auto e_adj_id:n1_adj_edges_on_front)
		{
			if ( e_adj_id != e_id
			    && var_front_edges_classification->value(e_adj_id) == 1)
			{
				m_EdgeInfo[e_adj_id].CORNER_n_face_created[n1_id] = true;
				m_EdgeInfo[e_adj_id].diag_next_node[n1_id] = n1_2_id;

				NodeNeighbourhoodOnFront_3D n1_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n1_id);
				n1_neighbourhood.execute();

				TCellID next_edge_id_0side = n1_neighbourhood.nextEdgeOfFace(e_front_faces[0], e_id) ;
				TCellID next_edge_id_1side = n1_neighbourhood.nextEdgeOfFace(e_front_faces[1], e_id) ;

				TCellID f_adj_0 = n1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e_adj_id, next_edge_id_0side, next_edge_id_1side);
				TCellID f_adj_1 = n1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e_adj_id, next_edge_id_1side, next_edge_id_0side);

				//TCellID f_adj_0 = AFront.adjacentFaceOnFront(m_meshH, e_front_faces[0], e_adj_id);
				//TCellID f_adj_1 = AFront.adjacentFaceOnFront(m_meshH, e_front_faces[1], e_adj_id);

				std::pair<TCellID, TCellID> pair_0(f_adj_0, n1_id) ;
				std::pair<TCellID, TCellID> pair_1(f_adj_1, n1_id) ;

				m_EdgeInfo[e_adj_id].CORNER_next_nodes[pair_0] = n1_1_id ;
				m_EdgeInfo[e_adj_id].CORNER_next_nodes[pair_1] = n1_3_id ;


				m_FaceInfo[f_adj_0].next_nodes[n1_id] = n1_1_id ;
				m_FaceInfo[f_adj_1].next_nodes[n1_id] = n1_3_id ;
			}
		}


		//---------------------------//
		// Update the face 			  //
		// information					  //
		//---------------------------//
		m_FaceInfo[e_front_faces[0]].next_nodes[n1_id] = n1_1_id ;
		m_FaceInfo[e_front_faces[1]].next_nodes[n1_id] = n1_3_id ;

	}

	Node n0_1 = m_meshH->get<Node>(n0_1_id);
	Node n0_2 = m_meshH->get<Node>(n0_2_id);
	Node n0_3 = m_meshH->get<Node>(n0_3_id);
	Node n1_1 = m_meshH->get<Node>(n1_1_id);
	Node n1_2 = m_meshH->get<Node>(n1_2_id);
	Node n1_3 = m_meshH->get<Node>(n1_3_id);

	// Create the hexa
	r_id = math::Utils::CreateHexaNConnectivities(m_meshH, e_nodes[0], n0_1, n0_2, n0_3, e_nodes[1], n1_1, n1_2, n1_3);

	//---------------------------//
	// Update the two new faces  //
	// layer id                  //
	//---------------------------//
	Variable<int>* var_face_couche_id = m_meshH->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");
	TCellID f_1_new_id = math::Utils::CommonFace3Nodes(m_meshH, n0_1_id, n0_2_id, n1_2_id);
	TCellID f_2_new_id = math::Utils::CommonFace3Nodes(m_meshH, n1_2_id,  n0_2.id(), n0_3.id());
	var_face_couche_id->set(f_1_new_id, AFront.getFrontID()+1);
	var_face_couche_id->set(f_2_new_id, AFront.getFrontID()+1);

	gmds::IGMeshIOService ioService(m_meshH);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("AeroExtrusion_3D_"+std::to_string(m_iteration)+".vtk");
	m_iteration++;

	return r_id;
}
/*------------------------------------------------------------------------*/
TCellID
AeroExtrusion_3D::TemplateEdgeEnd(Front_3D &AFront, TCellID e_id, double dc, int mark_edgesTreated, int mark_facesTreated)
{
	//std::cout << "Template Edge End" << std::endl;
	TCellID r_id(NullID);

	Edge e = m_meshH->get<Edge>(e_id);
	std::vector<Node> e_nodes = e.get<Node>();

	std::vector<TCellID> e_front_faces = AFront.edgeFacesOnFront(m_meshH, e_id) ;

	Variable<int>* var_node_couche_id = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* var_front_edges_classification = m_meshH->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");

	TCellID n0_id = e_nodes[0].id();
	TCellID n1_id = e_nodes[1].id();

	NodeNeighbourhoodOnFront_3D n0_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n0_id);
	n0_neighbourhood.execute();
	NodeNeighbourhoodOnFront_3D n1_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n1_id);
	n1_neighbourhood.execute();


	// First side
	TCellID e_00_id = n0_neighbourhood.nextEdgeOfFace(e_front_faces[0], e_id);
	Edge e_00 = m_meshH->get<Edge>(e_00_id);
	Node n_00 = e_00.getOppositeNode(m_meshH->get<Node>(n0_id));

	TCellID e_01_id = n0_neighbourhood.nextEdgeOfFace(e_front_faces[1], e_id);
	Edge e_01 = m_meshH->get<Edge>(e_01_id);
	Node n_01 = e_01.getOppositeNode(m_meshH->get<Node>(n0_id));

	Node n_02;
	if (!m_EdgeInfo[e_id].END_n_face_created[n0_id])
	{
		n_02 = m_meshH->get<Node>(m_FaceInfo[e_front_faces[0]].next_ideal_nodes[n0_id]);
		m_EdgeInfo[e_id].END_n_face_created[n0_id] = true;
		m_EdgeInfo[e_id].diag_next_node[n0_id] = n_02.id();
	}
	else
	{
		n_02 = m_meshH->get<Node>(m_EdgeInfo[e_id].diag_next_node[n0_id]);
	}

	// Second side
	TCellID e_10_id = n1_neighbourhood.nextEdgeOfFace(e_front_faces[0], e_id);
	Edge e_10 = m_meshH->get<Edge>(e_10_id);
	Node n_10 = e_10.getOppositeNode(m_meshH->get<Node>(n1_id));

	TCellID e_11_id = n1_neighbourhood.nextEdgeOfFace(e_front_faces[1], e_id);
	Edge e_11 = m_meshH->get<Edge>(e_11_id);
	Node n_11 = e_11.getOppositeNode(m_meshH->get<Node>(n1_id));

	Node n_12;
	if (!m_EdgeInfo[e_id].END_n_face_created[n1_id])
	{
		n_12 = m_meshH->get<Node>(m_FaceInfo[e_front_faces[0]].next_ideal_nodes[n1_id]);
		m_EdgeInfo[e_id].END_n_face_created[n1_id] = true;
		m_EdgeInfo[e_id].diag_next_node[n1_id] = n_12.id();
	}
	else
	{
		n_12 = m_meshH->get<Node>(m_EdgeInfo[e_id].diag_next_node[n1_id]);
	}

	// Create the hexa
	r_id = math::Utils::CreateHexaNConnectivities(m_meshH, e_nodes[0], n_00, n_02, n_01, e_nodes[1], n_10, n_12, n_11);

	//---------------------------//
	// Update the two new nodes  //
	// layer id                  //
	//---------------------------//
	var_node_couche_id->set(n_02.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n_12.id(), AFront.getFrontID()+1);

	// Update the FaceInfo structures
	Edge e_opp_0 = math::Utils::oppositeEdgeInFace(m_meshH, e_id, e_front_faces[0]);
	TCellID f_opp_0 = AFront.adjacentFaceOnFront(m_meshH, e_front_faces[0], e_opp_0.id());
	m_FaceInfo[f_opp_0].next_nodes[n_00.id()] = n_02.id() ;
	m_FaceInfo[f_opp_0].next_nodes[n_10.id()] = n_12.id() ;

	Edge e_opp_1 = math::Utils::oppositeEdgeInFace(m_meshH, e_id, e_front_faces[1]);
	TCellID f_opp_1 = AFront.adjacentFaceOnFront(m_meshH, e_front_faces[1], e_opp_1.id());
	m_FaceInfo[f_opp_1].next_nodes[n_01.id()] = n_02.id() ;
	m_FaceInfo[f_opp_1].next_nodes[n_11.id()] = n_12.id() ;

	// Mark the faces of the front as treated
	m_meshH->mark(m_meshH->get<Face>(e_front_faces[0]), mark_facesTreated);
	m_meshH->mark(m_meshH->get<Face>(e_front_faces[1]), mark_facesTreated);

	gmds::IGMeshIOService ioService(m_meshH);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("AeroExtrusion_3D_"+std::to_string(m_iteration)+".vtk");
	m_iteration++;

	return r_id;
}
/*------------------------------------------------------------------------*/
void
AeroExtrusion_3D::TemplateFace(TCellID f_id, Front_3D &Front_IN, std::map<TCellID, TCellID> map_new_nodes)
{
	//std::cout << "Template Face" << std::endl;

	Variable<int>* var_node_couche_id = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* var_face_couche_id = m_meshH->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");

	Face f = m_meshH->get<Face>(f_id);
	std::vector<Node> nodes = f.get<Node>();

	TCellID n0_id = m_FaceInfo[f_id].next_nodes[nodes[0].id()] ;
	TCellID n1_id = m_FaceInfo[f_id].next_nodes[nodes[1].id()] ;
	TCellID n2_id = m_FaceInfo[f_id].next_nodes[nodes[2].id()] ;
	TCellID n3_id = m_FaceInfo[f_id].next_nodes[nodes[3].id()] ;


	Node n0 = m_meshH->get<Node>(n0_id);
	Node n1 = m_meshH->get<Node>(n1_id);
	Node n2 = m_meshH->get<Node>(n2_id);
	Node n3 = m_meshH->get<Node>(n3_id);
	var_node_couche_id->set(n0_id, Front_IN.getFrontID()+1);
	var_node_couche_id->set(n1_id, Front_IN.getFrontID()+1);
	var_node_couche_id->set(n2_id, Front_IN.getFrontID()+1);
	var_node_couche_id->set(n3_id, Front_IN.getFrontID()+1);

	TCellID r_id = math::Utils::CreateHexaNConnectivities(m_meshH, nodes[0], nodes[1], nodes[2], nodes[3], n0, n1, n2, n3);

	TCellID f_new_layer_id = math::Utils::GetOrCreateQuadAndConnectivities(m_meshH, n0.id(), n1.id(), n2.id(), n3.id());
	var_face_couche_id->set(f_new_layer_id, Front_IN.getFrontID()+1);

	gmds::IGMeshIOService ioService(m_meshH);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("AeroExtrusion_3D_"+std::to_string(m_iteration)+".vtk");
	m_iteration++;

}
/*------------------------------------------------------------------------*/