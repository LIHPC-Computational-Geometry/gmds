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

AeroExtrusion_3D::AeroExtrusion_3D(Mesh *AMeshT, Mesh *AMeshH, ParamsAero& Aparams_aero, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField) :
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
	m_Patterns = m_meshH->getOrCreateVariable<int, GMDS_REGION>("GMDS_Patterns");
}


/*------------------------------------------------------------------------*/
AeroExtrusion_3D::STATUS
AeroExtrusion_3D::execute()
{
	//double pas_couche = 1.0/m_params_aero.nbr_couches ;

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

	double max_distance = ComputeMaxDistOnFront(Current_Front, m_meshT->getVariable<double,GMDS_NODE>("GMDS_Distance"));
	double pas_couche = (1.0-max_distance)/(m_params_aero.nbr_couches-1) ;

	// Compute the successive layers
	for (int i=2; i <= m_params_aero.nbr_couches; i++) {
		Current_Front = ComputeLayer(Current_Front, m_meshT->getVariable<double,GMDS_NODE>("GMDS_Distance"), max_distance+(i-1)*pas_couche,
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

	TInt mark_edgesTreated = m_meshH->newMark<Edge>();
	TInt mark_facesTreated = m_meshH->newMark<Face>();

	TInt mark_EdgesTemplates = Classification.getMarkEdgesTemplates();
	TInt mark_NodesTemplates = Classification.getMarkNodesTemplates();

	Variable<int>* var_TEST = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_TEST_NODES");

	// Create the HEXA on each NODE classified
	std::map<TCellID, int> singular_nodes = getSingularNodes(Front_IN, var_front_edges_classification);
	for (auto singu_node:singular_nodes)
	{
		TCellID n_id = singu_node.first ;
		int singu_type = singu_node.second;
		if (singu_type==1 && m_meshH->isMarked(m_meshH->get<Node>(n_id),mark_NodesTemplates))
		{
			TCellID r_id = TemplateNode3Corner(Front_IN, n_id, map_new_nodes, dist_cible);
			m_Patterns->set(r_id, 1);
		}
		else if (singu_type==2 && m_meshH->isMarked(m_meshH->get<Node>(n_id),mark_NodesTemplates))
		{
			TCellID r_id = TemplateNode2Corner1End(Front_IN, n_id, dist_cible, mark_edgesTreated, mark_facesTreated);
			m_Patterns->set(r_id, 2);
		}
		else if (singu_type==3 && m_meshH->isMarked(m_meshH->get<Node>(n_id),mark_NodesTemplates))
		{
			TCellID r_id = TemplateNode1Corner2End(Front_IN, n_id, dist_cible, mark_edgesTreated, mark_facesTreated);
			m_Patterns->set(r_id, 3);
		}
		else if (singu_type==4 && m_meshH->isMarked(m_meshH->get<Node>(n_id),mark_NodesTemplates))
		{
			TCellID r_id = TemplateNode3End(Front_IN, n_id, mark_edgesTreated, mark_facesTreated);
			m_Patterns->set(r_id, 4);
		}
		else if (singu_type==5 && m_meshH->isMarked(m_meshH->get<Node>(n_id),mark_NodesTemplates))
		{
			TCellID r_id = TemplateNode3Corner3End(Front_IN, n_id, dist_cible, mark_edgesTreated, mark_facesTreated);
			m_Patterns->set(r_id, 5);
		}
		else if (singu_type==6 && m_meshH->isMarked(m_meshH->get<Node>(n_id),mark_NodesTemplates))
		{
			std::vector<TCellID> r_id = TemplateNode2Corner2End(Front_IN, n_id, dist_cible, mark_edgesTreated, mark_facesTreated); // Create 2 hexas
			m_Patterns->set(r_id[0], 6);
			m_Patterns->set(r_id[1], 6);
		}
		else if (singu_type==7 && m_meshH->isMarked(m_meshH->get<Node>(n_id),mark_NodesTemplates))
		{
			std::vector<TCellID> r_id = TemplateNode2Corner1Reversal(Front_IN, n_id, dist_cible); // Create 2 hexas
			m_Patterns->set(r_id[0], 7);
			m_Patterns->set(r_id[1], 7);
		}
		else if (singu_type==8 && m_meshH->isMarked(m_meshH->get<Node>(n_id),mark_NodesTemplates))
		{
			std::vector<TCellID> r_id = TemplateNode2End1Reversal(Front_IN, n_id, dist_cible, mark_edgesTreated, mark_facesTreated); // Create 2 hexas
			m_Patterns->set(r_id[0], 8);
			m_Patterns->set(r_id[1], 8);
		}
		if (m_meshH->isMarked(m_meshH->get<Node>(n_id),mark_NodesTemplates))
		{
			//std::cout << "Noeud marqué pour template" << std::endl;
			var_TEST->set(n_id, 1);
		}
	}


	// Create the HEXA on each EDGE classified
	std::map<TCellID, int> singular_edges = getSingularEdges(Front_IN, var_front_edges_classification, mark_edgesTreated);
	for (auto singu_edge:singular_edges)
	{
		TCellID e_id = singu_edge.first ;
		int singu_type = singu_edge.second;
		if (singu_type==1
		    && !m_meshH->isMarked(m_meshH->get<Edge>(e_id), mark_edgesTreated)
		    && m_meshH->isMarked(m_meshH->get<Edge>(e_id), mark_EdgesTemplates))
		{
			TCellID r_id = TemplateEdgeCorner(Front_IN, e_id, dist_cible);
			m_Patterns->set(r_id, -1);
		}
		else if (singu_type==2
		         && !m_meshH->isMarked(m_meshH->get<Edge>(e_id), mark_edgesTreated)
		         && m_meshH->isMarked(m_meshH->get<Edge>(e_id), mark_EdgesTemplates))
		{
			TCellID r_id = TemplateEdgeEnd(Front_IN, e_id, dist_cible, mark_edgesTreated, mark_facesTreated);
			m_Patterns->set(r_id, -2);
		}
		else if (singu_type==3
		         && !m_meshH->isMarked(m_meshH->get<Edge>(e_id), mark_edgesTreated)
		         && m_meshH->isMarked(m_meshH->get<Edge>(e_id), mark_EdgesTemplates))
		{
			std::vector<TCellID> r_id = TemplateEdgeReversal(Front_IN, e_id, dist_cible, mark_edgesTreated, mark_facesTreated);
			m_Patterns->set(r_id[0], -3);
			m_Patterns->set(r_id[1], -3);
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

	for (auto n_id:m_meshH->nodes())
	{
		if (var_NODE_couche_id->value(n_id) == Front_OUT.getFrontID())
		{
			gmds::Cell::Data data = m_fl.find(m_meshH->get<Node>(n_id).point());
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
			exit(1);
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

	for (auto n_id:m_meshH->nodes())
	{
		if (var_NODE_couche_id->value(n_id)==1)
		{
			NodeNeighbourhoodOnFront_3D n_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &Front_OUT, n_id);
			std::vector<TCellID> n_edges = n_neighbourhood.getOrderedEdges() ;
			math::Point p({0,0,0});

		}
	}

	/*
	for (auto n_id:Front_OUT.getNodes())
	{
		Node n = m_meshH->get<Node>(n_id);
		TCellID r_id = m_fl.findTetra(n.point());
		Region r = m_meshT->get<Region>(r_id);
		std::vector<Node> r_nodes = r.get<Node>();
		double dist = math::Utils::linearInterpolation2D3Pt(r_nodes[0].point(), r_nodes[1].point(), r_nodes[2].point(), r_nodes[3].point(), n.point(),
		                                      m_DistanceField->value(r_nodes[0].id()), m_DistanceField->value(r_nodes[1].id()),
		                                      m_DistanceField->value(r_nodes[2].id()), m_DistanceField->value(r_nodes[3].id()));
	}
	 */


	return Front_OUT;

}
/*------------------------------------------------------------------------*/
double
AeroExtrusion_3D::ComputeMaxDistOnFront(Front_3D Front_IN, Variable<double>* A_distance)
{
	double max_distance(0.0);

	for (auto n_id:Front_IN.getNodes())
	{
		Node n = m_meshH->get<Node>(n_id) ;

		gmds::Cell::Data data = m_fl.find(n.point());
		TCellID n_closest_id = data.id;
		Node n_closest = m_meshT->get<Node>(n_closest_id);
		if (A_distance->value(n_closest_id) > max_distance)
		{
			max_distance = A_distance->value(n_closest_id);
		}
	}

	return max_distance;
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

	Front_3D Front_OUT = Front_3D(Front_IN.getFrontID()+1, new_front_nodes_id, new_front_faces_id);
	return Front_OUT;

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
		else if (compteur_corner==1 && compteur_end==2)
		{
			sing_nodes[n_id] = 3;
		}
		else if (compteur_end == 3 && n_ordered_edges.size()==3)
		{
			sing_nodes[n_id] = 4;
		}
		else if (n_ordered_edges.size() == 6
		         && ( ( front_edges_classification->value(n_ordered_edges[0]) == 1
		             && front_edges_classification->value(n_ordered_edges[2]) == 1
		             && front_edges_classification->value(n_ordered_edges[4]) == 1
		             && front_edges_classification->value(n_ordered_edges[1]) == 2
		             && front_edges_classification->value(n_ordered_edges[3]) == 2
		             && front_edges_classification->value(n_ordered_edges[5]) == 2 )
		         	|| ( front_edges_classification->value(n_ordered_edges[1]) == 1
		                 && front_edges_classification->value(n_ordered_edges[3]) == 1
		                 && front_edges_classification->value(n_ordered_edges[5]) == 1
		                 && front_edges_classification->value(n_ordered_edges[0]) == 2
		                 && front_edges_classification->value(n_ordered_edges[2]) == 2
		                 && front_edges_classification->value(n_ordered_edges[4]) == 2 ) ) )
		{
			sing_nodes[n_id] = 5;
		}
		else if (n_ordered_edges.size() == 6
		         && compteur_corner == 2 && compteur_end==2 && compteur_reversal==0)
		{
			sing_nodes[n_id] = 6;
		}
		else if (n_ordered_edges.size() == 3
		         && compteur_corner == 2 && compteur_end==0 && compteur_reversal==1)
		{
			sing_nodes[n_id] = 7;
		}
		else if (n_ordered_edges.size() == 6
		         && compteur_corner == 0 && compteur_end==2 && compteur_reversal==1)
		{
			sing_nodes[n_id] = 8;
		}

	}
	return sing_nodes;
}
/*------------------------------------------------------------------------*/
std::map<TCellID, int>
AeroExtrusion_3D::getSingularEdges(Front_3D &AFront, Variable<int>* front_edges_classification, TInt mark_singEdgesTreated)
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

	Edge ec_0 = m_meshH->get<Edge>(n_edges_corner[0]);
	Edge ec_1 = m_meshH->get<Edge>(n_edges_corner[1]);
	Edge ec_2 = m_meshH->get<Edge>(n_edges_corner[2]);

	Node n = m_meshH->get<Node>(n_id);
	Node n6 = m_meshH->get<Node>(map_new_nodes[n_id]);

	Node n_c0 = ec_0.getOppositeNode(n);
	Node n_c1 = ec_1.getOppositeNode(n);
	Node n_c2 = ec_2.getOppositeNode(n);

	/*
	math::Vector3d v1 = n_c0.point()-n.point() ;
	math::Vector3d v2 = n_c1.point()-n.point() ;
	math::Vector3d v3 = n_c2.point()-n.point() ;

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
	*/
	math::Vector3d v1 = AFront.outgoingNormal(m_meshH, n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec_1.id(), ec_2.id(), ec_0.id())) ;
	math::Vector3d v2 = AFront.outgoingNormal(m_meshH, n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec_0.id(), ec_2.id(), ec_1.id())) ;
	math::Vector3d v3 = AFront.outgoingNormal(m_meshH, n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec_0.id(), ec_1.id(), ec_2.id())) ;
	v1.normalize();
	v2.normalize();
	v3.normalize();
	math::Point p1 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, v1);
	math::Point p3 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, v2);
	math::Point p4 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, v3);
	math::Point p2 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, v1+v2);
	math::Point p5 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, v1+v3);
	math::Point p7 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, v2+v3);

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
	std::vector<TCellID> adj_faces = n_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec_0.id(), ec_1.id(), ec_2.id());
	for (auto f_adj_id:adj_faces)
	{
		m_FaceInfo[f_adj_id].next_nodes[n_id] = n4.id() ;
	}
	adj_faces = n_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec_1.id(), ec_2.id(), ec_0.id());
	for (auto f_adj_id:adj_faces)
	{
		m_FaceInfo[f_adj_id].next_nodes[n_id] = n1.id() ;
	}
	adj_faces = n_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec_0.id(), ec_2.id(), ec_1.id());
	for (auto f_adj_id:adj_faces)
	{
		m_FaceInfo[f_adj_id].next_nodes[n_id] = n3.id() ;
	}

	TCellID f_1 = math::Utils::CommonFace3Nodes(m_meshH, n_id, n_c0.id(), n_c1.id());
	TCellID f_2 = math::Utils::CommonFace3Nodes(m_meshH, n_id, n_c1.id(), n_c2.id());
	TCellID f_3 = math::Utils::CommonFace3Nodes(m_meshH, n_id, n_c2.id(), n_c0.id());
	/*
	m_FaceInfo[f_1].next_nodes[n_id] = n4.id() ;
	m_FaceInfo[f_2].next_nodes[n_id] = n1.id() ;
	m_FaceInfo[f_3].next_nodes[n_id] = n3.id() ;
	 */

	//---------------------------//
	// Update the corner edges   //
	// information					  //
	//---------------------------//
	TCellID e1_id = math::Utils::CommonEdge(m_meshH, n.id(), n_c0.id());
	TCellID e2_id = math::Utils::CommonEdge(m_meshH, n.id(), n_c1.id());
	TCellID e3_id = math::Utils::CommonEdge(m_meshH, n.id(), n_c2.id());

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
AeroExtrusion_3D::TemplateNode2Corner1End(Front_3D &AFront, TCellID n_id, double dc, TInt mark_edgesTreated, TInt mark_facesTreated)
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
AeroExtrusion_3D::TemplateNode1Corner2End(Front_3D &AFront, TCellID n_id, double dc, TInt mark_edgesTreated, TInt mark_facesTreated)
{
	//std::cout << "Template Node 1 Corner 2 End au noeud " << n_id << std::endl;
	TCellID r_id;

	Node n = m_meshH->get<Node>(n_id);

	NodeNeighbourhoodOnFront_3D n_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n_id);
	n_neighbourhood.execute();

	Variable<int>* var_front_edges_classification = m_meshH->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");

	std::vector<Edge> ee;	// For the two edge END
	Edge ec;						// For the edge CORNER

	// Get the local feature edges
	for (auto edge_id:n_neighbourhood.getOrderedEdges())
	{
		if (var_front_edges_classification->value(edge_id)==2)
		{
			ee.push_back(m_meshH->get<Edge>(edge_id));
		}
		else if (var_front_edges_classification->value(edge_id)==1)
		{
			ec = m_meshH->get<Edge>(edge_id);
		}
	}

	Node n_e0 = ee[0].getOppositeNode(n);
	Node n_e1 = ee[1].getOppositeNode(n);

	std::vector<TCellID> ee_faces = n_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ee[0].id(), ee[1].id(), ec.id()) ;
	if (ee_faces.size() != 3)
	{
		std::cout << "Attention AeroExtrusion_3D: Template Node 1 CORNER, 2 END, ne respecte pas la condition des 3 faces entre les deux Edges de type END." << std::endl;
	}

	Face f1 = m_meshH->get<Face>(ee_faces[1]);

	Edge e1 = m_meshH->get<Edge>(n_neighbourhood.nextEdgeOfFace(ee_faces[0], ee[0].id()));
	Edge e3 = m_meshH->get<Edge>(n_neighbourhood.nextEdgeOfFace(ee_faces[2], ee[1].id()));
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

	Node n4 = m_meshH->get<Node>(m_FaceInfo[f1.id()].next_ideal_nodes[n1.id()]);
	Node n5 = m_meshH->get<Node>(m_FaceInfo[f1.id()].next_ideal_nodes[n2.id()]);
	Node n6 = m_meshH->get<Node>(m_FaceInfo[f1.id()].next_ideal_nodes[n3.id()]);

	// Create the new hexa
	r_id = math::Utils::CreateHexaNConnectivities(m_meshH, n, n1, n2, n3, nc, n4, n5, n6);

	// Update the faces and edges that can't do actions anymore
	m_meshH->mark(ec, mark_edgesTreated);
	m_meshH->mark(f1, mark_facesTreated);

	// Update the EdgeInfo structure for the END edges ---->
	m_EdgeInfo[ee[0].id()].END_n_face_created[n_id] = true;
	m_EdgeInfo[ee[0].id()].diag_next_node[n_id] = n4.id();
	m_EdgeInfo[ee[1].id()].END_n_face_created[n_id] = true;
	m_EdgeInfo[ee[1].id()].diag_next_node[n_id] = n6.id();
	//<----

	// Update the EdgeInfo structure for the CORNER edge ---->
	NodeNeighbourhoodOnFront_3D nc_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, nc.id());
	nc_neighbourhood.execute();

	std::vector<TCellID> nc_front_edges = nc_neighbourhood.getOrderedEdges();
	Edge ec_next;
	for (auto e_loc_id:nc_front_edges)
	{
		if (var_front_edges_classification->value(e_loc_id) == 1
		    && e_loc_id != ec.id())
		{
			ec_next = m_meshH->get<Edge>(e_loc_id) ;
		}
	}
	m_EdgeInfo[ec_next.id()].CORNER_n_face_created[nc.id()] = true;
	m_EdgeInfo[ec_next.id()].diag_next_node[nc.id()] = n5.id();

	Face f_ce_0 = m_meshH->get<Face>( n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec.id(), ee[0].id(), ee[1].id()) );
	Face f_ce_1 = m_meshH->get<Face>( n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec.id(), ee[1].id(), ee[0].id()) );

	Edge e_e0_opp = math::Utils::oppositeEdgeInFace(m_meshH, ee[0].id(), f_ce_0.id()) ;
	Edge e_e1_opp = math::Utils::oppositeEdgeInFace(m_meshH, ee[1].id(), f_ce_1.id()) ;

	Face f_ec_next_0 = m_meshH->get<Face>(nc_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec_next.id(), e_e0_opp.id(), e_e1_opp.id()) );
	Face f_ec_next_1 = m_meshH->get<Face>(nc_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec_next.id(), e_e1_opp.id(), e_e0_opp.id()) );

	std::pair<TCellID, TCellID> pair_1(f_ec_next_0.id(), nc.id());
	std::pair<TCellID, TCellID> pair_2(f_ec_next_1.id(), nc.id());

	m_EdgeInfo[ec_next.id()].CORNER_next_nodes[pair_1] = n4.id() ;
	m_EdgeInfo[ec_next.id()].CORNER_next_nodes[pair_2] = n6.id() ;
	//<----

	// Update the FaceInfo structures ---->
	for (auto f_loc_id:nc_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec_next.id(), ec.id(), e_e1_opp.id()))
	{
		m_FaceInfo[f_loc_id].next_nodes[nc.id()] = n4.id();
	}
	for (auto f_loc_id:nc_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec_next.id(), ec.id(), e_e0_opp.id()))
	{
		m_FaceInfo[f_loc_id].next_nodes[nc.id()] = n6.id();
	}
	//<----

	// Write the new hexa in a VTK file <----
	gmds::IGMeshIOService ioService(m_meshH);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("AeroExtrusion_3D_"+std::to_string(m_iteration)+".vtk");
	m_iteration++;
	//<----

	return r_id;
}
/*------------------------------------------------------------------------*/
TCellID
AeroExtrusion_3D::TemplateNode3End(Front_3D &AFront, TCellID n_id, TInt mark_edgesTreated, TInt mark_facesTreated)
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
	for (auto const& n_loc:f0_nodes)
	{
		if (n_loc.id() != n_id && n_loc.id() != n_e0_opp.id() && n_loc.id() != n_e1_opp.id())
		{
			n_f0_diag = n_loc;
		}
	}
	for (auto const& n_loc:f1_nodes)
	{
		if (n_loc.id() != n_id && n_loc.id() != n_e1_opp.id() && n_loc.id() != n_e2_opp.id())
		{
			n_f1_diag = n_loc;
		}
	}
	for (auto const& n_loc:f2_nodes)
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
AeroExtrusion_3D::TemplateNode3Corner3End(Front_3D &AFront, TCellID n_id, double dc, TInt mark_edgesTreated, TInt mark_facesTreated)
{
	//std::cout << "Template Node 3 Corner 3 End au noeud " << n_id << std::endl;
	TCellID r_id;

	Node n = m_meshH->get<Node>(n_id);

	NodeNeighbourhoodOnFront_3D n_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n_id);
	n_neighbourhood.execute();
	std::vector<TCellID> n_ordered_edges = n_neighbourhood.getOrderedEdges();

	Variable<int>* var_front_edges_classification = m_meshH->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
	std::vector<Edge> ec;
	std::vector<Edge> ee;
	if (var_front_edges_classification->value(n_ordered_edges[0]) == 1)
	{
		ec.push_back(m_meshH->get<Edge>(n_ordered_edges[0]));
		ec.push_back(m_meshH->get<Edge>(n_ordered_edges[2]));
		ec.push_back(m_meshH->get<Edge>(n_ordered_edges[4]));
		ee.push_back(m_meshH->get<Edge>(n_ordered_edges[1]));
		ee.push_back(m_meshH->get<Edge>(n_ordered_edges[3]));
		ee.push_back(m_meshH->get<Edge>(n_ordered_edges[5]));
	}
	else if (var_front_edges_classification->value(n_ordered_edges[0])==2)
	{
		ee.push_back(m_meshH->get<Edge>(n_ordered_edges[0]));
		ee.push_back(m_meshH->get<Edge>(n_ordered_edges[2]));
		ee.push_back(m_meshH->get<Edge>(n_ordered_edges[4]));
		ec.push_back(m_meshH->get<Edge>(n_ordered_edges[5]));
		ec.push_back(m_meshH->get<Edge>(n_ordered_edges[1]));
		ec.push_back(m_meshH->get<Edge>(n_ordered_edges[3]));
	}

	Node nc0 = ec[0].getOppositeNode(n);
	Node nc1 = ec[1].getOppositeNode(n);
	Node nc2 = ec[2].getOppositeNode(n);

	Node n1 = m_meshH->get<Node>(m_FaceInfo[n.get<Face>()[0].id()].next_ideal_nodes[n_id]);

	math::Vector3d v01 = ( (nc0.point() - n.point()).normalize() + (nc1.point() - n.point()).normalize() + (n1.point() - n.point()).normalize() ).normalize() ;
	math::Vector3d v12 = ( (nc1.point() - n.point()).normalize() + (nc2.point() - n.point()).normalize() + (n1.point() - n.point()).normalize() ).normalize() ;
	math::Vector3d v02 = ( (nc0.point() - n.point()).normalize() + (nc2.point() - n.point()).normalize() + (n1.point() - n.point()).normalize() ).normalize() ;

	math::Point pe0 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, v01);
	math::Point pe1 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, v12);
	math::Point pe2 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, v02);

	Node ne0 = m_meshH->newNode(pe0);
	Node ne1 = m_meshH->newNode(pe1);
	Node ne2 = m_meshH->newNode(pe2);

	// Create the new hexa
	r_id = math::Utils::CreateHexaNConnectivities(m_meshH, n, nc0, ne0, nc1, nc2, ne2, n1, ne1);

	// Update the layer ids ---->
	Variable<int>* var_node_couche_id = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	var_node_couche_id->set(ne0.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(ne1.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(ne2.id(), AFront.getFrontID()+1);
	// <----

	// Update the edge struc info for the 3 end edges ---->
	m_EdgeInfo[ee[0].id()].END_n_face_created[n_id] = true;
	m_EdgeInfo[ee[1].id()].END_n_face_created[n_id] = true;
	m_EdgeInfo[ee[2].id()].END_n_face_created[n_id] = true;
	m_EdgeInfo[ee[0].id()].diag_next_node[n_id] = ne0.id();
	m_EdgeInfo[ee[1].id()].diag_next_node[n_id] = ne1.id();
	m_EdgeInfo[ee[2].id()].diag_next_node[n_id] = ne2.id();
	// <----

	TCellID f0_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec[0].id(), ee[0].id(), ee[2].id());
	TCellID f1_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ee[0].id(), ec[1].id(), ec[0].id());
	TCellID f2_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec[1].id(), ee[1].id(), ee[0].id());
	TCellID f3_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ee[1].id(), ec[2].id(), ec[1].id());
	TCellID f4_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec[2].id(), ee[2].id(), ee[1].id());
	TCellID f5_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ee[2].id(), ec[0].id(), ec[2].id());

	// Update the edge struc info for the first corner edge ---->
	NodeNeighbourhoodOnFront_3D n_c0_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, nc0.id()) ;
	n_c0_neighbourhood.execute();
	for (auto edge_id:n_c0_neighbourhood.getOrderedEdges()) {
		if (var_front_edges_classification->value(edge_id) == 1 && edge_id != ec[0].id()) {
			m_EdgeInfo[edge_id].CORNER_n_face_created[nc0.id()] = true;
			m_EdgeInfo[edge_id].diag_next_node[nc0.id()] = n1.id();

			TCellID edge_opp_f0_e0 = n_c0_neighbourhood.nextEdgeOfFace(f0_id, ec[0].id());
			for (auto faces_id : n_c0_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[0].id(), edge_id, edge_opp_f0_e0)) {
				m_FaceInfo[faces_id].next_nodes[nc0.id()] = ne2.id();
			}
			TCellID edge_opp_f5_ee2 = n_c0_neighbourhood.nextEdgeOfFace(f5_id, ec[0].id());
			for (auto faces_id : n_c0_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[0].id(), edge_id, edge_opp_f5_ee2)) {
				m_FaceInfo[faces_id].next_nodes[nc0.id()] = ne0.id();
			}

			std::pair<TCellID, TCellID> pair_1(n_c0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[0].id(), edge_opp_f5_ee2), nc0.id());
			std::pair<TCellID, TCellID> pair_2(n_c0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[0].id(), edge_opp_f0_e0), nc0.id());

			m_EdgeInfo[edge_id].CORNER_next_nodes[pair_1] = ne0.id();
			m_EdgeInfo[edge_id].CORNER_next_nodes[pair_2] = ne2.id();
		}
	}
	// <----

	// Update the edge struc info for the second corner edge ---->
	NodeNeighbourhoodOnFront_3D n_c1_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, nc1.id()) ;
	n_c1_neighbourhood.execute();
	for (auto edge_id:n_c1_neighbourhood.getOrderedEdges()) {
		if (var_front_edges_classification->value(edge_id) == 1 && edge_id != ec[1].id()) {
			m_EdgeInfo[edge_id].CORNER_n_face_created[nc1.id()] = true;
			m_EdgeInfo[edge_id].diag_next_node[nc1.id()] = n1.id();

			TCellID edge_opp_f1_e0 = n_c1_neighbourhood.nextEdgeOfFace(f1_id, ec[1].id());
			for (auto faces_id : n_c1_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[1].id(), edge_id, edge_opp_f1_e0)) {
				m_FaceInfo[faces_id].next_nodes[nc1.id()] = ne0.id();
			}
			TCellID edge_opp_f2_ee1 = n_c1_neighbourhood.nextEdgeOfFace(f2_id, ec[1].id());
			for (auto faces_id : n_c1_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[1].id(), edge_id, edge_opp_f2_ee1)) {
				m_FaceInfo[faces_id].next_nodes[nc0.id()] = ne0.id();
			}

			std::pair<TCellID, TCellID> pair_1(n_c1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[1].id(), edge_opp_f2_ee1), nc1.id());
			std::pair<TCellID, TCellID> pair_2(n_c1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[1].id(), edge_opp_f1_e0), nc1.id());

			m_EdgeInfo[edge_id].CORNER_next_nodes[pair_1] = ne0.id();
			m_EdgeInfo[edge_id].CORNER_next_nodes[pair_2] = ne1.id();
		}
	}
	// <----

	// Update the edge struc info for the thirst corner edge ---->
	NodeNeighbourhoodOnFront_3D n_c2_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, nc2.id()) ;
	n_c2_neighbourhood.execute();
	for (auto edge_id:n_c2_neighbourhood.getOrderedEdges()) {
		if (var_front_edges_classification->value(edge_id) == 1 && edge_id != ec[2].id()) {
			m_EdgeInfo[edge_id].CORNER_n_face_created[nc2.id()] = true;
			m_EdgeInfo[edge_id].diag_next_node[nc2.id()] = n1.id();

			TCellID edge_opp_f3_ee1 = n_c2_neighbourhood.nextEdgeOfFace(f3_id, ec[2].id());
			for (auto faces_id : n_c2_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[2].id(), edge_id, edge_opp_f3_ee1)) {
				m_FaceInfo[faces_id].next_nodes[nc2.id()] = ne2.id();
			}
			TCellID edge_opp_f4_ee2 = n_c2_neighbourhood.nextEdgeOfFace(f4_id, ec[2].id());
			for (auto faces_id : n_c2_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[2].id(), edge_id, edge_opp_f4_ee2)) {
				m_FaceInfo[faces_id].next_nodes[nc2.id()] = ne1.id();
			}

			std::pair<TCellID, TCellID> pair_1(n_c2_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[2].id(), edge_opp_f4_ee2), nc2.id());
			std::pair<TCellID, TCellID> pair_2(n_c2_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[2].id(), edge_opp_f3_ee1), nc2.id());

			m_EdgeInfo[edge_id].CORNER_next_nodes[pair_1] = ne1.id();
			m_EdgeInfo[edge_id].CORNER_next_nodes[pair_2] = ne2.id();
		}
	}
	// <----

	// Update the faces and edges that can't do actions anymore
	m_meshH->mark(ec[0], mark_edgesTreated);
	m_meshH->mark(ec[1], mark_edgesTreated);
	m_meshH->mark(ec[2], mark_edgesTreated);

	// Write the new hexa in a VTK file <----
	gmds::IGMeshIOService ioService(m_meshH);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("AeroExtrusion_3D_"+std::to_string(m_iteration)+".vtk");
	m_iteration++;
	//<----

	return r_id;
}
/*------------------------------------------------------------------------*/
std::vector<TCellID>
AeroExtrusion_3D::TemplateNode2Corner2End(Front_3D &AFront, TCellID n_id, double dc, TInt mark_edgesTreated, TInt mark_facesTreated)
{
	//std::cout << "Template Node 2 Corner 2 End au noeud " << n_id << std::endl;
	std::vector<TCellID> hexas_id;

	Node n = m_meshH->get<Node>(n_id);

	NodeNeighbourhoodOnFront_3D n_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n_id);
	n_neighbourhood.execute();
	std::vector<TCellID> n_ordered_edges = n_neighbourhood.getOrderedEdges();

	Variable<int>* var_front_edges_classification = m_meshH->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
	std::vector<Edge> ec;
	std::vector<Edge> ee;
	for (auto e_loc_id:n_ordered_edges)
	{
		if (var_front_edges_classification->value(e_loc_id)==1)
		{
			ec.push_back(m_meshH->get<Edge>(e_loc_id));
		}
	}

	// Store the edges in ee vector in the same side of the ones in ec
	// This way, ee[0] will be the closest end edge to ec[0], and ee[1] the closest to ec[1]
	// ---->
	std::vector<TCellID> ec_0_faces = n_neighbourhood.adjFacesToEdge(ec[0].id());
	if (var_front_edges_classification->value(n_neighbourhood.nextEdgeOfFace(ec_0_faces[0], ec[0].id()))==2)
	{
		ee.push_back(m_meshH->get<Edge>(n_neighbourhood.nextEdgeOfFace(ec_0_faces[0], ec[0].id())));
	}
	else
	{
		ee.push_back(m_meshH->get<Edge>(n_neighbourhood.nextEdgeOfFace(ec_0_faces[1], ec[0].id())));
	}

	std::vector<TCellID> ec_1_faces = n_neighbourhood.adjFacesToEdge(ec[1].id());
	if (var_front_edges_classification->value(n_neighbourhood.nextEdgeOfFace(ec_1_faces[0], ec[1].id()))==2)
	{
		ee.push_back(m_meshH->get<Edge>(n_neighbourhood.nextEdgeOfFace(ec_1_faces[0], ec[1].id())));
	}
	else
	{
		ee.push_back(m_meshH->get<Edge>(n_neighbourhood.nextEdgeOfFace(ec_1_faces[1], ec[1].id())));
	}
	// <----

	Node nc0 = ec[0].getOppositeNode(n);
	Node nc1 = ec[1].getOppositeNode(n);

	// Get the 5 faces around the node n
	TCellID f0_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec[0].id(), ee[0].id(), ee[1].id());
	TCellID f1_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec[0].id(), ee[1].id(), ee[0].id());
	TCellID f2_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ee[1].id(), ec[0].id(), ec[1].id());
	TCellID f3_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ee[1].id(), ec[1].id(), ec[0].id());
	TCellID f4_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec[1].id(), ee[0].id(), ee[1].id());
	TCellID f5_id = n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ee[0].id(), ec[1].id(), ec[0].id());

	// Get the 2 other nodes of Face f1 (the first ones are n and nc0) ---->
	TCellID f1_next_edge = n_neighbourhood.nextEdgeOfFace(f1_id, ec[0].id());
	TCellID n_f1_1_id = (m_meshH->get<Edge>(f1_next_edge)).getOppositeNodeId(n_id) ;
	Node n_f1_1 = m_meshH->get<Node>(n_f1_1_id);

	Edge f1_opp_edge = math::Utils::oppositeEdgeInFace(m_meshH, ec[0].id(), f1_id);
	Node n_f1_2 = f1_opp_edge.getOppositeNode(n_f1_1);
	// <----

	// Get the 2 other nodes of Face f4 (the first ones are n and nc1) ---->
	TCellID f4_next_edge = n_neighbourhood.nextEdgeOfFace(f4_id, ec[1].id());
	TCellID n_f4_1_id = (m_meshH->get<Edge>(f4_next_edge)).getOppositeNodeId(n_id) ;
	Node n_f4_1 = m_meshH->get<Node>(n_f4_1_id);

	Edge f4_opp_edge = math::Utils::oppositeEdgeInFace(m_meshH, ec[1].id(), f4_id);
	Node n_f4_2 = f4_opp_edge.getOppositeNode(n_f4_1);
	// <----

	Node n1 = m_meshH->get<Node>( m_FaceInfo[f0_id].next_ideal_nodes[n_id] ) ;
	Node n2 = m_meshH->get<Node>( m_FaceInfo[f1_id].next_ideal_nodes[n_f1_2.id()] );
	Node n3 = m_meshH->get<Node>( m_FaceInfo[f1_id].next_ideal_nodes[n_f1_1.id()] );
	Node n4 = m_meshH->get<Node>( m_FaceInfo[f4_id].next_ideal_nodes[n_f4_2.id()] );
	Node n5 = m_meshH->get<Node>( m_FaceInfo[f4_id].next_ideal_nodes[n_f4_1.id()] );

	// Create the two new hexas
	TCellID r_id = math::Utils::CreateHexaNConnectivities(m_meshH, n, nc0, n1, nc1, n_f1_1, n_f1_2, n2, n3);
	hexas_id.push_back(r_id);
	r_id = math::Utils::CreateHexaNConnectivities(m_meshH, n, nc0, n1, nc1, n_f4_1, n5, n4, n_f4_2);
	hexas_id.push_back(r_id);


	//===================//
	//      UPDATES      //
	//===================//

	// Update the faces and edges that can't create hexas anymore
	m_meshH->mark(ec[0], mark_edgesTreated);
	m_meshH->mark(ec[1], mark_edgesTreated);
	m_meshH->mark(m_meshH->get<Face>(f1_id), mark_facesTreated);
	m_meshH->mark(m_meshH->get<Face>(f4_id), mark_facesTreated);

	// Update the edge struc info for the 2 END edges ---->
	m_EdgeInfo[ee[0].id()].END_n_face_created[n_id] = true;
	m_EdgeInfo[ee[1].id()].END_n_face_created[n_id] = true;
	m_EdgeInfo[ee[0].id()].diag_next_node[n_id] = n5.id();
	m_EdgeInfo[ee[1].id()].diag_next_node[n_id] = n3.id();
	// <----

	// Update the edge struc info for the first CORNER edge ---->
	NodeNeighbourhoodOnFront_3D n_c0_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, nc0.id()) ;
	n_c0_neighbourhood.execute();
	for (auto edge_id:n_c0_neighbourhood.getOrderedEdges()) {
		if (var_front_edges_classification->value(edge_id) == 1 && edge_id != ec[0].id()) {
			m_EdgeInfo[edge_id].CORNER_n_face_created[nc0.id()] = true;
			m_EdgeInfo[edge_id].diag_next_node[nc0.id()] = n4.id();

			TCellID edge_f1_last = n_c0_neighbourhood.nextEdgeOfFace(f1_id, ec[0].id());
			for (auto faces_id : n_c0_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[0].id(), edge_id, edge_f1_last)) {
				m_FaceInfo[faces_id].next_nodes[nc0.id()] = n5.id();
			}
			TCellID edge_f0 = n_c0_neighbourhood.nextEdgeOfFace(f0_id, ec[0].id());
			for (auto faces_id : n_c0_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[0].id(), edge_id, edge_f0)) {
				m_FaceInfo[faces_id].next_nodes[nc0.id()] = n1.id();
			}

			std::pair<TCellID, TCellID> pair_1(n_c0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[0].id(), edge_f1_last), nc0.id());
			std::pair<TCellID, TCellID> pair_2(n_c0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[0].id(), edge_f0), nc0.id());

			m_EdgeInfo[edge_id].CORNER_next_nodes[pair_1] = n5.id();
			m_EdgeInfo[edge_id].CORNER_next_nodes[pair_2] = n1.id();
		}
	}
	// <----

	// Update the edge struc info for the second CORNER edge ---->
	NodeNeighbourhoodOnFront_3D n_c1_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, nc1.id()) ;
	n_c1_neighbourhood.execute();
	for (auto edge_id:n_c1_neighbourhood.getOrderedEdges()) {
		if (var_front_edges_classification->value(edge_id) == 1 && edge_id != ec[1].id()) {
			m_EdgeInfo[edge_id].CORNER_n_face_created[nc1.id()] = true;
			m_EdgeInfo[edge_id].diag_next_node[nc1.id()] = n2.id();

			TCellID edge_f4_last = n_c1_neighbourhood.nextEdgeOfFace(f4_id, ec[1].id());
			for (auto faces_id : n_c1_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[1].id(), edge_id, edge_f4_last)) {
				m_FaceInfo[faces_id].next_nodes[nc1.id()] = n3.id();
			}
			TCellID edge_f3 = n_c1_neighbourhood.nextEdgeOfFace(f3_id, ec[1].id());
			for (auto faces_id : n_c1_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(ec[1].id(), edge_id, edge_f3)) {
				m_FaceInfo[faces_id].next_nodes[nc1.id()] = n1.id();
			}

			std::pair<TCellID, TCellID> pair_1(n_c1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[1].id(), edge_f4_last), nc1.id());
			std::pair<TCellID, TCellID> pair_2(n_c1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(edge_id, ec[1].id(), edge_f3), nc1.id());

			m_EdgeInfo[edge_id].CORNER_next_nodes[pair_1] = n3.id();
			m_EdgeInfo[edge_id].CORNER_next_nodes[pair_2] = n1.id();
		}
	}
	// <----

	// Write the new hexa in a VTK file <----
	gmds::IGMeshIOService ioService(m_meshH);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("AeroExtrusion_3D_"+std::to_string(m_iteration)+".vtk");
	m_iteration++;
	//<----

	return hexas_id;
}
/*------------------------------------------------------------------------*/
std::vector<TCellID>
AeroExtrusion_3D::TemplateNode2Corner1Reversal(Front_3D &AFront, TCellID n_id, double dc)
{
	//std::cout << "Template Node 2 Corner 1 Reversal au noeud " << n_id << std::endl;
	std::vector<TCellID> hexas_id;

	Node n = m_meshH->get<Node>(n_id);
	NodeNeighbourhoodOnFront_3D n_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n_id);
	n_neighbourhood.execute();
	std::vector<TCellID> n_ordered_edges = n_neighbourhood.getOrderedEdges();
	std::vector<TCellID> n_ordered_faces = n_neighbourhood.getOrderedFaces();

	Variable<int>* var_front_edges_classification = m_meshH->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
	std::vector<Edge> ec;
	Edge er;
	for (auto e_loc_id:n_ordered_edges)
	{
		if (var_front_edges_classification->value(e_loc_id)==1)
		{
			ec.push_back(m_meshH->get<Edge>(e_loc_id));
		}
		else if (var_front_edges_classification->value(e_loc_id)==3)
		{
			er = m_meshH->get<Edge>(e_loc_id);
		}
	}

	Node n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11;
	n2 = m_meshH->get<Node>(m_FaceInfo[n_ordered_faces[0]].next_ideal_nodes[n_id]);

	// We build the 10 missing nodes ---->
	Node nr = er.getOppositeNode(n);
	Node nc0 = ec[0].getOppositeNode(n);
	Node nc1 = ec[1].getOppositeNode(n);

	math::Vector3d vr = (n.point()-nr.point()).normalize() ;
	math::Vector3d vc0 = (n.point()-nc0.point()).normalize() ;
	math::Vector3d vc1 = (n.point()-nc1.point()).normalize() ;

	Face f_rc0 = m_meshH->get<Face>(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(er.id(), ec[0].id(), ec[1].id()));
	Face f_rc1 = m_meshH->get<Face>(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(er.id(), ec[1].id(), ec[0].id()));
	Face f_c0c1 = m_meshH->get<Face>(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ec[0].id(), ec[1].id(), er.id()));

	math::Vector3d v_rc0 = f_rc0.normal();
	gmds::Cell::Data data_rc0 = m_fl.find(f_rc0.center());
	Node n_closest_f_rc0 = m_meshT->get<Node>(data_rc0.id);
	math::Vector3d v_test = m_VectorField->value(n_closest_f_rc0.id()).normalize() ;
	if (v_rc0.dot(v_test) <= 0)
	{
		v_rc0 = - v_rc0;
	}

	math::Vector3d v_rc1 = f_rc1.normal();
	gmds::Cell::Data data_rc1 = m_fl.find(f_rc1.center());
	Node n_closest_f_rc1 = m_meshT->get<Node>(data_rc1.id);
	v_test = m_VectorField->value(n_closest_f_rc1.id()).normalize() ;
	if (v_rc1.dot(v_test) <= 0)
	{
		v_rc1 = - v_rc1;
	}

	math::Vector3d v_c0c1 = f_c0c1.normal();
	gmds::Cell::Data data_c0c1 = m_fl.find(f_c0c1.center());
	Node n_closest_f_c0c1 = m_meshT->get<Node>(data_c0c1.id);
	v_test = m_VectorField->value(n_closest_f_c0c1.id()).normalize() ;
	if (v_c0c1.dot(v_test) <= 0)
	{
		v_c0c1 = - v_c0c1;
	}

	math::Point p1 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, vr);
	math::Point p3 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, (vc0+vc1).normalize());
	math::Point p7 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, v_rc0);
	math::Point p4 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, (vr+v_rc0).normalize());
	math::Point p6 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, (vc0+vc1+v_rc0).normalize());
	math::Point p5 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, (vr+vc0+vc1+v_rc0).normalize());

	math::Point p11 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, (v_rc1).normalize());
	math::Point p8 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, (vr+v_rc1).normalize());
	math::Point p10 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, (vc0+vc1+v_rc1).normalize());
	math::Point p9 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, (vr+vc0+vc1+v_rc1).normalize());

	n1 = m_meshH->newNode(p1);
	n3 = m_meshH->newNode(p3);
	n7 = m_meshH->newNode(p7);
	n4 = m_meshH->newNode(p4);
	n6 = m_meshH->newNode(p6);
	n5 = m_meshH->newNode(p5);
	n11 = m_meshH->newNode(p11);
	n8 = m_meshH->newNode(p8);
	n10 = m_meshH->newNode(p10);
	n9 = m_meshH->newNode(p9);
	// <----

	// Create the two new hexas
	TCellID r_id = math::Utils::CreateHexaNConnectivities(m_meshH, n, n1, n2, n3, n7, n4, n5, n6);
	hexas_id.push_back(r_id);
	r_id = math::Utils::CreateHexaNConnectivities(m_meshH, n, n1, n2, n3, n11,n8, n9, n10);
	hexas_id.push_back(r_id);

	// Update the layer ID for the 10 new nodes on the next layer ---->
	Variable<int>* var_node_couche_id = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	var_node_couche_id->set(n1.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n2.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n3.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n4.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n5.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n6.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n7.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n8.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n9.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n10.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n11.id(), AFront.getFrontID()+1);
	// <----

	// Update the layer ID for the 6 new faces on the next layer ---->
	Variable<int>* var_face_couche_id = m_meshH->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");
	TCellID f_1_id = math::Utils::CommonFace3Nodes(m_meshH, n4.id(), n5.id(), n6.id());
	TCellID f_2_id = math::Utils::CommonFace3Nodes(m_meshH, n1.id(),  n2.id(), n5.id());
	TCellID f_3_id = math::Utils::CommonFace3Nodes(m_meshH, n3.id(),  n2.id(), n5.id());
	TCellID f_4_id = math::Utils::CommonFace3Nodes(m_meshH, n3.id(),  n2.id(), n9.id());
	TCellID f_5_id = math::Utils::CommonFace3Nodes(m_meshH, n1.id(),  n2.id(), n9.id());
	TCellID f_6_id = math::Utils::CommonFace3Nodes(m_meshH, n8.id(),  n9.id(), n10.id());
	var_face_couche_id->set(f_1_id, AFront.getFrontID()+1);
	var_face_couche_id->set(f_2_id, AFront.getFrontID()+1);
	var_face_couche_id->set(f_3_id, AFront.getFrontID()+1);
	var_face_couche_id->set(f_4_id, AFront.getFrontID()+1);
	var_face_couche_id->set(f_5_id, AFront.getFrontID()+1);
	var_face_couche_id->set(f_6_id, AFront.getFrontID()+1);
	// <----

	// Update the edge structure info <----
	std::pair<TCellID, TCellID> pair_rc0(f_rc0.id(), n.id());
	std::pair<TCellID, TCellID> pair_rc1(f_rc1.id(), n.id());
	std::pair<TCellID, TCellID> pair_c0c1(f_c0c1.id(), n.id());
	m_EdgeInfo[er.id()].REVERSAL_n_faces_created[n_id] = true;
	m_EdgeInfo[er.id()].REVERSAL_diag_nodes[pair_rc0] = n6.id();
	m_EdgeInfo[er.id()].REVERSAL_adj_nodes[pair_rc0] = n7.id();
	m_EdgeInfo[er.id()].REVERSAL_diag_nodes[pair_rc1] = n10.id();
	m_EdgeInfo[er.id()].REVERSAL_adj_nodes[pair_rc1] = n11.id();
	m_EdgeInfo[er.id()].REVERSAL_medium_node[n_id] = n3.id() ;

	m_EdgeInfo[ec[0].id()].CORNER_n_face_created[n_id] = true;
	m_EdgeInfo[ec[0].id()].diag_next_node[n_id] = n4.id() ;
	m_EdgeInfo[ec[0].id()].CORNER_next_nodes[pair_rc0] = n7.id() ;
	m_EdgeInfo[ec[0].id()].CORNER_next_nodes[pair_c0c1] = n1.id() ;

	m_EdgeInfo[ec[1].id()].CORNER_n_face_created[n_id] = true;
	m_EdgeInfo[ec[1].id()].diag_next_node[n_id] = n8.id() ;
	m_EdgeInfo[ec[1].id()].CORNER_next_nodes[pair_rc1] = n11.id() ;
	m_EdgeInfo[ec[1].id()].CORNER_next_nodes[pair_c0c1] = n1.id() ;
	// ---->

	// Update the faces info ---->
	m_FaceInfo[f_rc0.id()].next_nodes[n_id] = n7.id();
	m_FaceInfo[f_rc1.id()].next_nodes[n_id] = n11.id();
	m_FaceInfo[f_c0c1.id()].next_nodes[n_id] = n1.id();
	// <----

	// Write the new hexa in a VTK file <----
	gmds::IGMeshIOService ioService(m_meshH);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("AeroExtrusion_3D_"+std::to_string(m_iteration)+".vtk");
	m_iteration++;
	//<----

	return hexas_id;
}
/*------------------------------------------------------------------------*/
std::vector<TCellID>
AeroExtrusion_3D::TemplateNode2End1Reversal(Front_3D &AFront, TCellID n_id, double dc, TInt mark_edgesTreated, TInt mark_facesTreated)
{
	//std::cout << "Template Node 2 End 1 Reversal au noeud " << n_id << std::endl;
	std::vector<TCellID> hexas_id;

	Node n = m_meshH->get<Node>(n_id);
	NodeNeighbourhoodOnFront_3D n_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n_id);
	n_neighbourhood.execute();
	std::vector<TCellID> n_ordered_edges = n_neighbourhood.getOrderedEdges();
	std::vector<TCellID> n_ordered_faces = n_neighbourhood.getOrderedFaces();

	Variable<int>* var_front_edges_classification = m_meshH->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
	std::vector<Edge> ee;
	Edge er;
	for (auto e_loc_id:n_ordered_edges)
	{
		if (var_front_edges_classification->value(e_loc_id)==2)
		{
			ee.push_back(m_meshH->get<Edge>(e_loc_id));
		}
		else if (var_front_edges_classification->value(e_loc_id)==3)
		{
			er = m_meshH->get<Edge>(e_loc_id);
		}
	}

	// Get the nodes for the 2 new hexas
	//Node n2 = m_meshH->get<Node>(m_FaceInfo[n_ordered_faces[0]].next_ideal_nodes[n_id]);
	Node nr = er.getOppositeNode(n);

	Face f_re0 = m_meshH->get<Face>(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(er.id(), ee[0].id(), ee[1].id()));
	Face f_re1 = m_meshH->get<Face>(n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(er.id(), ee[1].id(), ee[0].id()));

	Face f_0 = m_meshH->get<Face>( n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ee[0].id(), ee[1].id(), er.id()));
	Face f_1 = m_meshH->get<Face>( n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(ee[1].id(), ee[0].id(), er.id()));

	Edge e_0 = m_meshH->get<Edge>( n_neighbourhood.nextEdgeOfFace(f_0.id(), ee[0].id()));
	Edge e_1 = m_meshH->get<Edge>( n_neighbourhood.nextEdgeOfFace(f_1.id(), ee[1].id()));

	Face f_2 = m_meshH->get<Face>( n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e_0.id(), e_1.id(), ee[0].id()) );
	Face f_3 = m_meshH->get<Face>( n_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(e_1.id(), e_0.id(), ee[1].id()) );

	Edge e_2 = m_meshH->get<Edge>(n_neighbourhood.nextEdgeOfFace(f_2.id(), e_0.id()));

	Node n1 = e_2.getOppositeNode(n);
	Node n4 = e_0.getOppositeNode(n);
	Node n8 = e_1.getOppositeNode(n);

	Edge e_2_opp_f2 = math::Utils::oppositeEdgeInFace(m_meshH, e_2.id(), f_2.id()) ;
	Edge e_2_opp_f3 = math::Utils::oppositeEdgeInFace(m_meshH, e_2.id(), f_3.id()) ;

	Node n5 = e_2_opp_f2.getOppositeNode(n4);
	Node n9 = e_2_opp_f3.getOppositeNode(n8);

	Node n2 = m_meshH->get<Node>(m_FaceInfo[f_2.id()].next_ideal_nodes[n1.id()]);
	Node n3 = m_meshH->get<Node>(m_FaceInfo[f_2.id()].next_ideal_nodes[n4.id()]);
	Node n6 = m_meshH->get<Node>(m_FaceInfo[f_2.id()].next_ideal_nodes[n5.id()]);
	Node n7 = m_meshH->get<Node>(m_FaceInfo[f_3.id()].next_ideal_nodes[n8.id()]);
	Node n10 = m_meshH->get<Node>(m_FaceInfo[f_3.id()].next_ideal_nodes[n9.id()]);
	// <----

	// Create the two new hexas
	TCellID r_id = math::Utils::CreateHexaNConnectivities(m_meshH, nr, n, n1, n2, n3, n4, n5, n6);
	hexas_id.push_back(r_id);
	r_id = math::Utils::CreateHexaNConnectivities(m_meshH, nr, n, n1, n2, n7, n8, n9, n10);
	hexas_id.push_back(r_id);

	// Mark the faces and the edge of the front where the hexa insertion is not possible anymore
	m_meshH->mark(f_2, mark_facesTreated);
	m_meshH->mark(f_3, mark_facesTreated);
	m_meshH->mark(er, mark_edgesTreated);

	// Update the EdgeInfo for the two END edges ---->
	m_EdgeInfo[ee[0].id()].END_n_face_created[n_id] = true;
	m_EdgeInfo[ee[1].id()].END_n_face_created[n_id] = true;
	m_EdgeInfo[ee[0].id()].diag_next_node[n_id] = n3.id() ;
	m_EdgeInfo[ee[1].id()].diag_next_node[n_id] = n7.id() ;
	// <----

	// Update the EdgeInfo for the next REVERSAL edge ---->
	NodeNeighbourhoodOnFront_3D nr_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, nr.id());
	nr_neighbourhood.execute();
	Edge er_next;
	for (auto e_loc_id:nr_neighbourhood.getOrderedEdges())
	{
		if (e_loc_id != er.id()
		    && var_front_edges_classification->value(e_loc_id) == 3)
		{
			er_next = m_meshH->get<Edge>(e_loc_id);
		}
	}
	Edge e_opp_e0_fre0 = math::Utils::oppositeEdgeInFace(m_meshH, ee[0].id(), f_re0.id());
	Edge e_opp_e1_fre1 = math::Utils::oppositeEdgeInFace(m_meshH, ee[1].id(), f_re1.id());

	Face f_0_er_next = m_meshH->get<Face>( nr_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(er_next.id(), e_opp_e0_fre0.id(), e_opp_e1_fre1.id()) ) ;
	Face f_1_er_next = m_meshH->get<Face>( nr_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(er_next.id(), e_opp_e1_fre1.id(), e_opp_e0_fre0.id()) ) ;

	std::pair<TCellID, TCellID> pair_0(f_0_er_next.id(), nr.id());
	std::pair<TCellID, TCellID> pair_1(f_1_er_next.id(), nr.id());

	m_EdgeInfo[er_next.id()].REVERSAL_n_faces_created[nr.id()] = true;
	m_EdgeInfo[er_next.id()].REVERSAL_adj_nodes[pair_0] = n3.id();
	m_EdgeInfo[er_next.id()].REVERSAL_adj_nodes[pair_1] = n7.id();
	m_EdgeInfo[er_next.id()].REVERSAL_diag_nodes[pair_0] = n6.id();
	m_EdgeInfo[er_next.id()].REVERSAL_diag_nodes[pair_1] = n10.id();
	m_EdgeInfo[er_next.id()].REVERSAL_medium_node[nr.id()] = n2.id();
	// <----

	// Update the FaceInfo for the faces around the next reversal edge ---->
	for (auto f_loc:nr_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(er_next.id(), er.id(), e_opp_e1_fre1.id()))
	{
		m_FaceInfo[f_loc].next_nodes[nr.id()] = n3.id();
	}
	for (auto f_loc:nr_neighbourhood.facesBtwEdge1nEdge2AvoidingEdge3(er_next.id(), er.id(), e_opp_e0_fre0.id()))
	{
		m_FaceInfo[f_loc].next_nodes[nr.id()] = n7.id();
	}
	// <----

	// Write the new hexa in a VTK file <----
	gmds::IGMeshIOService ioService(m_meshH);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("AeroExtrusion_3D_"+std::to_string(m_iteration)+".vtk");
	m_iteration++;
	//<----

	return hexas_id;
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
AeroExtrusion_3D::TemplateEdgeEnd(Front_3D &AFront, TCellID e_id, double dc, TInt mark_edgesTreated, TInt mark_facesTreated)
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
std::vector<TCellID>
AeroExtrusion_3D::TemplateEdgeReversal(Front_3D &AFront, TCellID e_id, double dc, TInt mark_edgesTreated, TInt mark_facesTreated)
{
	//std::cout << "Template Edge Reversal " << e_id << std::endl;
	std::vector<TCellID> hexas_id;

	Node n0 = ((m_meshH->get<Edge>(e_id)).get<Node>())[0];
	Node n1 = ((m_meshH->get<Edge>(e_id)).get<Node>())[1];

	Node n2, n3, n4, n5, n6, n7, n8, n9, n10, n11;

	std::vector<TCellID> e_faces = AFront.edgeFacesOnFront(m_meshH, e_id);

	// First side
	if (m_EdgeInfo[e_id].REVERSAL_n_faces_created[n0.id()])	// The faces are already created
	{
		std::pair<TCellID, TCellID> pair_1(e_faces[0], n0.id());
		n4 = m_meshH->get<Node>(m_EdgeInfo[e_id].REVERSAL_adj_nodes[pair_1]);
		n7 = m_meshH->get<Node>(m_EdgeInfo[e_id].REVERSAL_diag_nodes[pair_1]);
		std::pair<TCellID, TCellID> pair_2(e_faces[1], n0.id());
		n8 = m_meshH->get<Node>(m_EdgeInfo[e_id].REVERSAL_adj_nodes[pair_2]);
		n11 = m_meshH->get<Node>(m_EdgeInfo[e_id].REVERSAL_diag_nodes[pair_2]);

		n3 = m_meshH->get<Node>(m_EdgeInfo[e_id].REVERSAL_medium_node[n0.id()]);
	}
	else	// Create the faces
	{
		n3 = m_meshH->get<Node>(m_FaceInfo[e_faces[0]].next_ideal_nodes[n0.id()]) ;
		math::Vector3d v = (n3.point()-n0.point()).normalize() ;

		math::Point f0_barycentre = m_meshH->get<Face>(e_faces[0]).center() ;
		math::Vector3d f0_normale = m_meshH->get<Face>(e_faces[0]).normal();
		gmds::Cell::Data data = m_fl.find(f0_barycentre);
		Node n_closest = m_meshT->get<Node>(data.id);
		math::Vector3d v1 = m_VectorField->value(n_closest.id()).normalize() ;
		if (f0_normale.dot(v1) <= 0)
		{
			f0_normale = - f0_normale;
		}

		math::Point p4 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n0.point(), dc, m_DistanceField, f0_normale);
		math::Point p7 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n0.point(), dc, m_DistanceField, (f0_normale+v).normalize());

		n4 = m_meshH->newNode(p4);
		n7 = m_meshH->newNode(p7);

		math::Point f1_barycentre = m_meshH->get<Face>(e_faces[1]).center() ;
		math::Vector3d f1_normale = m_meshH->get<Face>(e_faces[1]).normal();
		gmds::Cell::Data data_f1 = m_fl.find(f1_barycentre);
		Node n_closest_f1 = m_meshT->get<Node>(data_f1.id);
		math::Vector3d v_f1 = m_VectorField->value(n_closest_f1.id()).normalize() ;
		if (f1_normale.dot(v_f1) <= 0)
		{
			f1_normale = - f1_normale;
		}

		math::Point p8 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n0.point(), dc, m_DistanceField, f1_normale);
		math::Point p11 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n0.point(), dc, m_DistanceField, (f1_normale+v).normalize());

		n8 = m_meshH->newNode(p8);
		n11 = m_meshH->newNode(p11);

		//=====================//
		// Update the REVERSAL //
		// edge info for the	  //
		// next edge			  //
		//=====================//
		Variable<int>* edge_classification = m_meshH->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
		NodeNeighbourhoodOnFront_3D n0_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n0.id());
		n0_neighbourhood.execute();
		Edge next_edge;
		for (auto e_loc_id:n0_neighbourhood.getOrderedEdges())
		{
			if ((e_loc_id != e_id)
			    && edge_classification->value(e_loc_id)==3 )
			{
				next_edge = m_meshH->get<Edge>(e_loc_id);
			}
		}

		Edge e_f0 = m_meshH->get<Edge>(n0_neighbourhood.nextEdgeOfFace(e_faces[0], e_id));
		Edge e_f1 = m_meshH->get<Edge>(n0_neighbourhood.nextEdgeOfFace(e_faces[1], e_id));

		Face f0_next = m_meshH->get<Face>(n0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(next_edge.id(), e_f0.id(), e_f1.id())) ;
		Face f1_next = m_meshH->get<Face>(n0_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(next_edge.id(), e_f1.id(), e_f0.id())) ;

		std::pair<TCellID, TCellID> pair_0(f0_next.id(), n0.id());
		std::pair<TCellID, TCellID> pair_1(f1_next.id(), n0.id());

		m_EdgeInfo[next_edge.id()].REVERSAL_n_faces_created[n0.id()] = true;
		m_EdgeInfo[next_edge.id()].REVERSAL_medium_node[n0.id()] = n3.id();
		m_EdgeInfo[next_edge.id()].REVERSAL_adj_nodes[pair_0] = n4.id();
		m_EdgeInfo[next_edge.id()].REVERSAL_adj_nodes[pair_1] = n8.id();
		m_EdgeInfo[next_edge.id()].REVERSAL_diag_nodes[pair_0] = n7.id();
		m_EdgeInfo[next_edge.id()].REVERSAL_diag_nodes[pair_1] = n11.id();
		//==============//
		// END UPDATES  //
		//==============//

	}

	// Second side
	if (m_EdgeInfo[e_id].REVERSAL_n_faces_created[n1.id()])	// The faces are already created
	{
		std::pair<TCellID, TCellID> pair_1(e_faces[0], n1.id());
		n5 = m_meshH->get<Node>(m_EdgeInfo[e_id].REVERSAL_adj_nodes[pair_1]);
		n6 = m_meshH->get<Node>(m_EdgeInfo[e_id].REVERSAL_diag_nodes[pair_1]);
		std::pair<TCellID, TCellID> pair_2(e_faces[1], n1.id());
		n9 = m_meshH->get<Node>(m_EdgeInfo[e_id].REVERSAL_adj_nodes[pair_2]);
		n10 = m_meshH->get<Node>(m_EdgeInfo[e_id].REVERSAL_diag_nodes[pair_2]);

		n2 = m_meshH->get<Node>(m_EdgeInfo[e_id].REVERSAL_medium_node[n1.id()]);
	}
	else	// Create the faces
	{
		n2 = m_meshH->get<Node>(m_FaceInfo[e_faces[0]].next_ideal_nodes[n1.id()]) ;
		math::Vector3d v = (n2.point()-n1.point()).normalize() ;

		math::Point f0_barycentre = m_meshH->get<Face>(e_faces[0]).center() ;
		math::Vector3d f0_normale = m_meshH->get<Face>(e_faces[0]).normal();
		gmds::Cell::Data data = m_fl.find(f0_barycentre);
		Node n_closest = m_meshT->get<Node>(data.id);
		math::Vector3d v1 = m_VectorField->value(n_closest.id()).normalize() ;
		if (f0_normale.dot(v1) <= 0)
		{
			f0_normale = - f0_normale;
		}

		math::Point p5 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n1.point(), dc, m_DistanceField, f0_normale);
		math::Point p6 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n1.point(), dc, m_DistanceField, (f0_normale+v).normalize());

		n5 = m_meshH->newNode(p5);
		n6 = m_meshH->newNode(p6);

		math::Point f1_barycentre = m_meshH->get<Face>(e_faces[1]).center() ;
		math::Vector3d f1_normale = m_meshH->get<Face>(e_faces[1]).normal();
		gmds::Cell::Data data_f1 = m_fl.find(f1_barycentre);
		Node n_closest_f1 = m_meshT->get<Node>(data_f1.id);
		math::Vector3d v_f1 = m_VectorField->value(n_closest_f1.id()).normalize() ;
		if (f1_normale.dot(v_f1) <= 0)
		{
			f1_normale = - f1_normale;
		}

		math::Point p9 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n1.point(), dc, m_DistanceField, f1_normale);
		math::Point p10 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n1.point(), dc, m_DistanceField, (f1_normale+v).normalize());

		n9 = m_meshH->newNode(p9);
		n10 = m_meshH->newNode(p10);

		//=====================//
		// Update the REVERSAL //
		// edge info for the	  //
		// next edge			  //
		//=====================//
		Variable<int>* edge_classification = m_meshH->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
		NodeNeighbourhoodOnFront_3D n1_neighbourhood = NodeNeighbourhoodOnFront_3D(m_meshH, &AFront, n1.id());
		n1_neighbourhood.execute();
		Edge next_edge;
		for (auto e_loc_id:n1_neighbourhood.getOrderedEdges())
		{
			if ((e_loc_id != e_id)
			    && edge_classification->value(e_loc_id)==3 )
			{
				next_edge = m_meshH->get<Edge>(e_loc_id);
			}
		}

		Edge e_f0 = m_meshH->get<Edge>(n1_neighbourhood.nextEdgeOfFace(e_faces[0], e_id));
		Edge e_f1 = m_meshH->get<Edge>(n1_neighbourhood.nextEdgeOfFace(e_faces[1], e_id));

		Face f0_next = m_meshH->get<Face>(n1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(next_edge.id(), e_f0.id(), e_f1.id())) ;
		Face f1_next = m_meshH->get<Face>(n1_neighbourhood.adjFaceToEdge1InEdge2SideAvoidingEdge3(next_edge.id(), e_f1.id(), e_f0.id())) ;

		std::pair<TCellID, TCellID> pair_0(f0_next.id(), n1.id());
		std::pair<TCellID, TCellID> pair_1(f1_next.id(), n1.id());

		m_EdgeInfo[next_edge.id()].REVERSAL_n_faces_created[n1.id()] = true;
		m_EdgeInfo[next_edge.id()].REVERSAL_medium_node[n1.id()] = n2.id();
		m_EdgeInfo[next_edge.id()].REVERSAL_adj_nodes[pair_0] = n5.id();
		m_EdgeInfo[next_edge.id()].REVERSAL_adj_nodes[pair_1] = n9.id();
		m_EdgeInfo[next_edge.id()].REVERSAL_diag_nodes[pair_0] = n6.id();
		m_EdgeInfo[next_edge.id()].REVERSAL_diag_nodes[pair_1] = n10.id();
		//==============//
		// END UPDATES //
		//==============//
	}

	TCellID r_id = math::Utils::CreateHexaNConnectivities(m_meshH, n0, n1, n2, n3, n4, n5, n6, n7);
	hexas_id.push_back(r_id);
	r_id = math::Utils::CreateHexaNConnectivities(m_meshH, n0, n1, n2, n3, n8, n9, n10, n11);
	hexas_id.push_back(r_id);

	// Update the layer id of the new nodes ---->
	Variable<int>* var_node_couche_id = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	var_node_couche_id->set(n4.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n7.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n3.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n11.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n8.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n5.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n6.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n2.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n10.id(), AFront.getFrontID()+1);
	var_node_couche_id->set(n9.id(), AFront.getFrontID()+1);
	// <----

	// Update the layer ID for the 6 new faces on the next layer ---->
	Variable<int>* var_face_couche_id = m_meshH->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");
	TCellID f_1_id = math::Utils::CommonFace3Nodes(m_meshH, n4.id(), n5.id(), n6.id());
	TCellID f_2_id = math::Utils::CommonFace3Nodes(m_meshH, n2.id(), n3.id(), n7.id());
	TCellID f_3_id = math::Utils::CommonFace3Nodes(m_meshH, n2.id(), n3.id(), n11.id());
	TCellID f_4_id = math::Utils::CommonFace3Nodes(m_meshH, n9.id(), n10.id(), n11.id());
	var_face_couche_id->set(f_1_id, AFront.getFrontID()+1);
	var_face_couche_id->set(f_2_id, AFront.getFrontID()+1);
	var_face_couche_id->set(f_3_id, AFront.getFrontID()+1);
	var_face_couche_id->set(f_4_id, AFront.getFrontID()+1);
	// <----

	// Update the Face Info structures ---->
	m_FaceInfo[e_faces[0]].next_nodes[n0.id()] = n4.id();
	m_FaceInfo[e_faces[0]].next_nodes[n1.id()] = n5.id();
	m_FaceInfo[e_faces[1]].next_nodes[n0.id()] = n8.id();
	m_FaceInfo[e_faces[1]].next_nodes[n1.id()] = n9.id();
	// <----

	// Write the new hexa in a VTK file <----
	gmds::IGMeshIOService ioService(m_meshH);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	vtkWriter.write("AeroExtrusion_3D_"+std::to_string(m_iteration)+".vtk");
	m_iteration++;
	//<----

	return hexas_id;
}
/*------------------------------------------------------------------------*/
void
AeroExtrusion_3D::TemplateFace(TCellID f_id, Front_3D &Front_IN, std::map<TCellID, TCellID>& map_new_nodes)
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