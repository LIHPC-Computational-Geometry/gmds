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
#include <gmds/claire/LayerStructureManager_3D.h>
#include <gmds/claire/PatternFace.h>
#include <gmds/claire/PatternNode3Corner.h>
#include <gmds/claire/PatternNode2Corner1End.h>
#include <gmds/claire/PatternNode1Corner2End.h>
#include <gmds/claire/PatternNode3End.h>
#include <gmds/claire/PatternNode3Corner3End.h>
#include <gmds/claire/PatternNode2Corner2End.h>
#include <gmds/claire/PatternNode2Corner1Reversal.h>
#include <gmds/claire/PatternNode2End1Reversal.h>
#include <gmds/claire/PatternEdgeCorner.h>
#include <gmds/claire/PatternEdgeEnd.h>
#include <gmds/claire/PatternEdgeReversal.h>

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

	Front_3D Current_Front = Front_Geom;
	if (m_params_aero.insertions_allowed_on_first_layer)
	{
		// Compute the first layer of blocks with patterns
		Current_Front = ComputeLayer(Front_Geom, m_meshT->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"), m_params_aero.delta_cl,
		                             m_VectorField);
	}
	else
	{
		// Compute the first layer of blocks without patterns
		Current_Front = Compute1stLayer(Front_Geom, m_meshT->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"), m_VectorField);
	}

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

	// Variables
	Variable<int>* var_NODE_couche_id = m_meshH->getVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* var_front_edges_classification = m_meshH->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
	Variable<int>*var_front_NbrFeatureEdgesAroundNode = m_meshH->getOrCreateVariable<int, GMDS_NODE>("NbrFeatureEdgesAroundNode");


	// Mise à jour de l'indice de couche
	for (auto n_id:Front_IN.getNodes()){
		var_NODE_couche_id->set(map_new_nodes[n_id], Front_IN.getFrontID()+1);
	}

	// Front edges and nodes classification
	FrontEdgesNodesClassification_3D Classification = FrontEdgesNodesClassification_3D(m_meshH,
	                                                                                   &Front_IN,
	                                                                                   m_params_aero,
	                                                                                   var_front_edges_classification,
	                                                                                   var_front_NbrFeatureEdgesAroundNode,
	                                                                                   m_meshT,
	                                                                                   &m_fl,
	                                                                                   m_VectorField);
	Classification.execute();

	/*
	FrontPaths_3D Paths = FrontPaths_3D(m_meshH, &Front_IN, &Classification);
	Paths.execute();
	TInt mark_EdgesTemplates = Paths.getMarkEdgesForTemplates();
	TInt mark_NodesTemplates = Paths.getMarkNodesForTemplates();
	*/

	TInt mark_edgesTreated = m_meshH->newMark<Edge>();
	TInt mark_facesTreated = m_meshH->newMark<Face>();

	TInt mark_EdgesTemplates = Classification.getMarkEdgesTemplates();
	TInt mark_NodesTemplates = Classification.getMarkNodesTemplates();

	Variable<int>* var_TEST = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_TEST_NODES");

	LayerStructureManager_3D StructManager = LayerStructureManager_3D(m_meshH, &Front_IN, map_new_nodes);

	// Create the HEXA on each NODE classified
	if (m_params_aero.insertions_allowed)
	{
		std::map<TCellID, int> singular_nodes = getSingularNodes(Front_IN, var_front_edges_classification);
		for (auto singu_node : singular_nodes) {
			TCellID n_id = singu_node.first;
			int singu_type = singu_node.second;
			if (singu_type == 1 && m_meshH->isMarked(m_meshH->get<Node>(n_id), mark_NodesTemplates)) {
				std::cout << "Template NODE: 3 Corner" << std::endl;
				PatternNode3Corner p(m_meshH, &Front_IN, n_id, &StructManager, m_meshT, &m_fl, dist_cible, A_distance, A_vectors);
				p.execute();
				std::vector<Region> r_newHex = p.getNewHex();
				m_Patterns->set(r_newHex[0].id(), 1);
				if (m_params_aero.with_debug_files) {
					WriteVTKforDebug();
				}
				// TCellID r_id = TemplateNode3Corner(Front_IN, n_id, map_new_nodes, dist_cible, A_distance);
				// m_Patterns->set(r_id, 1);
			}
			else if (singu_type == 2 && m_meshH->isMarked(m_meshH->get<Node>(n_id), mark_NodesTemplates)) {
				// std::cout << "Template NODE: 2 Corner, 1 End" << std::endl;
				// TCellID r_id = TemplateNode2Corner1End(Front_IN, n_id, dist_cible, A_distance, mark_edgesTreated, mark_facesTreated);
				// m_Patterns->set(r_id, 2);
				PatternNode2Corner1End p(m_meshH, &Front_IN, n_id, &StructManager, m_meshT, &m_fl, dist_cible, A_distance, A_vectors);
				p.execute();
				std::vector<Region> r_newHex = p.getNewHex();
				m_Patterns->set(r_newHex[0].id(), 2);
				if (m_params_aero.with_debug_files) {
					WriteVTKforDebug();
				}
			}
			else if (singu_type == 3 && m_meshH->isMarked(m_meshH->get<Node>(n_id), mark_NodesTemplates)) {
				std::cout << "Template NODE: 1 Corner, 2 End" << std::endl;
				// TCellID r_id = TemplateNode1Corner2End(Front_IN, n_id, dist_cible, mark_edgesTreated, mark_facesTreated);
				// m_Patterns->set(r_id, 3);
				PatternNode1Corner2End p(m_meshH, &Front_IN, n_id, &StructManager, m_meshT, &m_fl, dist_cible, A_distance, A_vectors);
				p.execute();
				std::vector<Region> r_newHex = p.getNewHex();
				m_Patterns->set(r_newHex[0].id(), 3);
				if (m_params_aero.with_debug_files) {
					WriteVTKforDebug();
				}
			}
			else if (singu_type == 4 && m_meshH->isMarked(m_meshH->get<Node>(n_id), mark_NodesTemplates)) {
				std::cout << "Template NODE: 3 End" << std::endl;
				// TCellID r_id = TemplateNode3End(Front_IN, n_id, mark_edgesTreated, mark_facesTreated);
				// m_Patterns->set(r_id, 4);
				PatternNode3End p(m_meshH, &Front_IN, n_id, &StructManager, m_meshT, &m_fl, dist_cible, A_distance, A_vectors);
				p.execute();
				std::vector<Region> r_newHex = p.getNewHex();
				m_Patterns->set(r_newHex[0].id(), 4);
				if (m_params_aero.with_debug_files) {
					WriteVTKforDebug();
				}
			}
			else if (singu_type == 5 && m_meshH->isMarked(m_meshH->get<Node>(n_id), mark_NodesTemplates)) {
				std::cout << "Template NODE: 3 Corner, 3 End" << std::endl;
				// TCellID r_id = TemplateNode3Corner3End(Front_IN, n_id, dist_cible, A_distance, mark_edgesTreated, mark_facesTreated);
				// m_Patterns->set(r_id, 5);
				PatternNode3Corner3End p(m_meshH, &Front_IN, n_id, &StructManager, m_meshT, &m_fl, dist_cible, A_distance, A_vectors);
				p.execute();
				std::vector<Region> r_newHex = p.getNewHex();
				m_Patterns->set(r_newHex[0].id(), 5);
				if (m_params_aero.with_debug_files) {
					WriteVTKforDebug();
				}
			}
			else if (singu_type == 6 && m_meshH->isMarked(m_meshH->get<Node>(n_id), mark_NodesTemplates)) {
				std::cout << "Template NODE: 2 Corner, 2 End" << std::endl;
				// std::vector<TCellID> r_id = TemplateNode2Corner2End(Front_IN, n_id, dist_cible, mark_edgesTreated, mark_facesTreated); // Create 2 hexas
				// m_Patterns->set(r_id[0], 6);
				// m_Patterns->set(r_id[1], 6);
				PatternNode2Corner2End p(m_meshH, &Front_IN, n_id, &StructManager, m_meshT, &m_fl, dist_cible, A_distance, A_vectors);
				p.execute();
				std::vector<Region> r_newHex = p.getNewHex();
				m_Patterns->set(r_newHex[0].id(), 6);
				m_Patterns->set(r_newHex[1].id(), 6);
				if (m_params_aero.with_debug_files) {
					WriteVTKforDebug();
				}
			}
			else if (singu_type == 7 && m_meshH->isMarked(m_meshH->get<Node>(n_id), mark_NodesTemplates)) {
				std::cout << "Template NODE: 2 Corner, 1 Reversal" << std::endl;
				// std::vector<TCellID> r_id = TemplateNode2Corner1Reversal(Front_IN, n_id, dist_cible, A_distance); // Create 2 hexas
				// m_Patterns->set(r_id[0], 7);
				// m_Patterns->set(r_id[1], 7);
				PatternNode2Corner1Reversal p(m_meshH, &Front_IN, n_id, &StructManager, m_meshT, &m_fl, dist_cible, A_distance, A_vectors);
				p.execute();
				std::vector<Region> r_newHex = p.getNewHex();
				m_Patterns->set(r_newHex[0].id(), 7);
				m_Patterns->set(r_newHex[1].id(), 7);
				if (m_params_aero.with_debug_files) {
					WriteVTKforDebug();
				}
			}
			else if (singu_type == 8 && m_meshH->isMarked(m_meshH->get<Node>(n_id), mark_NodesTemplates)) {
				std::cout << "Template NODE: 2 End, 1 Reversal" << std::endl;
				// std::vector<TCellID> r_id = TemplateNode2End1Reversal(Front_IN, n_id, dist_cible, mark_edgesTreated, mark_facesTreated); // Create 2 hexas
				// m_Patterns->set(r_id[0], 8);
				// m_Patterns->set(r_id[1], 8);
				PatternNode2End1Reversal p(m_meshH, &Front_IN, n_id, &StructManager, m_meshT, &m_fl, dist_cible, A_distance, A_vectors);
				p.execute();
				std::vector<Region> r_newHex = p.getNewHex();
				m_Patterns->set(r_newHex[0].id(), 8);
				m_Patterns->set(r_newHex[1].id(), 8);
				if (m_params_aero.with_debug_files) {
					WriteVTKforDebug();
				}
			}
			if (m_meshH->isMarked(m_meshH->get<Node>(n_id), mark_NodesTemplates)) {
				// std::cout << "Noeud marqué pour template" << std::endl;
				var_TEST->set(n_id, 1);
			}
		}
	}

	// Create the HEXA on each EDGE classified
	if (m_params_aero.insertions_allowed)
	{
		std::map<TCellID, int> singular_edges = getSingularEdges(Front_IN, var_front_edges_classification, mark_edgesTreated);
		for (auto singu_edge : singular_edges) {
			TCellID e_id = singu_edge.first;
			int singu_type = singu_edge.second;
			if (singu_type == 1 && !StructManager.isEdgeTreated(e_id) && m_meshH->isMarked(m_meshH->get<Edge>(e_id), mark_EdgesTemplates)) {
				// std::cout << "Template EDGE: Corner" << std::endl;
				PatternEdgeCorner p(m_meshH, &Front_IN, e_id, &StructManager, m_meshT, &m_fl, dist_cible, A_distance, A_vectors);
				p.execute();
				std::vector<Region> r_newHex = p.getNewHex();
				m_Patterns->set(r_newHex[0].id(), -1);
				if (m_params_aero.with_debug_files) {
					WriteVTKforDebug();
				}
				// TCellID r_id = TemplateEdgeCorner(Front_IN, e_id, dist_cible, A_distance);
				// m_Patterns->set(r_id, -1);
			}
			else if (singu_type == 2 && !StructManager.isEdgeTreated(e_id) && m_meshH->isMarked(m_meshH->get<Edge>(e_id), mark_EdgesTemplates)) {
				// TCellID r_id = TemplateEdgeEnd(Front_IN, e_id, dist_cible, mark_edgesTreated, mark_facesTreated);
				// m_Patterns->set(r_id, -2);
				// std::cout << "Template EDGE: END" << std::endl;
				PatternEdgeEnd p(m_meshH, &Front_IN, e_id, &StructManager, m_meshT, &m_fl, dist_cible, A_distance, A_vectors);
				p.execute();
				std::vector<Region> r_newHex = p.getNewHex();
				m_Patterns->set(r_newHex[0].id(), -2);
				if (m_params_aero.with_debug_files) {
					WriteVTKforDebug();
				}
			}
			else if (singu_type == 3 && !StructManager.isEdgeTreated(e_id) && m_meshH->isMarked(m_meshH->get<Edge>(e_id), mark_EdgesTemplates)) {
				// std::vector<TCellID> r_id = TemplateEdgeReversal(Front_IN, e_id, dist_cible, A_distance, mark_edgesTreated, mark_facesTreated);
				// m_Patterns->set(r_id[0], -3);
				// m_Patterns->set(r_id[1], -3);
				// std::cout << "Template EDGE: Reversal" << std::endl;
				PatternEdgeReversal p(m_meshH, &Front_IN, e_id, &StructManager, m_meshT, &m_fl, dist_cible, A_distance, A_vectors);
				p.execute();
				std::vector<Region> r_newHex = p.getNewHex();
				m_Patterns->set(r_newHex[0].id(), -3);
				m_Patterns->set(r_newHex[1].id(), -3);
				if (m_params_aero.with_debug_files) {
					WriteVTKforDebug();
				}
			}
		}
	}

	// Create the HEXA on each FACE of the front
	for (auto f_id:Front_IN.getFaces()){
		Face f = m_meshH->get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		if (!StructManager.isFaceTreated(f_id))
		{
			PatternFace p_face(m_meshH, &Front_IN, f_id, &StructManager);
			p_face.execute();
			if (m_params_aero.with_debug_files)
			{
				WriteVTKforDebug();
			}
			//TemplateFace(f_id, Front_IN, map_new_nodes);
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
			std::cout << "ATTENTION: the layer is not valid. At least one edge of the front is not adjacent to 2 faces of the front." << std::endl;
			std::cout << "Edge... " << e_id << std::endl;
			std::cout << "Nodes... " << e.get<Node>()[0].id() << ", and " << e.get<Node>()[1].id() << std::endl;
			exit(1);
		}
	}

	return Front_OUT;
}
/*------------------------------------------------------------------------*/
Front_3D
AeroExtrusion_3D::Compute1stLayer(Front_3D Front_IN, Variable<double>* A_distance, Variable<math::Vector3d>* A_vectors)
{
	std::cout << "---------> build layer: " << Front_IN.getFrontID()+1 << std::endl;

	std::map<TCellID, TCellID> map_new_nodes = ComputeIdealPositions(Front_IN, m_params_aero.delta_cl, A_distance, A_vectors);

	LayerStructureManager_3D StructManager = LayerStructureManager_3D(m_meshH, &Front_IN, map_new_nodes);

	// Mise à jour de l'indice de couche
	Variable<int>*var_NODE_couche_id = m_meshH->getVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	for (auto n_id: Front_IN.getNodes()){
		var_NODE_couche_id->set(map_new_nodes[n_id], Front_IN.getFrontID()+1);
	}

	// Create the HEXA on each FACE of the front
	for (auto f_id: Front_IN.getFaces()){
		Face f = m_meshH->get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		//TemplateFace(f_id, Front_IN, map_new_nodes);
		PatternFace p_face(m_meshH, &Front_IN, f_id, &StructManager);
		p_face.execute();
		if (m_params_aero.with_debug_files)
		{
			WriteVTKforDebug();
		}
	}

	// Erase the nodes connected to nothing
	math::Utils::MeshCleaner(m_meshH);

	// Init the Front OUT
	Front_3D Front_OUT = InitFrontOUT(Front_IN);

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
void AeroExtrusion_3D::WriteVTKforDebug()
{
	gmds::IGMeshIOService ioService(m_meshH);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::R);
	vtkWriter.setDataOptions(gmds::N|gmds::R);
	//vtkWriter.write("AeroExtrusion_3D_"+std::to_string(m_iteration)+".vtk");
	m_iteration++;
}
/*------------------------------------------------------------------------*/