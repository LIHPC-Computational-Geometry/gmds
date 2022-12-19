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
#include <gmds/claire/AeroMeshQuality.h>

#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
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
	//Front_3D Current_Front = Front_Geom;

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

	// Init the multiple_info type
	for (auto f_id:Front_IN.getFaces())
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

	// Mise à jour de l'indice de couche
	Variable<int>* couche_id = m_meshH->getVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	for (auto n_id:Front_IN.getNodes()){
		couche_id->set(map_new_nodes[n_id], Front_IN.getFrontID()+1);
	}

	std::vector<TCellID> front_nodes = Front_IN.getNodes();
	std::vector<TCellID> front_faces = Front_IN.getFaces();

	Variable<int>* var_front_edges_classification = FrontEdgesClassification(Front_IN);

	int mark_edgesTreated = m_meshH->newMark<Edge>();
	int mark_facesTreated = m_meshH->newMark<Face>();

	std::map<TCellID, int> singular_nodes = getSingularNodes(Front_IN, var_front_edges_classification);

	for (auto singu_node:singular_nodes)
	{
		TCellID n_id = singu_node.first ;
		int singu_type = singu_node.second;
		if (singu_type==1)
		{
			TCellID r_id = TemplateNode3Corner(Front_IN, n_id, map_new_nodes, dist_cible);
		}
		else if (singu_type==4)
		{
			TCellID r_id = TemplateNode3End(Front_IN, n_id, mark_edgesTreated, mark_facesTreated);
		}
	}

	// Insert Hexa on Singular Edges
	std::map<TCellID, int> singular_edges = getSingularEdges(Front_IN, var_front_edges_classification, mark_edgesTreated);

	for (auto singu_edge:singular_edges)
	{
		TCellID e_id = singu_edge.first ;
		int singu_type = singu_edge.second;
		if (singu_type==1)
		{
			TCellID r_id = TemplateEdgeCorner(Front_IN, e_id, dist_cible);
		}
		else
		{
			//std::cout << "Edge singularity not implemented yet" << std::endl;
		}
	}

	m_meshH->unmarkAll<Edge>(mark_edgesTreated);
	m_meshH->freeMark<Edge>(mark_edgesTreated);
	m_meshH->unmarkAll<Face>(mark_facesTreated);
	m_meshH->freeMark<Face>(mark_facesTreated);


	// Ajout des hex restants
	for (auto f_id:front_faces){
		Face f = m_meshH->get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		TemplateFace(f_id, Front_IN, map_new_nodes);
	}

	// Supression des noeuds non utilisés
	math::Utils::MeshCleaner(m_meshH);



	// Initialisation du front de sortie
	std::vector<TCellID> new_front_nodes_id;
	std::vector<TCellID> new_front_faces_id;

	Variable<int>* var_node_couche_id = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* var_face_couche_id = m_meshH->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");


	for (auto n_id:m_meshH->nodes())
	{
		if (var_node_couche_id->value(n_id) == Front_IN.getFrontID()+1)
		{
			new_front_nodes_id.push_back(n_id);
		}
	}
	for (auto f_id:m_meshH->faces())
	{
		if (var_face_couche_id->value(f_id) == Front_IN.getFrontID()+1)
		{
			new_front_faces_id.push_back(f_id);
		}
	}
	Front_3D Front_OUT = Front_3D(Front_IN.getFrontID()+1, new_front_nodes_id, new_front_faces_id);



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
Front_3D
AeroExtrusion_3D::Compute1stLayer(Front_3D A_Front_IN, Variable<double>* A_distance, Variable<math::Vector3d>* A_vectors)
{
	std::cout << "---------> build layer: " << A_Front_IN.getFrontID()+1 << std::endl;

	std::map<TCellID, TCellID> map_new_nodes = ComputeIdealPositions(A_Front_IN, m_params_aero.delta_cl, A_distance, A_vectors);

	// Init the multiple_info type
	for (auto f_id:A_Front_IN.getFaces())
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

	// Mise à jour de l'indice de couche
	Variable<int>* couche_id = m_meshH->getVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	for (auto n_id:A_Front_IN.getNodes()){
		couche_id->set(map_new_nodes[n_id], A_Front_IN.getFrontID()+1);
	}

	std::vector<TCellID> front_nodes = A_Front_IN.getNodes();
	std::vector<TCellID> front_faces = A_Front_IN.getFaces();

	// Ajout des hex restants
	for (auto f_id:front_faces){
		Face f = m_meshH->get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		TemplateFace(f_id, A_Front_IN, map_new_nodes);
	}

	// Supression des noeuds non utilisés
	math::Utils::MeshCleaner(m_meshH);

	//----------------------//
	// Initialisation du 	//
	// front de sortie		//
	//----------------------//
	std::vector<TCellID> new_front_nodes_id;
	std::vector<TCellID> new_front_faces_id;

	Variable<int>* var_node_couche_id = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* var_face_couche_id = m_meshH->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");

	for (auto n_id:m_meshH->nodes())
	{
		if (var_node_couche_id->value(n_id) == A_Front_IN.getFrontID()+1)
		{
			new_front_nodes_id.push_back(n_id);
		}
	}
	for (auto f_id:m_meshH->faces())
	{
		if (var_face_couche_id->value(f_id) == A_Front_IN.getFrontID()+1)
		{
			new_front_faces_id.push_back(f_id);
		}
	}
	Front_3D Front_OUT = Front_3D(A_Front_IN.getFrontID()+1, new_front_nodes_id, new_front_faces_id);

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
int
AeroExtrusion_3D::SingleEdgeClassification(TCellID e_id, Front_3D &Front)
{
	int edge_classification(0);

	Variable<int>* var_node_couche_id = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* var_face_couche_id = m_meshH->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");

	Edge e = m_meshH->get<Edge>(e_id);

	std::vector<Face> e_faces = e.get<Face>();
	std::vector<Face> e_faces_on_front;
	for (auto f:e_faces)
	{
		if (var_face_couche_id->value(f.id()) == Front.getFrontID())
		{
			e_faces_on_front.push_back(f);
		}
	}
	// Each edge of the front is connected to two faces of the front.
	// So the size of the vector e_faces_on_front is 2.

	math::Vector3d n1 = e_faces_on_front[0].normal();
	math::Vector3d n2 = e_faces_on_front[1].normal();

	// Each face of the front is connected to a unique region
	// So the size of the vector e_faces_on_front[0].get<Region>() is 1.
	Region r1 = (e_faces_on_front[0].get<Region>())[0];
	Region r2 = (e_faces_on_front[1].get<Region>())[0];

	math::Point r1_center = r1.center();
	math::Point r2_center = r2.center();

	math::Point f1_center = e_faces_on_front[0].center();
	math::Point f2_center = e_faces_on_front[1].center();

	math::Vector3d v1 = (f1_center-r1.center()).normalize() ;
	math::Vector3d v2 = (f2_center-r2.center()).normalize() ;

	// Compute if the vectors n1 and n2 are well oriented.
	if ( n1.dot(v1) < 0 )
	{
		n1=-n1;
	}
	if ( n2.dot(v2) < 0)
	{
		n2=-n2;
	}

	/*
	double angle(0);
	std::vector<Region> e_regions = e.get<Region>() ;
	for (auto r:e_regions)
	{
		std::vector<Face> two_adj_faces = math::Utils::getFacesAdjToEdgeInHexa(m_meshH, e.id(), r.id()) ;
		Edge e_opp_0 = math::Utils::oppositeEdgeInFace(m_meshH, e.id(), two_adj_faces[0].id());
		Edge e_opp_1 = math::Utils::oppositeEdgeInFace(m_meshH, e.id(), two_adj_faces[1].id());

		math::Point center_e = e.center();
		math::Vector3d v0 = (e_opp_0.center()-center_e).normalize() ;
		math::Vector3d v1 = (e_opp_1.center()-center_e).normalize() ;

		angle += acos(v0.dot(v1));
	}
	 */

	// Edge classification.
	// 0 : Side
	// 1 : Corner
	// 2 : End
	// 3 : Reversal

	double angle = acos(n1.dot(n2));

	if (M_PI/4.0 <= angle && angle < 3.0*M_PI/4.0)
	{
		edge_classification = 1;
	}
	else if (5.0*M_PI/4.0 <= angle && angle < 7.0*M_PI/4.0)
	{
		edge_classification = 2;
	}
	else if (3.0*M_PI/4.0 <= angle && angle < 5.0*M_PI/4.0)
	{
		edge_classification = 3;
	}

	return edge_classification;
}
/*------------------------------------------------------------------------*/
Variable<int>*
AeroExtrusion_3D::FrontEdgesClassification(Front_3D &Front)
{
	Variable<int>* var_front_edges_classification = m_meshH->getOrCreateVariable<int, GMDS_EDGE>("Edges_Classification");
	Variable<int>* var_node_couche_id = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* var_face_couche_id = m_meshH->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");

	for (auto e_id:m_meshH->edges())
	{
		Edge e = m_meshH->get<Edge>(e_id);
		std::vector<Node> e_nodes = e.get<Node>() ;
		if (var_node_couche_id->value(e_nodes[0].id()) == Front.getFrontID()
		    && var_node_couche_id->value(e_nodes[1].id()) == Front.getFrontID())
		{
			// So, the edge is on the front.
			int edge_classification = SingleEdgeClassification(e_id, Front);
			var_front_edges_classification->set(e_id, edge_classification);
			if (edge_classification==1)
			{
				// Init Corner Edge information
				m_EdgeInfo[e_id].singularity_type = 1;
				m_EdgeInfo[e_id].CORNER_n_face_created[e_nodes[0].id()] = false;
				m_EdgeInfo[e_id].CORNER_n_face_created[e_nodes[1].id()] = false;
			}
		}
	}

	return var_front_edges_classification;

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

		if ( (compteur_corner == 0 && compteur_end == 0 && compteur_reversal == 0)
		    || (compteur_corner==2 && n_ordered_edges.size()==4)
		    || (compteur_end==2 && n_ordered_edges.size()==4))
		{
			// The node is regular
		}
		else if (compteur_corner == 3 && n_ordered_edges.size()==3)
		{
			sing_nodes[n_id] = 1;
		}
		else if (compteur_end == 3 && n_ordered_edges.size()==3)
		{
			sing_nodes[n_id] = 4;
			std::cout << "Test" << std::endl;
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

	std::vector<TCellID> n_ordered_edges = AFront.orderedFrontEdgesAroundNode(m_meshH, n_id);

	Edge e_adj_1 = m_meshH->get<Edge>(n_ordered_edges[0]);
	Edge e_adj_2 = m_meshH->get<Edge>(n_ordered_edges[1]);
	Edge e_adj_3 = m_meshH->get<Edge>(n_ordered_edges[2]);

	Node n = m_meshH->get<Node>(n_id);
	Node n6 = m_meshH->get<Node>(map_new_nodes[n_id]);

	Node n_adj_1 = e_adj_1.getOppositeNode(n);
	Node n_adj_2 = e_adj_2.getOppositeNode(n);
	Node n_adj_3 = e_adj_3.getOppositeNode(n);

	math::Vector3d v1 = n_adj_1.point()-n.point() ;
	math::Vector3d v2 = n_adj_2.point()-n.point() ;
	math::Vector3d v3 = n_adj_3.point()-n.point() ;

	/*
	 * TEST 1
	math::Point p1 = n.point() + -v1 ;
	math::Point p3 = n.point() + -v2 ;
	math::Point p4 = n.point() + -v3 ;
	math::Point p2 = n.point() + (-v1-v2) ;
	math::Point p5 = n.point() + (p4-n.point()) + (p1-n.point()) ;
	math::Point p7 = n.point() + (p4-n.point()) + (p3-n.point()) ;
	*/

	/*
	 * TEST 2
	double dist_new = (1.0/m_params_aero.nbr_couches)*(AFront.getFrontID() + 1.0/sqrt(3.0) ) ;
	std::cout << "Dist new: " << dist_new << std::endl;
	math::Point p1 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dist_new, m_DistanceField, -v1);
	math::Point p3 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dist_new, m_DistanceField, -v2);
	math::Point p4 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dist_new, m_DistanceField, -v3);
	math::Point p2 = n.point() + (p3-n.point()) + (p1-n.point()) ;
	math::Point p5 = p4 + (p1-n.point()) ;
	math::Point p7 = p4 + (p3-n.point()) ;
	 */

	// TEST 3
	/*
	double diag = (n.point()-n6.point()).norm() ;
	double cote = diag/sqrt(3.0);
	math::Point p1 = n.point() + cote*(-v1.normalize()) ;
	math::Point p3 = n.point() + cote*(-v2.normalize()) ;
	math::Point p4 = n.point() + cote*(-v3.normalize()) ;
	math::Point p2 = n.point() + (p3-n.point()) + (p1-n.point()) ;
	math::Point p5 = p4 + (p1-n.point()) ;
	math::Point p7 = p4 + (p3-n.point()) ;
	 */
	v1.normalize();
	v2.normalize();
	v3.normalize();
	math::Point p1 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, -v1);
	math::Point p3 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, -v2);
	math::Point p4 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, -v3);
	math::Point p2 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, -v1-v2); // (p3-n.point()) + (p1-n.point()) ;
	math::Point p5 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, -v1-v3); // p4 + (p1-n.point()) ;
	math::Point p7 = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, n.point(), dc, m_DistanceField, -v2-v3); //p4 + (p3-n.point()) ;


	//math::Point p6 = p4 + (p7-p4) + (p5-p4) ;
	//n6.setPoint(p6);

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
	TCellID f_1 = math::Utils::CommonFace3Nodes(m_meshH, n_id, n_adj_1.id(), n_adj_2.id());
	TCellID f_2 = math::Utils::CommonFace3Nodes(m_meshH, n_id,  n_adj_2.id(), n_adj_3.id());
	TCellID f_3 = math::Utils::CommonFace3Nodes(m_meshH, n_id,  n_adj_3.id(), n_adj_1.id());

	m_FaceInfo[f_1].next_nodes[n_id] = n4.id() ;
	m_FaceInfo[f_2].next_nodes[n_id] = n1.id() ;
	m_FaceInfo[f_3].next_nodes[n_id] = n3.id() ;

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

	std::pair<TCellID, TCellID> pair_1(f_1, n.id());
	std::pair<TCellID, TCellID> pair_2(f_3, n.id());
	std::pair<TCellID, TCellID> pair_3(f_2, n.id());

	m_EdgeInfo[e1_id].CORNER_next_nodes[pair_1] = n4.id();
	m_EdgeInfo[e1_id].CORNER_next_nodes[pair_2] = n3.id();

	m_EdgeInfo[e2_id].CORNER_next_nodes[pair_1] = n4.id();
	m_EdgeInfo[e2_id].CORNER_next_nodes[pair_3] = n1.id();

	m_EdgeInfo[e3_id].CORNER_next_nodes[pair_3] = n1.id();
	m_EdgeInfo[e3_id].CORNER_next_nodes[pair_2] = n3.id();

	m_EdgeInfo[e1_id].CORNER_diag_next_node[n.id()] = n7.id() ;
	m_EdgeInfo[e2_id].CORNER_diag_next_node[n.id()] = n5.id() ;
	m_EdgeInfo[e3_id].CORNER_diag_next_node[n.id()] = n2.id() ;

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
AeroExtrusion_3D::TemplateNode3End(Front_3D &AFront, TCellID n_id, int mark_edgesTreated, int mark_facesTreated)
{
	//std::cout << "Template Node 3 End au noeud " << n_id << std::endl;
	TCellID r_id;

	Node n = m_meshH->get<Node>(n_id);

	std::vector<TCellID> n_ordered_edges = AFront.orderedFrontEdgesAroundNode(m_meshH, n_id);

	Edge e_0 = m_meshH->get<Edge>(n_ordered_edges[0]);
	Edge e_1 = m_meshH->get<Edge>(n_ordered_edges[1]);
	Edge e_2 = m_meshH->get<Edge>(n_ordered_edges[2]);

	Node n_e0_opp = e_0.getOppositeNode(n);
	Node n_e1_opp = e_0.getOppositeNode(n);
	Node n_e2_opp = e_0.getOppositeNode(n);

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
		n0_2_id = m_EdgeInfo[e_id].CORNER_diag_next_node[n0_id];
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

		/*
		 * TEST 1
		Node n0_1_new = m_meshH->newNode(p_n0 + f0_normal );
		Node n0_2_new = m_meshH->newNode(p_n0 + f0_normal + f1_normal );
		Node n0_3_new = m_meshH->newNode(p_n0 + f1_normal );
		*/

		// TEST 2
		/*
		double dist_new = (1.0/m_params_aero.nbr_couches)*(AFront.getFrontID() + 1.0/sqrt(3.0) ) ;
		math::Point p0_1_new  = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, p_n0, dist_new, m_DistanceField, f0_normal);
		math::Point p0_3_new  = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, p_n0, dist_new, m_DistanceField, f1_normal);

		Node n0_1_new = m_meshH->newNode(p0_1_new);
		Node n0_3_new = m_meshH->newNode(p0_3_new);
		Node n0_2_new = m_meshH->newNode(p_n0 + (p0_1_new-p_n0) + (p0_3_new-p_n0));
		*/

		// TEST 3
		/*
		n0_2_id = m_FaceInfo[e_front_faces[0]].next_ideal_nodes[e_nodes[0].id()];
		Node n0_2_new = m_meshH->get<Node>(n0_2_id);
		math::Point p0_2_new = n0_2_new.point() ;

		double diag = (e_nodes[0].point()-p0_2_new).norm() ;
		double cote = diag/sqrt(3.0);

		Node n0_1_new = m_meshH->newNode(p_n0 + cote*f0_normal.normalize() );
		Node n0_3_new = m_meshH->newNode(p_n0 + cote*f1_normal.normalize() );
		 */

		// TEST 4
		n0_2_id = m_FaceInfo[e_front_faces[0]].next_ideal_nodes[e_nodes[0].id()];
		Node n0_2_new = m_meshH->get<Node>(n0_2_id);
		math::Point p0_2_new = n0_2_new.point() ;
		math::Point p0_1_new = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, e_nodes[0].point(), dc, m_DistanceField, f0_normal.normalize());
		math::Point p0_3_new = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, e_nodes[0].point(), dc, m_DistanceField, f1_normal.normalize());

		Node n0_1_new = m_meshH->newNode(p0_1_new);
		Node n0_3_new = m_meshH->newNode(p0_3_new);

		//Node n0_1_new = m_meshH->newNode(p0_1_new);
		//Node n0_2_new = m_meshH->newNode(p0_2_new);
		//Node n0_3_new = m_meshH->newNode(p0_3_new);

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
				m_EdgeInfo[e_adj_id].CORNER_diag_next_node[n0_id] = n0_2_id;

				TCellID f_adj_0 = AFront.adjacentFaceOnFront(m_meshH, e_front_faces[0], e_adj_id);
				TCellID f_adj_1 = AFront.adjacentFaceOnFront(m_meshH, e_front_faces[1], e_adj_id);

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
		n1_2_id = m_EdgeInfo[e_id].CORNER_diag_next_node[n1_id];
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

		/*
		// TEST 1
		Node n1_1_new = m_meshH->newNode(p_n1 + f0_normal );
		Node n1_2_new = m_meshH->newNode(p_n1 + f0_normal + f1_normal );
		Node n1_3_new = m_meshH->newNode(p_n1 + f1_normal );
		 */

		/*
		// TEST 2
		double dist_new = (1.0/m_params_aero.nbr_couches)*(AFront.getFrontID() + 1.0/sqrt(3.0) ) ;
		math::Point p1_1_new  = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, p_n1, dist_new, m_DistanceField, f0_normal);
		math::Point p1_3_new  = math::Utils::AdvectedPointRK4_UniqVector_3D(m_meshT, &m_fl, p_n1, dist_new, m_DistanceField, f1_normal);

		Node n1_1_new = m_meshH->newNode(p1_1_new);
		Node n1_3_new = m_meshH->newNode(p1_3_new);
		Node n1_2_new = m_meshH->newNode(p_n1 + (p1_1_new-p_n1) + (p1_3_new-p_n1));
		 */

		// TEST 3
		/*
		n1_2_id = m_FaceInfo[e_front_faces[0]].next_ideal_nodes[e_nodes[1].id()];
		Node n1_2_new = m_meshH->get<Node>(n1_2_id);
		math::Point p1_2_new = n1_2_new.point() ;

		double diag = (e_nodes[1].point()-p1_2_new).norm() ;
		double cote = diag/sqrt(3.0);

		Node n1_1_new = m_meshH->newNode(p_n1 + cote*f0_normal.normalize() );
		Node n1_3_new = m_meshH->newNode(p_n1 + cote*f1_normal.normalize() );
		*/

		// TEST 4
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
				m_EdgeInfo[e_adj_id].CORNER_diag_next_node[n1_id] = n1_2_id;

				TCellID f_adj_0 = AFront.adjacentFaceOnFront(m_meshH, e_front_faces[0], e_adj_id);
				TCellID f_adj_1 = AFront.adjacentFaceOnFront(m_meshH, e_front_faces[1], e_adj_id);

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

	//std::cout << "t1" << std::endl;
	//std::cout << "n0: " << n0_id << ", n1: " << n1_id << std::endl;
	// Create the hexa
	r_id = math::Utils::CreateHexaNConnectivities(m_meshH, e_nodes[0], n0_1, n0_2, n0_3, e_nodes[1], n1_1, n1_2, n1_3);
	//std::cout << "t2" << std::endl;

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
void
AeroExtrusion_3D::TemplateFace(TCellID f_id, Front_3D &Front_IN, std::map<TCellID, TCellID> map_new_nodes)
{
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
std::map<int, std::vector<Edge>>
AeroExtrusion_3D::ComputeGlobalFeatureEdges(Front_3D &Front_IN, std::map<TCellID, int> singular_nodes, std::map<TCellID, int> singular_edges)
{
	std::map<int, std::vector<Edge>> feature_edges;

	int mark_isEdgeUsed = m_meshH->newMark<Edge>();

	for (auto sing_node:singular_nodes)
	{
		Node n = m_meshH->get<Node>(sing_node.first);
		std::vector<Edge> local_feature_edges;
		for (auto sing_edge:singular_edges)
		{
			Edge e_loc = m_meshH->get<Edge>(sing_edge.first) ;
			std::vector<Node> e_loc_nodes = e_loc.get<Node>();
			if (e_loc_nodes[0].id() == n.id() || e_loc_nodes[1].id() == n.id())
			{
				local_feature_edges.push_back(e_loc);
			}
		}

		for (auto e_loc_feature:local_feature_edges)
		{
			if(!m_meshH->isMarked(e_loc_feature, mark_isEdgeUsed))
			{
				std::vector<Edge> newGlobalFeature;
				Node n_loc = m_meshH->get<Node>(n.id());
				/*
				{
					newGlobalFeature.push_back(e_loc_feature);
					Node n_opp = e_loc_feature.getOppositeNode(n_loc) ;
				 	// Get next feature edge
				   std::vector<Edge> local_feature_edges_around_n_opp;
					for (auto sing_edge:singular_edges)
					{
						Edge e_loc = m_meshH->get<Edge>(sing_edge.first) ;
						std::vector<Node> e_loc_nodes = e_loc.get<Node>();
						if ( (e_loc_nodes[0].id() == n_opp.id() && e_loc_nodes[1].id() != n_loc.id() )
						    || (e_loc_nodes[1].id() == n_opp.id() && e_loc_nodes[0].id() != n_loc.id()) )
						{
							local_feature_edges.push_back(e_loc);
						}
					}
					if (local_feature_edges_around_n_opp.size() > 0)
					{
						Node n_loc = n_opp;
					}
				}
				*/
			}
		}
	}

	m_meshH->unmarkAll<Edge>(mark_isEdgeUsed);
	m_meshH->freeMark<Edge>(mark_isEdgeUsed);

	return feature_edges;
}
/*------------------------------------------------------------------------*/