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

	std::map<TCellID, int> test = getSingularNodes(Front_IN, var_front_edges_classification);

	// Ajout des hex restants
	for (auto f_id:front_faces){
		Face f = m_meshH->get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		CreateNormalHexa(f_id, Front_IN, map_new_nodes);
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
		CreateNormalHexa(f_id, A_Front_IN, map_new_nodes);
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

	return Front_OUT;

}
/*------------------------------------------------------------------------*/
void
AeroExtrusion_3D::CreateNormalHexa(TCellID f_id, Front_3D &Front_IN, std::map<TCellID, TCellID> map_new_nodes)
{
	Variable<int>* var_node_couche_id = m_meshH->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	Variable<int>* var_face_couche_id = m_meshH->getOrCreateVariable<int, GMDS_FACE>("GMDS_FACE_Couche_Id");

	Face f = m_meshH->get<Face>(f_id);
	std::vector<Node> nodes = f.get<Node>();

	/*
	TCellID n0_id = map_new_nodes[nodes[0].id()] ;
	TCellID n1_id = map_new_nodes[nodes[1].id()] ;
	TCellID n2_id = map_new_nodes[nodes[2].id()] ;
	TCellID n3_id = map_new_nodes[nodes[3].id()] ;
	 */
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
			//std::cout << "Edge: " << e_id << ", Classified: " << edge_classification << std::endl;
			var_front_edges_classification->set(e_id, edge_classification);
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
		//std::cout << "-----------------" << std::endl;
		//std::cout << "Node " << n_id << std::endl;
		for (auto e_id:n_ordered_edges)
		{
			Edge e = m_meshH->get<Edge>(e_id);
			//std::cout << "edge " << e.get<Node>()[0].id() << ", " << e.get<Node>()[1].id() << std::endl;
		}
	}
	return sing_nodes;
}
/*------------------------------------------------------------------------*/