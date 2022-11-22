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

AeroExtrusion_3D::AeroExtrusion_3D(Mesh *AMeshT, Mesh *AMeshH, ParamsAero Aparams_aero, Variable<double>* A_DistanceField, Variable<math::Vector3d>* A_VectorField) {
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
		AdvectedPointRK4_3D advpoint(m_meshT, M, dist_cible, A_distance, A_vectors);
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

	// Mise à jour de l'indice de couche
	Variable<int>* couche_id = m_meshH->getVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	for (auto n_id:Front_IN.getNodes()){
		couche_id->set(map_new_nodes[n_id], Front_IN.getFrontID()+1);
	}

	std::vector<TCellID> front_nodes = Front_IN.getNodes();
	std::vector<TCellID> front_faces = Front_IN.getFaces();


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

	TCellID n0_id = map_new_nodes[nodes[0].id()] ;
	TCellID n1_id = map_new_nodes[nodes[1].id()] ;
	TCellID n2_id = map_new_nodes[nodes[2].id()] ;
	TCellID n3_id = map_new_nodes[nodes[3].id()] ;

	Node n0 = m_meshH->get<Node>(n0_id);
	Node n1 = m_meshH->get<Node>(n1_id);
	Node n2 = m_meshH->get<Node>(n2_id);
	Node n3 = m_meshH->get<Node>(n3_id);
	var_node_couche_id->set(n0_id, Front_IN.getFrontID()+1);
	var_node_couche_id->set(n1_id, Front_IN.getFrontID()+1);
	var_node_couche_id->set(n2_id, Front_IN.getFrontID()+1);
	var_node_couche_id->set(n3_id, Front_IN.getFrontID()+1);

	// Create the hex associated to the face
	Region r = m_meshH->newHex(nodes[0], nodes[1], nodes[2], nodes[3], n0, n1, n2, n3);	// R->N (x8)
	nodes[0].add<Region>(r);
	nodes[1].add<Region>(r);
	nodes[2].add<Region>(r);
	nodes[3].add<Region>(r);
	n0.add<Region>(r);
	n1.add<Region>(r);
	n2.add<Region>(r);
	n3.add<Region>(r);

	// Create the face on the new front
	Face f_new_front = m_meshH->newQuad(n0, n1, n2, n3);	// F->N (x4)
	n0.add<Face>(f_new_front);
	n1.add<Face>(f_new_front);
	n2.add<Face>(f_new_front);
	n3.add<Face>(f_new_front);
	var_face_couche_id->set(f_new_front.id(), Front_IN.getFrontID()+1);

	// The faces between the two layers are not created.

	TCellID f0_id = math::Utils::CommonFace(m_meshH, nodes[0].id(), nodes[1].id(), n0_id, n1_id);
	if (f0_id == NullID)		// Then, the face doesn't exist yet
	{
		Face f0 = m_meshH->newQuad(nodes[0], nodes[1], n1, n0);	// F->N (x4)
		nodes[0].add<Face>(f0);	// N->F
		nodes[1].add<Face>(f0);	// N->F
		n0.add<Face>(f0);			// N->F
		n1.add<Face>(f0);			// N->F
	}

	TCellID f1_id = math::Utils::CommonFace(m_meshH, nodes[1].id(), nodes[2].id(), n1_id, n2_id);
	if (f1_id == NullID)		// Then, the face doesn't exist yet
	{
		Face f1 = m_meshH->newQuad(nodes[1], nodes[2], n2, n1);	// F->N (x4)
		nodes[1].add<Face>(f1);	// N->F
		nodes[2].add<Face>(f1);	// N->F
		n1.add<Face>(f1);			// N->F
		n2.add<Face>(f1);			// N->F
	}

	TCellID f2_id = math::Utils::CommonFace(m_meshH, nodes[2].id(), nodes[3].id(), n2_id, n3_id);
	if (f2_id == NullID)		// Then, the face doesn't exist yet
	{
		Face f2 = m_meshH->newQuad(nodes[2], nodes[3], n3, n2);	// F->N (x4)
		nodes[2].add<Face>(f2);	// N->F
		nodes[3].add<Face>(f2);	// N->F
		n2.add<Face>(f2);			// N->F
		n3.add<Face>(f2);			// N->F
	}

	TCellID f3_id = math::Utils::CommonFace(m_meshH, nodes[0].id(), nodes[3].id(), n0_id, n3_id);
	if (f3_id == NullID)		// Then, the face doesn't exist yet
	{
		Face f3 = m_meshH->newQuad(nodes[0], nodes[3], n3, n0);	// F->N (x4)
		nodes[0].add<Face>(f3);	// N->F
		nodes[3].add<Face>(f3);	// N->F
		n0.add<Face>(f3);			// N->F
		n3.add<Face>(f3);			// N->F
	}

	/*
	TCellID e0_id = math::Utils::CommonEdge(m_meshQ, nodes[0].id(), n0.id());
	TCellID e1_id = math::Utils::CommonEdge(m_meshQ, nodes[1].id(), n1.id());
	Edge e0, e1;
	if ( e0_id == NullID ){
		// Si l'arête n'existe pas, on la créé et on initialise les connectivités N->E
		e0 = m_meshQ->newEdge(nodes[0].id(), n0.id());		// E->N (x2)
		nodes[0].add<Edge>(e0);												// N->E
		n0.add<Edge>(e0);														// N->E
	}
	else{
		e0 = m_meshQ->get<Edge>(e0_id);
	}

	// Idem pour l'arête e1
	if ( e1_id == NullID ){
		e1 = m_meshQ->newEdge(nodes[1].id(), n1.id());		// E->N (x2)
		nodes[1].add<Edge>(e1);												// N->E
		n1.add<Edge>(e1);														// N->E
	}
	else{
		e1 = m_meshQ->get<Edge>(e1_id);
	}

	// Connectivités F->E
	f.add<Edge>(e);		// F->E
	f.add<Edge>(e0);		// F->E
	f.add<Edge>(e1);		// F->E
	f.add<Edge>(e_opp);	// F->E

	// Connectivités E->F
	e.add<Face>(f);		// E->F
	e0.add<Face>(f);		// E->F
	e1.add<Face>(f);		// E->F
	e_opp.add<Face>(f);	// E->F
	 */

}
/*------------------------------------------------------------------------*/