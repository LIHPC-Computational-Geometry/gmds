//
// Created by rochec on 14/04/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/claire/Utils.h>
#include <gmds/claire/AeroException.h>
#include <gmds/claire/AeroExtrusion_2D.h>
#include <gmds/claire/AdvectedPointRK4_2D.h>
#include <gmds/claire/AeroMeshQuality.h>
#include <gmds/claire/SmoothingPaving_2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AeroExtrusion_2D::AeroExtrusion_2D(Mesh *AMeshT, Mesh *AMeshQ, ParamsAero Aparams_aero, Variable<math::Vector3d>* A_VectorField) :
	m_meshT(AMeshT),
  	m_meshQ(AMeshQ),
  	m_fl(m_meshT)
{
	m_params_aero = Aparams_aero;
	m_VectorField = A_VectorField;
}


/*------------------------------------------------------------------------*/
AeroExtrusion_2D::STATUS
AeroExtrusion_2D::execute()
{
	// Exemple exception
	//if(m_mesh==NULL)
	//	throw AeroException("ERROR: Invalid mesh pointer");

	Front Current_Front = Compute1stLayer(m_meshT->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"), m_params_aero.delta_cl,
	                m_VectorField);
	/*
	Current_Front = ComputeLayer(Current_Front, m_meshT->getVariable<double,GMDS_NODE>("GMDS_Distance_2"), 0.25,
	                             m_meshT->getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	Current_Front = ComputeLayer(Current_Front, m_meshT->getVariable<double,GMDS_NODE>("GMDS_Distance_2"), 0.5,
	                             m_meshT->getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	Current_Front = ComputeLayer(Current_Front, m_meshT->getVariable<double,GMDS_NODE>("GMDS_Distance_2"), 0.75,
	                             m_meshT->getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	Current_Front = ComputeLayer(Current_Front, m_meshT->getVariable<double,GMDS_NODE>("GMDS_Distance_2"), 1,
	                             m_meshT->getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	                             */

	double pas_couche = 1.0/m_params_aero.nbr_couches ;

	for (int i=2; i <= m_params_aero.nbr_couches; i++)
	{
		Current_Front = ComputeLayer(Current_Front, m_meshT->getVariable<double,GMDS_NODE>("GMDS_Distance"), i*pas_couche,
	                             m_VectorField);
		std::cout << "distance : " << i*pas_couche << std::endl;
	}

	return AeroExtrusion_2D::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
std::map<TCellID, TCellID>
   AeroExtrusion_2D::ComputeIdealPositions(Front AFront, double dist_cible, Variable<double>* A_distance, Variable<math::Vector3d>* A_vectors)
{
	std::map<TCellID, TCellID> map_optnexpoint;
	std::vector<TCellID> front_nodes = AFront.getNodes();

	for (auto n_id:front_nodes){
		Node n = m_meshQ->get<Node>(n_id);
		math::Point M = n.point();
		AdvectedPointRK4_2D advpoint(m_meshT, &m_fl, M, dist_cible, A_distance, A_vectors);
		advpoint.execute();
		//math::Point P = advpoint.getPend();
		Node n_new = m_meshQ->newNode(advpoint.getPend());
		map_optnexpoint[n_id] = n_new.id() ;
		}

	return map_optnexpoint;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
Front
AeroExtrusion_2D::Compute1stLayer(Variable<double>* A_distance, double dist_cible, Variable<math::Vector3d>* A_vectors){

	std::cout << "--------- couche " << 1 << std::endl;

	Front Front_Paroi;
	Front_Paroi.setFrontID(0);
	Front_Paroi.initializeFromLayerId(m_meshQ, 0);

	//Front First_Front;
	//First_Front.setFrontID(Front_Paroi.getFrontID()+1);

	std::map<TCellID, TCellID> map_new_nodes = ComputeIdealPositions(Front_Paroi, dist_cible, A_distance, A_vectors);

	// Mise à jour de l'indice de couche
	Variable<int>* couche_id = m_meshQ->getVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	for (auto n_id:Front_Paroi.getNodes()){
		couche_id->set(map_new_nodes[n_id], Front_Paroi.getFrontID()+1);
	}

	Front_Paroi.initializeNodeType(m_meshQ, map_new_nodes);

	std::vector<TCellID> front_nodes = Front_Paroi.getNodes();
	std::vector<TCellID> front_edges = Front_Paroi.getEdges();

	/*
	// Remarque : on pourait initialiser les arêtes e0 et e1 qui relient les deux couches dans cette boucle.
	for (auto n_id:front_nodes){
		// Ajout de la nouvelle arête du front dans l'objet concerné
		First_Front.addNodeId(map_new_nodes[n_id]);
		couche_id->set(map_new_nodes[n_id], First_Front.getFrontID());
	}
	 */

	/*
	for (auto e_id:front_edges){

		Edge e = m_meshQ->get<Edge>(e_id);
		std::vector<Node> nodes = e.get<Node>();

		Node n0 = m_meshQ->get<Node>(map_new_nodes[nodes[0].id()]);
		Node n1 = m_meshQ->get<Node>(map_new_nodes[nodes[1].id()]);


		// On créé la face associée à l'arête
		Face f = m_meshQ->newQuad(nodes[0], nodes[1], n1, n0) ; // F->N (x4)
		nodes[0].add<Face>(f);	// N->F
		nodes[1].add<Face>(f);	// N->F
		n0.add<Face>(f);			// N->F
		n1.add<Face>(f);			// N->F

		// Création de l'arête sur le nouveau front
		Edge e_opp = m_meshQ->newEdge(n0, n1);	// E->N (x2)
		n0.add<Edge>(e_opp);										// N->E
		n1.add<Edge>(e_opp);										// N->E

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


		// Ajout de la nouvelle arête du front dans l'objet concerné
		First_Front.addEdgeId(e_opp.id());

	}
	 */

	for (auto n_id:front_nodes){

		std::vector<TCellID> neighbor_nodes_id = Front_Paroi.getNeighbors(n_id) ;
		Node n_neighbor_1 = m_meshQ->get<Node>(neighbor_nodes_id[0]);
		Node n_neighbor_2 = m_meshQ->get<Node>(neighbor_nodes_id[1]);
		Node n = m_meshQ->get<Node>(n_id);

		math::Vector3d v1 = (n_neighbor_1.point() - n.point()).normalize();
		math::Vector3d v2 = (n_neighbor_2.point() - n.point()).normalize();

		double angle = (acos(v1.dot(v2))*180/M_PI) ;

		if (abs(angle) < 40)
		{
			std::cout << "Angle : " << angle << std::endl;
			//Insertion_Double(Front_Paroi, n_id, A_distance, dist_cible, A_vectors);
		}

	}

	// Ajout des quad restants
	for (auto e_id:front_edges){
		Edge e = m_meshQ->get<Edge>(e_id);
		std::vector<Node> nodes = e.get<Node>();
			CreateNormalQuad(e_id, Front_Paroi);
	}

	// Initialisation du front de sortie
	Front First_Front;
	First_Front.initializeFromLayerId(m_meshQ, Front_Paroi.getFrontID()+1);

	return First_Front;

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
Front
AeroExtrusion_2D::ComputeLayer(Front Front_IN, Variable<double>* A_distance, double dist_cible, Variable<math::Vector3d>* A_vectors){

	std::cout << "--------- couche " << Front_IN.getFrontID()+1 << std::endl;

	std::map<TCellID, TCellID> map_new_nodes = ComputeIdealPositions(Front_IN, dist_cible, A_distance, A_vectors);

	// Mise à jour de l'indice de couche
	Variable<int>* couche_id = m_meshQ->getVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	for (auto n_id:Front_IN.getNodes()){
		couche_id->set(map_new_nodes[n_id], Front_IN.getFrontID()+1);
	}

	Front_IN.initializeNodeType(m_meshQ, map_new_nodes);

	std::vector<TCellID> front_nodes = Front_IN.getNodes();
	std::vector<TCellID> front_edges = Front_IN.getEdges();

	Front_IN.setMultorFusFromLimits(m_meshQ, m_params_aero.x_lim, m_params_aero.y_lim, m_params_aero.z_lim);

	// Boucle d'insertion d'éléments
	bool do_smooth(true);

	if (Front_IN.getFrontID()+1 != m_params_aero.nbr_couches ) {
		while (do_smooth) {
			do_smooth = false;
			// A FAIRE : lissage de la couche

			TCellID node_id;
			int type_node(0);
			getSingularNode(Front_IN, node_id, type_node);

			if (type_node == 1) {
				// Insertion
				std::cout << "INSERTION QUAD AU NOEUD : " << node_id << std::endl;
				Insertion(Front_IN, node_id, A_distance, dist_cible, A_vectors);
				do_smooth = true;
			}
			else if (type_node == 2) {
				// Fusion
				std::cout << "FUSION QUAD AU NOEUD : " << node_id << std::endl;
				Fusion(Front_IN, node_id);
				do_smooth = true;
			}
		}
	}

	// Ajout des quad restants
	for (auto e_id:front_edges){
		Edge e = m_meshQ->get<Edge>(e_id);
		std::vector<Node> nodes = e.get<Node>();
		if (Front_IN.getNodeType(nodes[0].id()) != 2
		    && Front_IN.getNodeType(nodes[1].id()) != 2 ){
			CreateNormalQuad(e_id, Front_IN);
		}

	}

	// Supression des noeuds non utilisés
	math::Utils::MeshCleaner(m_meshQ);

	// Initialisation du front de sortie
	Front Front_OUT;
	Front_OUT.initializeFromLayerId(m_meshQ, Front_IN.getFrontID()+1);

	// Lissage de la couche
	/*
	if ( abs(dist_cible - 1.0) > pow(10,-6)) {
		SmoothingPaving_2D smoother(m_meshQ, Front_OUT);
		smoother.execute();
	}
	 */

	return Front_OUT;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroExtrusion_2D::getSingularNode(Front Front_IN, TCellID &node_id, int &type){

	std::vector<TCellID> front_nodes = Front_IN.getNodes();
	std::vector<TCellID> front_edges = Front_IN.getEdges();

	bool singu_not_found(true);

	node_id = NullID;
	type = 0;

	for (auto n_id:front_nodes) {

		Node n = m_meshQ->get<Node>(n_id);
		std::vector<TCellID> neighbors_nodes = Front_IN.getNeighbors(n_id);

		std::vector<Node> nodes_quad_1;
		std::vector<Node> nodes_quad_2;

		nodes_quad_1.push_back(m_meshQ->get<Node>(n_id));
		nodes_quad_1.push_back(m_meshQ->get<Node>(neighbors_nodes[0]));
		nodes_quad_1.push_back(m_meshQ->get<Node>(Front_IN.getNextNode(neighbors_nodes[0], n_id)));
		nodes_quad_1.push_back(m_meshQ->get<Node>(Front_IN.getNextNode(n_id, neighbors_nodes[0])));

		nodes_quad_2.push_back(m_meshQ->get<Node>(n_id));
		nodes_quad_2.push_back(m_meshQ->get<Node>(neighbors_nodes[1]));
		nodes_quad_2.push_back(m_meshQ->get<Node>(Front_IN.getNextNode(neighbors_nodes[1], n_id)));
		nodes_quad_2.push_back(m_meshQ->get<Node>(Front_IN.getNextNode(n_id, neighbors_nodes[1])));

		// Les tests pour la FUSION
		/*
		if (Front_IN.isFusionable(n_id)) {

			// test arête écrasée premier quad
			double min_lenght_1 =
			   math::AeroMeshQuality::minlenghtedge(nodes_quad_1[0].point(), nodes_quad_1[1].point(), nodes_quad_1[2].point(), nodes_quad_1[3].point());

			// test arête écrasée second quad
			double min_lenght_2 =
			   math::AeroMeshQuality::minlenghtedge(nodes_quad_2[0].point(), nodes_quad_2[1].point(), nodes_quad_2[2].point(), nodes_quad_2[3].point());

			if (singu_not_found && (min_lenght_1 < 1.0 && min_lenght_2 < 1.0)) {
				node_id = n_id;
				type = 2;
				singu_not_found = false;
			}

		}
		 */


		// Les tests pour l'INSERSION
		if (singu_not_found && Front_IN.isMultiplicable(n_id)) {

			// Test angle ouverture

			Node n_ideal = m_meshQ->get<Node>(Front_IN.getIdealNode(n_id));
			Node n_neighbor_1 = m_meshQ->get<Node>(neighbors_nodes[0]);
			Node n_neighbor_2 = m_meshQ->get<Node>(neighbors_nodes[1]);

			/*
			double angle_ouverture = math::AeroMeshQuality::AngleOuverture(n.point(), n_ideal.point(),
			                                                               n_neighbor_1.point(), n_neighbor_2.point());

			if (singu_not_found && (angle_ouverture > 7.0*M_PI/6.0)) {
				node_id = n_id;
				type = 1;
				singu_not_found = false;
			}
			*/


			// Alignement avec le flow
			//math::Vector3d v = (m_meshQ->get<Node>(Front_IN.getNextNode(n_id, neighbors_nodes[1])).point() - m_meshQ->get<Node>(n_id).point()) ;
			math::Vector3d v = (n_ideal.point() - n.point()) ;
			v.normalize();
			math::Vector3d v_flow({cos(m_params_aero.angle_attack*M_PI/180.0), sin(m_params_aero.angle_attack*M_PI/180.0), 0.0}) ;
			double angle = v.dot(v_flow)*(180.0/M_PI);
			double angle_on_layer = math::AeroMeshQuality::AngleOuverture(n.point(), n_ideal.point(),
			                                                               n_neighbor_1.point(), n_neighbor_2.point());
			double theta = M_PI/4.0 ;
			bool test(angle_on_layer > 3.0*M_PI/2.0);
			if (singu_not_found && angle_on_layer > 3.0*M_PI/2.0 - theta && abs(angle) < 55 && abs(angle) > 35) //
			{
				node_id = n_id;
				type = 1;
				singu_not_found = false;
			}



			/*
			double internal_angle_1 = math::AeroMeshQuality::InternalAngleDeviationQUAD(nodes_quad_1[0].point(), nodes_quad_1[1].point(),
			                                                                            nodes_quad_1[2].point(), nodes_quad_1[3].point()) ;
			double internal_angle_2 = math::AeroMeshQuality::InternalAngleDeviationQUAD(nodes_quad_2[0].point(), nodes_quad_2[1].point(),
			                                                  									nodes_quad_2[2].point(), nodes_quad_2[3].point()) ;

			   if (singu_not_found && (internal_angle_1 < 40 || internal_angle_2 < 44.5)) {
			      node_id = n_id;
			      type = 1;
			      singu_not_found = false;
			   }
         */
			}

	}


}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroExtrusion_2D::Insertion(Front &Front_IN, TCellID n_id,
                            Variable<double>* A_distance, double dist_cible, Variable<math::Vector3d>* A_vectors){

	std::vector<TCellID> neighbors_nodes = Front_IN.getNeighbors(n_id);
	Node n = m_meshQ->get<Node>(n_id);

	// Contruction du premier noeud n1
	Node n_neighbor = m_meshQ->get<Node>(neighbors_nodes[0]);
	//math::Point P1_construction = n.point() + 0.33*(n_neighbor.point()-n.point()) ;
	//AdvectedPointRK4_2D advpoint_n1(m_meshT, P1_construction, dist_cible, A_distance, A_vectors);
	//advpoint_n1.execute();
	//Node n1 = m_meshQ->newNode(advpoint_n1.getPend());

	// TEST 2
	/*
	TCellID n_next_ideal_id = Front_IN.getIdealNode(n_id);
	TCellID n_next_id = Front_IN.getNextNode(neighbors_nodes[0],n_id);
	Node n_next_ideal = m_meshQ->get<Node>(n_next_ideal_id) ;
	Node n_next = m_meshQ->get<Node>(n_next_id) ;
	math::Vector3d P1_construction = 0.33*(n_next.point()-n_next_ideal.point()) ;
	Node n1 = m_meshQ->newNode(n_next_ideal.point() + P1_construction);
	*/



	// TEST 3
	/*
	math::Vector3d v_flow({cos(m_params_aero.angle_attack*M_PI/180.0), sin(m_params_aero.angle_attack*M_PI/180.0), 0.0}) ;
	math::Vector3d v_flow_ortho({-v_flow.Y(), v_flow.X(), 0.0});
	v_flow.normalize();
	v_flow_ortho.normalize();

	math::Vector3d v10 = (n.point() - n_neighbor.point()).normalize() ;
	math::Vector3d v10_ortho({-v10.Y(), v10.X(), 0.0}) ;

	math::Vector3d v_follow;
	double min_angle(2.0*M_PI);
	if ( v10_ortho.dot(v_flow) < min_angle )
	{
		min_angle = v10_ortho.dot(v_flow);
		v_follow = v_flow;
	}
	if ( v10_ortho.dot(-v_flow) < min_angle )
	{
		min_angle = v10_ortho.dot(-v_flow);
		v_follow = -v_flow;
	}
	if ( v10_ortho.dot(v_flow_ortho) < min_angle )
	{
		min_angle = v10_ortho.dot(v_flow_ortho);
		v_follow = v_flow_ortho;
	}
	if ( v10_ortho.dot(-v_flow_ortho) < min_angle )
	{
		min_angle = v10_ortho.dot(-v_flow_ortho);
		v_follow = -v_flow_ortho;
	}

	Variable<math::Vector3d>* var_flow = m_meshT->newVariable<math::Vector3d , GMDS_NODE>("Flow") ;
	for (auto n_id:m_meshT->nodes())
	{
		var_flow->set(n_id, v_follow);
	}
	AdvectedPointRK4_2D advpoint_n1(m_meshT, n.point(), dist_cible, A_distance, var_flow);
	advpoint_n1.execute();
	Node n1 = m_meshQ->newNode(advpoint_n1.getPend());
	m_meshT->deleteVariable(GMDS_NODE, "Flow") ;
	*/



	// TEST 4
	math::Vector3d v_follow = (m_meshQ->get<Node>(Front_IN.getIdealNode(n_neighbor.id())).point() - n_neighbor.point()).normalize() ;
	v_follow += (m_meshQ->get<Node>(Front_IN.getIdealNode(n_id)).point() - n.point() ).normalize() ;
	v_follow.normalize();
	Variable<math::Vector3d>* var_flow = m_meshT->newVariable<math::Vector3d , GMDS_NODE>("Flow") ;
	for (auto n_id:m_meshT->nodes())
	{
		var_flow->set(n_id, v_follow);
	}
	AdvectedPointRK4_2D advpoint_n1(m_meshT, &m_fl, n.point(), dist_cible, A_distance, var_flow);
	advpoint_n1.execute();
	Node n1 = m_meshQ->newNode(advpoint_n1.getPend());
	m_meshT->deleteVariable(GMDS_NODE, "Flow") ;



	// Contruction du premier noeud n2
	n_neighbor = m_meshQ->get<Node>(neighbors_nodes[1]);
	/*
	P1_construction = n.point() + 0.33*(n_neighbor.point()-n.point()) ;
	AdvectedPointRK4_2D advpoint_n2(m_meshT, P1_construction, dist_cible, A_distance, A_vectors);
	advpoint_n2.execute();
	Node n2 = m_meshQ->newNode(advpoint_n2.getPend());
	 */

	// TEST 2
	/*
	n_next_ideal_id = Front_IN.getIdealNode(n_id);
	n_next_id = Front_IN.getNextNode(neighbors_nodes[1],n_id);
	n_next_ideal = m_meshQ->get<Node>(n_next_ideal_id) ;
	n_next = m_meshQ->get<Node>(n_next_id) ;
	P1_construction = 0.33*(n_next.point()-n_next_ideal.point()) ;
	Node n2 = m_meshQ->newNode(n_next_ideal.point() + P1_construction);
	*/

	// TEST 3
	/*
	math::Vector3d v20 = (n.point() - n_neighbor.point()).normalize() ;
	math::Vector3d v20_ortho({-v10.Y(), v10.X(), 0.0}) ;

	min_angle = 2.0*M_PI;
	if ( v20_ortho.dot(v_flow) < min_angle )
	{
		min_angle = v20_ortho.dot(v_flow);
		v_follow = v_flow;
	}
	if ( v20_ortho.dot(-v_flow) < min_angle )
	{
		min_angle = v20_ortho.dot(-v_flow);
		v_follow = -v_flow;
	}
	if ( v20_ortho.dot(v_flow_ortho) < min_angle )
	{
		min_angle = v20_ortho.dot(v_flow_ortho);
		v_follow = v_flow_ortho;
	}
	if ( v20_ortho.dot(-v_flow_ortho) < min_angle )
	{
		min_angle = v20_ortho.dot(-v_flow_ortho);
		v_follow = -v_flow_ortho;
	}

	var_flow = m_meshT->newVariable<math::Vector3d , GMDS_NODE>("Flow") ;
	for (auto n_id:m_meshT->nodes())
	{
		var_flow->set(n_id, v_follow);
	}
	AdvectedPointRK4_2D advpoint_n2(m_meshT, n.point(), dist_cible, A_distance, var_flow);
	advpoint_n1.execute();
	Node n2 = m_meshQ->newNode(advpoint_n1.getPend());
	m_meshT->deleteVariable(GMDS_NODE, "Flow") ;
	*/


	// TEST 4
	v_follow = m_meshQ->get<Node>(Front_IN.getIdealNode(n_neighbor.id())).point() - n_neighbor.point() ;
	var_flow = m_meshT->newVariable<math::Vector3d , GMDS_NODE>("Flow") ;
	for (auto n_id:m_meshT->nodes())
	{
		var_flow->set(n_id, v_follow);
	}
	AdvectedPointRK4_2D advpoint_n2(m_meshT, &m_fl, n.point(), dist_cible, A_distance, var_flow);
	advpoint_n2.execute();
	Node n2 = m_meshQ->newNode(advpoint_n2.getPend());
	m_meshT->deleteVariable(GMDS_NODE, "Flow") ;






	//=================================================
	//			QUAD CREATION & CONNECTIVITIES
	//=================================================


	// Création du quad
	TCellID n0_id = Front_IN.getIdealNode(n_id);
	Node n0 = m_meshQ->get<Node>(n0_id);

	// On créé la face associée à l'arête
	Face f = m_meshQ->newQuad(n, n1, n0, n2) ; // F->N (x4)
	n.add<Face>(f);		// N->F
	n0.add<Face>(f);		// N->F
	n1.add<Face>(f);		// N->F
	n2.add<Face>(f);		// N->F

	// Création des arêtes
	Edge e0 = m_meshQ->newEdge(n0, n1);	// E->N (x2)
	n0.add<Edge>(e0);										// N->E
	n1.add<Edge>(e0);										// N->E

	Edge e1 = m_meshQ->newEdge(n, n1);		// E->N (x2)
	n.add<Edge>(e1);										// N->E
	n1.add<Edge>(e1);										// N->E

	Edge e2 = m_meshQ->newEdge(n, n2);		// E->N (x2)
	n.add<Edge>(e2);										// N->E
	n2.add<Edge>(e2);										// N->E

	Edge e3 = m_meshQ->newEdge(n0, n2);	// E->N (x2)
	n0.add<Edge>(e3);										// N->E
	n2.add<Edge>(e3);										// N->E

	// Connectivités F->E
	f.add<Edge>(e0);		// F->E
	f.add<Edge>(e1);		// F->E
	f.add<Edge>(e2);		// F->E
	f.add<Edge>(e3);		// F->E

	// Connectivités E->F
	e0.add<Face>(f);		// E->F
	e1.add<Face>(f);		// E->F
	e2.add<Face>(f);		// E->F
	e3.add<Face>(f);		// E->F


	// Mise à jour du Front_IN
	Front_IN.setMultipleNode(n_id);
	Front_IN.setNextNode( n_id, neighbors_nodes[0], n1.id());
	Front_IN.setNextNode( n_id, neighbors_nodes[1], n2.id());
	Front_IN.setNonMultiplicable(n_id);
	Front_IN.setNonFusionable(n_id);
	Front_IN.setNonFusionable(neighbors_nodes[0]);
	Front_IN.setNonFusionable(neighbors_nodes[1]);


	// Mise à jour de l'indice de couche
	Variable<int>* couche_id = m_meshQ->getVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	couche_id->set(n1.id(), Front_IN.getFrontID()+1);
	couche_id->set(n2.id(), Front_IN.getFrontID()+1);

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroExtrusion_2D::Insertion_Double(Front &Front_IN, TCellID n_id,
                            Variable<double>* A_distance, double dist_cible, Variable<math::Vector3d>* A_vectors){

	std::cout << "INSERTION DOUBLE QUAD AU NOEUD : " << n_id << std::endl;

	std::vector<TCellID> neighbors_nodes = Front_IN.getNeighbors(n_id);
	Node n = m_meshQ->get<Node>(n_id);

	TCellID n_id_next = Front_IN.getIdealNode(n_id);
	Node n_next = m_meshQ->get<Node>(n_id_next);

	// Build the first quad
	// Build the first 2 nodes
	Node n_neighbor_1 = m_meshQ->get<Node>(neighbors_nodes[0]);
	TCellID n_id_next_neighbor_1 = Front_IN.getNextNode(neighbors_nodes[0],n_id);
	Node n_next_neighbor_1 = m_meshQ->get<Node>(n_id_next_neighbor_1);

	// Compute the positions of n1,Q1 and n2,Q1
	math::Vector3d v1 = m_params_aero.delta_cl  *(n_next_neighbor_1.point() - n_neighbor_1.point()).normalize() ;
	math::Point P1_Q1 = n.point() + v1 ;

	math::Vector3d v2 = (n_next.point() - n.point()).normalize() + (P1_Q1 - n.point()).normalize() ;
	v2 = ( n_next.point() - P1_Q1 ).norm() *v2.normalize();
	math::Point P2_Q1 = n.point() + v2 ;

	Node n1_Q1 = m_meshQ->newNode(P1_Q1);
	Node n2_Q1 = m_meshQ->newNode(P2_Q1);

	// Build the quad
	Face f = m_meshQ->newQuad(n, n_next, n2_Q1, n1_Q1) ; // F->N (x4)
	n.add<Face>(f);			// N->F
	n_next.add<Face>(f);		// N->F
	n2_Q1.add<Face>(f);		// N->F
	n1_Q1.add<Face>(f);		// N->F

	// Création des arêtes
	Edge e0 = m_meshQ->newEdge(n, n_next);	// E->N (x2)
	n.add<Edge>(e0);											// N->E
	n_next.add<Edge>(e0);									// N->E

	Edge e1 = m_meshQ->newEdge(n_next, n2_Q1);	// E->N (x2)
	n_next.add<Edge>(e1);										// N->E
	n2_Q1.add<Edge>(e1);											// N->E

	Edge e2 = m_meshQ->newEdge(n2_Q1, n1_Q1);	// E->N (x2)
	n2_Q1.add<Edge>(e2);											// N->E
	n1_Q1.add<Edge>(e2);											// N->E

	Edge e3 = m_meshQ->newEdge(n1_Q1, n);		// E->N (x2)
	n1_Q1.add<Edge>(e3);										// N->E
	n.add<Edge>(e3);											// N->E

	// Connectivités F->E
	f.add<Edge>(e0);		// F->E
	f.add<Edge>(e1);		// F->E
	f.add<Edge>(e2);		// F->E
	f.add<Edge>(e3);		// F->E

	// Connectivités E->F
	e0.add<Face>(f);		// E->F
	e1.add<Face>(f);		// E->F
	e2.add<Face>(f);		// E->F
	e3.add<Face>(f);		// E->F





	// Build the second quad
	// Build the first 2 nodes
	Node n_neighbor_2 = m_meshQ->get<Node>(neighbors_nodes[1]);
	TCellID n_id_next_neighbor_2 = Front_IN.getNextNode(neighbors_nodes[1],n_id);
	Node n_next_neighbor_2 = m_meshQ->get<Node>(n_id_next_neighbor_2);

	// Compute the positions of n1,Q1 and n2,Q1
	v1 = m_params_aero.delta_cl  *(n_next_neighbor_2.point() - n_neighbor_2.point()).normalize() ;
	math::Point P1_Q2 = n.point() + v1 ;

	v2 = (n_next.point() - n.point()).normalize() + (P1_Q2 - n.point()).normalize() ;
	v2 = ( n_next.point() - P1_Q2 ).norm() *v2.normalize();
	math::Point P2_Q2 = n.point() + v2 ;

	//math::Point P1_Q2 = n_next.point() + 0.5*(n_next_neighbor_2.point() - n_next.point() ) ;
	//math::Point P2_Q2 = n_next.point() + 0.25*(n_next_neighbor_2.point() - n_next.point() ) ;

	Node n1_Q2 = m_meshQ->newNode(P1_Q2);
	Node n2_Q2 = m_meshQ->newNode(P2_Q2);

	// Build the quad
	Face f2 = m_meshQ->newQuad(n, n_next, n2_Q2, n1_Q2) ; // F->N (x4)
	n.add<Face>(f2);			// N->F
	n_next.add<Face>(f2);		// N->F
	n2_Q2.add<Face>(f2);		// N->F
	n1_Q2.add<Face>(f2);		// N->F

	// Création des arêtes
	// The edge e0 is already known, it was created for the first quad.

	e1 = m_meshQ->newEdge(n_next, n2_Q2);	// E->N (x2)
	n_next.add<Edge>(e1);								// N->E
	n2_Q2.add<Edge>(e1);									// N->E

	e2 = m_meshQ->newEdge(n2_Q2, n1_Q2);	// E->N (x2)
	n2_Q2.add<Edge>(e2);									// N->E
	n1_Q2.add<Edge>(e2);									// N->E

	e3 = m_meshQ->newEdge(n1_Q2, n);		// E->N (x2)
	n1_Q2.add<Edge>(e3);									// N->E
	n.add<Edge>(e3);										// N->E

	// Connectivités F->E
	f2.add<Edge>(e0);		// F->E
	f2.add<Edge>(e1);		// F->E
	f2.add<Edge>(e2);		// F->E
	f2.add<Edge>(e3);		// F->E

	// Connectivités E->F
	e0.add<Face>(f2);		// E->F
	e1.add<Face>(f2);		// E->F
	e2.add<Face>(f2);		// E->F
	e3.add<Face>(f2);		// E->F




	// Update the Front_IN
	Front_IN.setMultipleNode(n_id);
	Front_IN.setNextNode( n_id, neighbors_nodes[0], n1_Q1.id());
	Front_IN.setNextNode( n_id, neighbors_nodes[1], n1_Q2.id());
	Front_IN.setNonMultiplicable(n_id);
	Front_IN.setNonFusionable(n_id);
	Front_IN.setNonFusionable(neighbors_nodes[0]);
	Front_IN.setNonFusionable(neighbors_nodes[1]);

	// Update the layer index
	Variable<int>* couche_id = m_meshQ->getVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	couche_id->set(n1_Q1.id(), Front_IN.getFrontID()+1);
	couche_id->set(n2_Q1.id(), Front_IN.getFrontID()+1);
	couche_id->set(n1_Q2.id(), Front_IN.getFrontID()+1);
	couche_id->set(n2_Q2.id(), Front_IN.getFrontID()+1);

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroExtrusion_2D::Fusion(Front &Front_IN, TCellID n_id){

	std::vector<TCellID> neighbors_nodes = Front_IN.getNeighbors(n_id);
	Node n = m_meshQ->get<Node>(n_id);

	// Création du quad
	TCellID n0_id = Front_IN.getIdealNode(n_id);
	Node n0 = m_meshQ->get<Node>(n0_id);

	Node n1 = m_meshQ->get<Node>(neighbors_nodes[0]);
	Node n2 = m_meshQ->get<Node>(neighbors_nodes[1]);

	// On créé la face
	Face f = m_meshQ->newQuad(n, n1, n0, n2) ; // F->N (x4)
	n.add<Face>(f);		// N->F
	n0.add<Face>(f);		// N->F
	n1.add<Face>(f);		// N->F
	n2.add<Face>(f);		// N->F

	// Récupération des 2 arêtes existantes
	TCellID e1_id = math::Utils::CommonEdge(m_meshQ, n_id, n1.id()) ;
	Edge e1 = m_meshQ->get<Edge>(e1_id);

	TCellID e2_id = math::Utils::CommonEdge(m_meshQ, n_id, n2.id()) ;
	Edge e2 = m_meshQ->get<Edge>(e2_id);

	// Création des deux arêtes manquantes
	Edge e0 = m_meshQ->newEdge(n0, n1);	// E->N (x2)
	n0.add<Edge>(e0);										// N->E
	n1.add<Edge>(e0);										// N->E

	Edge e3 = m_meshQ->newEdge(n0, n2);	// E->N (x2)
	n0.add<Edge>(e3);										// N->E
	n2.add<Edge>(e3);										// N->E

	// Connectivités F->E
	f.add<Edge>(e0);		// F->E
	f.add<Edge>(e1);		// F->E
	f.add<Edge>(e2);		// F->E
	f.add<Edge>(e3);		// F->E

	// Connectivités E->F
	e0.add<Face>(f);		// E->F
	e1.add<Face>(f);		// E->F
	e2.add<Face>(f);		// E->F
	e3.add<Face>(f);		// E->F


	// Mise à jour du Front_IN
	Front_IN.setContractedNode(n_id);
	std::vector<TCellID> n1_neighbors = Front_IN.getNeighbors(neighbors_nodes[0]);
	Front_IN.setNextNode( neighbors_nodes[0], n1_neighbors[0], n0.id());
	Front_IN.setNextNode( neighbors_nodes[0], n1_neighbors[1], n0.id());
	std::vector<TCellID> n2_neighbors = Front_IN.getNeighbors(neighbors_nodes[1]);
	Front_IN.setNextNode( neighbors_nodes[1], n2_neighbors[0], n0.id());
	Front_IN.setNextNode( neighbors_nodes[1], n2_neighbors[1], n0.id());
	Front_IN.setNonMultiplicable(n_id);
	Front_IN.setNonFusionable(n_id);
	Front_IN.setNonMultiplicable(neighbors_nodes[0]);
	Front_IN.setNonMultiplicable(neighbors_nodes[1]);
	Front_IN.setNonFusionable(neighbors_nodes[0]);
	Front_IN.setNonFusionable(neighbors_nodes[1]);
	Front_IN.setNonFusionable(n1_neighbors[0]);
	Front_IN.setNonFusionable(n1_neighbors[1]);
	Front_IN.setNonFusionable(n2_neighbors[0]);
	Front_IN.setNonFusionable(n2_neighbors[1]);


	/*
	// Mise à jour du Front_OUT
	Front_OUT.addNodeId(n0.id());
	 */

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroExtrusion_2D::CreateNormalQuad(TCellID e_id, Front &Front_IN){

	Edge e = m_meshQ->get<Edge>(e_id);
	std::vector<Node> nodes = e.get<Node>();

	TCellID n0_id = Front_IN.getNextNode(nodes[0].id(), nodes[1].id());
	TCellID n1_id = Front_IN.getNextNode(nodes[1].id(), nodes[0].id());
	Node n0 = m_meshQ->get<Node>(n0_id);
	Node n1 = m_meshQ->get<Node>(n1_id);


	// On créé la face associée à l'arête
	Face f = m_meshQ->newQuad(nodes[0], nodes[1], n1, n0) ; // F->N (x4)
	nodes[0].add<Face>(f);	// N->F
	nodes[1].add<Face>(f);	// N->F
	n0.add<Face>(f);			// N->F
	n1.add<Face>(f);			// N->F

	// Création de l'arête sur le nouveau front
	Edge e_opp = m_meshQ->newEdge(n0, n1);	// E->N (x2)
	n0.add<Edge>(e_opp);										// N->E
	n1.add<Edge>(e_opp);										// N->E

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

}
/*------------------------------------------------------------------------*/