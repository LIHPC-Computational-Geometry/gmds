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

AeroExtrusion_2D::AeroExtrusion_2D(Mesh *AMeshT, Mesh *AMeshQ, ParamsAero Aparams_aero) {
	m_meshT = AMeshT;
	m_meshQ = AMeshQ;
	m_params_aero = Aparams_aero;
}


/*------------------------------------------------------------------------*/
AeroExtrusion_2D::STATUS
AeroExtrusion_2D::execute()
{
	// Exemple exception
	//if(m_mesh==NULL)
	//	throw AeroException("ERROR: Invalid mesh pointer");

	Front Current_Front = Compute1stLayer(m_meshT->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"), m_params_aero.delta_cl,
	                m_meshT->getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
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
		Current_Front = ComputeLayer(Current_Front, m_meshT->getVariable<double,GMDS_NODE>("GMDS_Distance_3"), i*pas_couche,
	                             m_meshT->getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
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
		AdvectedPointRK4_2D advpoint(m_meshT, M, dist_cible, A_distance, A_vectors);
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

	Front First_Front;
	First_Front.setFrontID(Front_Paroi.getFrontID()+1);

	std::map<TCellID, TCellID> map_new_nodes = ComputeIdealPositions(Front_Paroi, dist_cible, A_distance, A_vectors);

	std::vector<TCellID> front_nodes = Front_Paroi.getNodes();
	std::vector<TCellID> front_edges = Front_Paroi.getEdges();

	Variable<int>* couche_id = m_meshQ->getVariable<int, GMDS_NODE>("GMDS_Couche_Id");

	// Remarque : on pourait initialiser les ar??tes e0 et e1 qui relient les deux couches dans cette boucle.
	for (auto n_id:front_nodes){
		// Ajout de la nouvelle ar??te du front dans l'objet concern??
		First_Front.addNodeId(map_new_nodes[n_id]);
		couche_id->set(map_new_nodes[n_id], First_Front.getFrontID());
	}

	for (auto e_id:front_edges){

		Edge e = m_meshQ->get<Edge>(e_id);
		std::vector<Node> nodes = e.get<Node>();

		Node n0 = m_meshQ->get<Node>(map_new_nodes[nodes[0].id()]);
		Node n1 = m_meshQ->get<Node>(map_new_nodes[nodes[1].id()]);


		// On cr???? la face associ??e ?? l'ar??te
		Face f = m_meshQ->newQuad(nodes[0], nodes[1], n1, n0) ; // F->N (x4)
		nodes[0].add<Face>(f);	// N->F
		nodes[1].add<Face>(f);	// N->F
		n0.add<Face>(f);			// N->F
		n1.add<Face>(f);			// N->F

		// Cr??ation de l'ar??te sur le nouveau front
		Edge e_opp = m_meshQ->newEdge(n0, n1);	// E->N (x2)
		n0.add<Edge>(e_opp);										// N->E
		n1.add<Edge>(e_opp);										// N->E

		TCellID e0_id = math::Utils::CommonEdge(m_meshQ, nodes[0].id(), n0.id());
		TCellID e1_id = math::Utils::CommonEdge(m_meshQ, nodes[1].id(), n1.id());
		Edge e0, e1;
		if ( e0_id == NullID ){
			// Si l'ar??te n'existe pas, on la cr???? et on initialise les connectivit??s N->E
			e0 = m_meshQ->newEdge(nodes[0].id(), n0.id());		// E->N (x2)
			nodes[0].add<Edge>(e0);												// N->E
			n0.add<Edge>(e0);														// N->E
		}
		else{
			e0 = m_meshQ->get<Edge>(e0_id);
		}

		// Idem pour l'ar??te e1
		if ( e1_id == NullID ){
			e1 = m_meshQ->newEdge(nodes[1].id(), n1.id());		// E->N (x2)
			nodes[1].add<Edge>(e1);												// N->E
			n1.add<Edge>(e1);														// N->E
		}
		else{
			e1 = m_meshQ->get<Edge>(e1_id);
		}

		// Connectivit??s F->E
		f.add<Edge>(e);		// F->E
		f.add<Edge>(e0);		// F->E
		f.add<Edge>(e1);		// F->E
		f.add<Edge>(e_opp);	// F->E

		// Connectivit??s E->F
		e.add<Face>(f);		// E->F
		e0.add<Face>(f);		// E->F
		e1.add<Face>(f);		// E->F
		e_opp.add<Face>(f);	// E->F


		// Ajout de la nouvelle ar??te du front dans l'objet concern??
		First_Front.addEdgeId(e_opp.id());

	}

	return First_Front;

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
Front
AeroExtrusion_2D::ComputeLayer(Front Front_IN, Variable<double>* A_distance, double dist_cible, Variable<math::Vector3d>* A_vectors){

	std::cout << "--------- couche " << Front_IN.getFrontID()+1 << std::endl;

	std::map<TCellID, TCellID> map_new_nodes = ComputeIdealPositions(Front_IN, dist_cible, A_distance, A_vectors);

	// Mise ?? jour de l'indice de couche
	Variable<int>* couche_id = m_meshQ->getVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	for (auto n_id:Front_IN.getNodes()){
		couche_id->set(map_new_nodes[n_id], Front_IN.getFrontID()+1);
	}

	Front_IN.initializeNodeType(m_meshQ, map_new_nodes);

	std::vector<TCellID> front_nodes = Front_IN.getNodes();
	std::vector<TCellID> front_edges = Front_IN.getEdges();

	Front_IN.setMultorFusFromLimits(m_meshQ, m_params_aero.x_lim, m_params_aero.y_lim, m_params_aero.z_lim);

	// Boucle d'insertion d'??l??ments
	bool do_smooth(true);

	while (do_smooth){
		do_smooth = false;
		// A FAIRE : lissage de la couche

		TCellID node_id;
		int type_node(0);
		getSingularNode(Front_IN, node_id, type_node);

		if (type_node == 1){
			// Insertion
			std::cout << "INSERTION QUAD AU NOEUD : " << node_id << std::endl;
			Insertion(Front_IN, node_id, A_distance, dist_cible, A_vectors);
			do_smooth = true;
		}
		else if (type_node == 2){
			// Fusion
			std::cout << "FUSION QUAD AU NOEUD : " << node_id << std::endl;
			Fusion(Front_IN, node_id);
			do_smooth = true;
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

	// Supression des noeuds non utilis??s
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

	for (auto n_id:front_nodes){

		// Les tests pour la FUSION
		std::vector<TCellID> neighbors_nodes = Front_IN.getNeighbors(n_id);
		if (Front_IN.isFusionable(n_id)){

			// test ar??te ??cras??e premier quad
			double min_lenght_1 = math::AeroMeshQuality::minlenghtedge(m_meshQ, n_id, neighbors_nodes[0],
			                                  Front_IN.getNextNode(neighbors_nodes[0],n_id),
			                                  Front_IN.getNextNode(n_id,neighbors_nodes[0]));

			// test ar??te ??cras??e second quad
			double min_lenght_2 = math::AeroMeshQuality::minlenghtedge(m_meshQ, n_id, neighbors_nodes[1],
			                                                           Front_IN.getNextNode(neighbors_nodes[1],n_id),
			                                                           Front_IN.getNextNode(n_id,neighbors_nodes[1]));

			if (singu_not_found && (min_lenght_1 < 0.001 || min_lenght_2 < 0.001) ){
				node_id = n_id;
				type = 2;
				singu_not_found = false;
			}

			// test angles
			double angle = math::AeroMeshQuality::angleouverture(m_meshQ, n_id,
			                                                     Front_IN.getIdealNode(n_id),
			                                                     neighbors_nodes[0], neighbors_nodes[1]) ;
			//std::cout << "angle : " << angle << std::endl;
			if (singu_not_found
			    && ( angle < 3.0*M_PI/4.0 )){
				node_id = n_id;
				type = 2;
				singu_not_found = false;
			}

		}

		// Les tests pour l'INSERSION
		if(singu_not_found && Front_IN.isMultiplicable(n_id)) {

			// test angles
			{
				double angle = math::AeroMeshQuality::angleouverture(m_meshQ, n_id,
				                                                     Front_IN.getIdealNode(n_id), neighbors_nodes[0], neighbors_nodes[1]);
				if (singu_not_found && (angle > 7.0*M_PI/6.0)) {
					node_id = n_id;
					type = 1;
					singu_not_found = false;
				}
			}

			// Internal Angle Deviation
			/*
			{
			   double iad1 = math::AeroMeshQuality::InternalAngleDeviationQUAD(m_meshQ, n_id, neighbors_nodes[0], Front_IN.getNextNode(neighbors_nodes[0], n_id),
			                                                                   Front_IN.getNextNode(n_id, neighbors_nodes[0]));
			   double iad2 = math::AeroMeshQuality::InternalAngleDeviationQUAD(m_meshQ, n_id, neighbors_nodes[1], Front_IN.getNextNode(neighbors_nodes[1], n_id),
			                                                                   Front_IN.getNextNode(n_id, neighbors_nodes[1]));

			   // std::cout << "angle : " << angle << std::endl;
			   if (singu_not_found && (iad1 > 40.0 || iad2 > 40.0)) {
			      node_id = n_id;
			      type = 1;
			      singu_not_found = false;
			   }
			}
			*/

			// Aspect Ratio
			/*
			{
				double ar1 = math::AeroMeshQuality::AspectRatioQUAD(m_meshQ, n_id, neighbors_nodes[0], Front_IN.getNextNode(neighbors_nodes[0], n_id),
																					 Front_IN.getNextNode(n_id, neighbors_nodes[0]));
				double ar2 = math::AeroMeshQuality::AspectRatioQUAD(m_meshQ, n_id, neighbors_nodes[1], Front_IN.getNextNode(neighbors_nodes[1], n_id),
																					 Front_IN.getNextNode(n_id, neighbors_nodes[1]));
				if (singu_not_found && (ar1 < 0.88 && ar2 < 0.88)) {
					node_id = n_id;
					type = 1;
					singu_not_found = false;
				}
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
	TCellID n_next_ideal_id = Front_IN.getIdealNode(n_id);
	TCellID n_next_id = Front_IN.getNextNode(neighbors_nodes[0],n_id);
	Node n_next_ideal = m_meshQ->get<Node>(n_next_ideal_id) ;
	Node n_next = m_meshQ->get<Node>(n_next_id) ;
	math::Vector3d P1_construction = 0.33*(n_next.point()-n_next_ideal.point()) ;
	Node n1 = m_meshQ->newNode(n_next_ideal.point() + P1_construction);

	// Contruction du premier noeud n2
	n_neighbor = m_meshQ->get<Node>(neighbors_nodes[1]);
	/*
	P1_construction = n.point() + 0.33*(n_neighbor.point()-n.point()) ;
	AdvectedPointRK4_2D advpoint_n2(m_meshT, P1_construction, dist_cible, A_distance, A_vectors);
	advpoint_n2.execute();
	Node n2 = m_meshQ->newNode(advpoint_n2.getPend());
	 */
	n_next_ideal_id = Front_IN.getIdealNode(n_id);
	n_next_id = Front_IN.getNextNode(neighbors_nodes[1],n_id);
	n_next_ideal = m_meshQ->get<Node>(n_next_ideal_id) ;
	n_next = m_meshQ->get<Node>(n_next_id) ;
	P1_construction = 0.33*(n_next.point()-n_next_ideal.point()) ;
	Node n2 = m_meshQ->newNode(n_next_ideal.point() + P1_construction);

	// Cr??ation du quad
	TCellID n0_id = Front_IN.getIdealNode(n_id);
	Node n0 = m_meshQ->get<Node>(n0_id);

	// On cr???? la face associ??e ?? l'ar??te
	Face f = m_meshQ->newQuad(n, n1, n0, n2) ; // F->N (x4)
	n.add<Face>(f);		// N->F
	n0.add<Face>(f);		// N->F
	n1.add<Face>(f);		// N->F
	n2.add<Face>(f);		// N->F

	// Cr??ation des ar??tes
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

	// Connectivit??s F->E
	f.add<Edge>(e0);		// F->E
	f.add<Edge>(e1);		// F->E
	f.add<Edge>(e2);		// F->E
	f.add<Edge>(e3);		// F->E

	// Connectivit??s E->F
	e0.add<Face>(f);		// E->F
	e1.add<Face>(f);		// E->F
	e2.add<Face>(f);		// E->F
	e3.add<Face>(f);		// E->F


	// Mise ?? jour du Front_IN
	Front_IN.setMultipleNode(n_id);
	Front_IN.setNextNode( n_id, neighbors_nodes[0], n1.id());
	Front_IN.setNextNode( n_id, neighbors_nodes[1], n2.id());
	Front_IN.setNonMultiplicable(n_id);
	Front_IN.setNonFusionable(n_id);
	Front_IN.setNonFusionable(neighbors_nodes[0]);
	Front_IN.setNonFusionable(neighbors_nodes[1]);


	/*
	// Mise ?? jour du Front_OUT
	Front_OUT.addEdgeId(e0.id());
	Front_OUT.addEdgeId(e3.id());
	Front_OUT.addNodeId(n0.id());
	Front_OUT.addNodeId(n1.id());
	Front_OUT.addNodeId(n2.id());
	 */


	// Mise ?? jour de l'indice de couche
	Variable<int>* couche_id = m_meshQ->getVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	couche_id->set(n1.id(), Front_IN.getFrontID()+1);
	couche_id->set(n2.id(), Front_IN.getFrontID()+1);

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroExtrusion_2D::Fusion(Front &Front_IN, TCellID n_id){

	std::vector<TCellID> neighbors_nodes = Front_IN.getNeighbors(n_id);
	Node n = m_meshQ->get<Node>(n_id);

	// Cr??ation du quad
	TCellID n0_id = Front_IN.getIdealNode(n_id);
	Node n0 = m_meshQ->get<Node>(n0_id);

	Node n1 = m_meshQ->get<Node>(neighbors_nodes[0]);
	Node n2 = m_meshQ->get<Node>(neighbors_nodes[1]);

	// On cr???? la face
	Face f = m_meshQ->newQuad(n, n1, n0, n2) ; // F->N (x4)
	n.add<Face>(f);		// N->F
	n0.add<Face>(f);		// N->F
	n1.add<Face>(f);		// N->F
	n2.add<Face>(f);		// N->F

	// R??cup??ration des 2 ar??tes existantes
	TCellID e1_id = math::Utils::CommonEdge(m_meshQ, n_id, n1.id()) ;
	Edge e1 = m_meshQ->get<Edge>(e1_id);

	TCellID e2_id = math::Utils::CommonEdge(m_meshQ, n_id, n2.id()) ;
	Edge e2 = m_meshQ->get<Edge>(e2_id);

	// Cr??ation des deux ar??tes manquantes
	Edge e0 = m_meshQ->newEdge(n0, n1);	// E->N (x2)
	n0.add<Edge>(e0);										// N->E
	n1.add<Edge>(e0);										// N->E

	Edge e3 = m_meshQ->newEdge(n0, n2);	// E->N (x2)
	n0.add<Edge>(e3);										// N->E
	n2.add<Edge>(e3);										// N->E

	// Connectivit??s F->E
	f.add<Edge>(e0);		// F->E
	f.add<Edge>(e1);		// F->E
	f.add<Edge>(e2);		// F->E
	f.add<Edge>(e3);		// F->E

	// Connectivit??s E->F
	e0.add<Face>(f);		// E->F
	e1.add<Face>(f);		// E->F
	e2.add<Face>(f);		// E->F
	e3.add<Face>(f);		// E->F


	// Mise ?? jour du Front_IN
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
	// Mise ?? jour du Front_OUT
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


	// On cr???? la face associ??e ?? l'ar??te
	Face f = m_meshQ->newQuad(nodes[0], nodes[1], n1, n0) ; // F->N (x4)
	nodes[0].add<Face>(f);	// N->F
	nodes[1].add<Face>(f);	// N->F
	n0.add<Face>(f);			// N->F
	n1.add<Face>(f);			// N->F

	// Cr??ation de l'ar??te sur le nouveau front
	Edge e_opp = m_meshQ->newEdge(n0, n1);	// E->N (x2)
	n0.add<Edge>(e_opp);										// N->E
	n1.add<Edge>(e_opp);										// N->E

	TCellID e0_id = math::Utils::CommonEdge(m_meshQ, nodes[0].id(), n0.id());
	TCellID e1_id = math::Utils::CommonEdge(m_meshQ, nodes[1].id(), n1.id());
	Edge e0, e1;
	if ( e0_id == NullID ){
		// Si l'ar??te n'existe pas, on la cr???? et on initialise les connectivit??s N->E
		e0 = m_meshQ->newEdge(nodes[0].id(), n0.id());		// E->N (x2)
		nodes[0].add<Edge>(e0);												// N->E
		n0.add<Edge>(e0);														// N->E
	}
	else{
		e0 = m_meshQ->get<Edge>(e0_id);
	}

	// Idem pour l'ar??te e1
	if ( e1_id == NullID ){
		e1 = m_meshQ->newEdge(nodes[1].id(), n1.id());		// E->N (x2)
		nodes[1].add<Edge>(e1);												// N->E
		n1.add<Edge>(e1);														// N->E
	}
	else{
		e1 = m_meshQ->get<Edge>(e1_id);
	}

	// Connectivit??s F->E
	f.add<Edge>(e);		// F->E
	f.add<Edge>(e0);		// F->E
	f.add<Edge>(e1);		// F->E
	f.add<Edge>(e_opp);	// F->E

	// Connectivit??s E->F
	e.add<Face>(f);		// E->F
	e0.add<Face>(f);		// E->F
	e1.add<Face>(f);		// E->F
	e_opp.add<Face>(f);	// E->F

}
/*------------------------------------------------------------------------*/