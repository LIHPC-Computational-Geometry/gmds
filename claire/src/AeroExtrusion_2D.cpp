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
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AeroExtrusion_2D::AeroExtrusion_2D(Mesh *AMeshT, Mesh *AMeshQ) {
	m_meshT = AMeshT;
	m_meshQ = AMeshQ;
}


/*------------------------------------------------------------------------*/
AeroExtrusion_2D::STATUS
AeroExtrusion_2D::execute()
{
	// Exemple exception
	//if(m_mesh==NULL)
	//	throw AeroException("ERROR: Invalid mesh pointer");

	Front Current_Front = Compute1stLayer(m_meshT->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"), 0.5,
	                m_meshT->getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	Current_Front = ComputeLayer(Current_Front, m_meshT->getVariable<double,GMDS_NODE>("GMDS_Distance"), 1,
	                             m_meshT->getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));

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

	Front Front_Paroi;
	Front_Paroi.initializeFromLayerId(m_meshQ, 0);
	Front First_Front;

	std::map<TCellID, TCellID> map_new_nodes = ComputeIdealPositions(Front_Paroi, dist_cible, A_distance, A_vectors);

	std::vector<TCellID> front_nodes = Front_Paroi.getNodes();
	std::vector<TCellID> front_edges = Front_Paroi.getEdges();

	Variable<int>* couche_id = m_meshQ->getVariable<int, GMDS_NODE>("GMDS_Couche_Id");

	// Remarque : on pourait initialiser les arêtes e0 et e1 qui relient les deux couches dans cette boucle.
	for (auto n_id:front_nodes){
		// Ajout de la nouvelle arête du front dans l'objet concerné
		First_Front.addNodeId(map_new_nodes[n_id]);
		couche_id->set(map_new_nodes[n_id], 1);
	}

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

	return First_Front;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
Front
AeroExtrusion_2D::ComputeLayer(Front Front_IN, Variable<double>* A_distance, double dist_cible, Variable<math::Vector3d>* A_vectors){
	Front Front_OUT;

	std::map<TCellID, TCellID> map_new_nodes = ComputeIdealPositions(Front_IN, dist_cible, A_distance, A_vectors);

	Front_IN.initializeNodeType(m_meshQ, map_new_nodes);

	std::vector<TCellID> front_nodes = Front_IN.getNodes();
	std::vector<TCellID> front_edges = Front_IN.getEdges();

	// Boucle d'insertion d'éléments
	bool do_smooth(true);

	while (do_smooth){
		do_smooth = false;
		// A FAIRE : lissage de la couche

		TCellID node_id;
		int type_node(0);

		getSingularNode(Front_IN, node_id, type_node);

		if (type_node == 1){
			// Insertion
			Insertion(Front_IN, node_id, Front_OUT,
			          A_distance, dist_cible, A_vectors);
			do_smooth = true;
		}
		else if (type_node == 2){
			// Fusion
			//do_smooth = true;
		}
	}

	return Front_OUT;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroExtrusion_2D::getSingularNode(Front Front_IN, TCellID &node_id, int &type){

	std::vector<TCellID> front_nodes = Front_IN.getNodes();
	std::vector<TCellID> front_edges = Front_IN.getEdges();

	node_id = NullID;
	type = 0;

	for (auto n_id:front_nodes){
		// si le noeud est normal
		if (Front_IN.getNodeType(n_id) == 0){
			std::vector<TCellID> neighbors_nodes = Front_IN.getNeighbors(n_id);
			// test premier quad
			double r1 = math::AeroMeshQuality::oppositeedgeslenghtratio(m_meshT, n_id, neighbors_nodes[0],
			                                                Front_IN.getNextNode(neighbors_nodes[0],n_id),
			                                                Front_IN.getNextNode(n_id,neighbors_nodes[0]));
			// test second quad
			double r2 = math::AeroMeshQuality::oppositeedgeslenghtratio(m_meshT, n_id, neighbors_nodes[1],
			                                                            Front_IN.getNextNode(neighbors_nodes[1],n_id),
			                                                            Front_IN.getNextNode(n_id,neighbors_nodes[1]));
			std::cout << "r1 et r2 : " << r1 << " " << r2 << std::endl;
			if (r1 < 0.95 || r2 < 0.95){
				node_id = n_id;
				type = 1;
			}
		}
	}


}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroExtrusion_2D::Insertion(Front &Front_IN, TCellID n_id, Front &Front_OUT,
                            Variable<double>* A_distance, double dist_cible, Variable<math::Vector3d>* A_vectors){

	std::vector<TCellID> neighbors_nodes = Front_IN.getNeighbors(n_id);

	// Contruction du premier noeud n1
	Node n = m_meshQ->get<Node>(n_id);
	Node n_neighbor = m_meshQ->get<Node>(neighbors_nodes[0]);
	std::cout << "Premier point : " << n_neighbor.point() << " " << n_neighbor.id() << std::endl;
	std::cout << "Second point : " << n.point() << " " << n.id() << std::endl;
	math::Point P1_construction = n.point() + 0.3*(n_neighbor.point()-n.point()) ;
	std::cout << "P1 : " << P1_construction << std::endl;
	AdvectedPointRK4_2D advpoint_n1(m_meshT, P1_construction, dist_cible, A_distance, A_vectors);
	advpoint_n1.execute();
	Node n1 = m_meshQ->newNode(advpoint_n1.getPend());
	// Contruction du premier noeud n2
	n_neighbor = m_meshQ->get<Node>(neighbors_nodes[1]);
	P1_construction = n.point() + 0.3*(n_neighbor.point()-n.point()) ;
	AdvectedPointRK4_2D advpoint_n2(m_meshT, P1_construction, dist_cible, A_distance, A_vectors);
	advpoint_n2.execute();
	Node n2 = m_meshQ->newNode(advpoint_n2.getPend());

	// Mise à jour du Front_IN
	Front_IN.setMultipleNode(n_id);
	Front_IN.setNextNode(neighbors_nodes[0], n_id, n1.id());
	Front_IN.setNextNode(neighbors_nodes[1], n_id, n2.id());

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


	// Mise à jour du Front_OUT
	Front_OUT.addEdgeId(e0.id());
	Front_OUT.addEdgeId(e3.id());
	Front_OUT.addNodeId(n0.id());
	Front_OUT.addNodeId(n1.id());
	Front_OUT.addNodeId(n2.id());

}
/*------------------------------------------------------------------------*/