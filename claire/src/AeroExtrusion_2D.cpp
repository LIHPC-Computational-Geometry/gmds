//
// Created by rochec on 14/04/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/claire/Utils.h>
#include <gmds/claire/AeroException.h>
#include <gmds/claire/AeroExtrusion_2D.h>
#include <gmds/claire/AdvectedPointRK4_2D.h>
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

	Compute1stLayer(m_meshT->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"), 0.25);

	return AeroExtrusion_2D::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
Front
AeroExtrusion_2D::Compute1stLayer(Variable<double>* A_distance, double dist_cible){

	Front Front_Paroi;
	Front_Paroi.initializeFromLayerId(m_meshQ, 0);
	Front First_Front;

	std::map<TCellID, TCellID> map_new_nodes = ComputeIdealPositions(Front_Paroi, dist_cible,
	                                                                 m_meshT->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
	                                                                    m_meshT->getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));

	std::vector<TCellID> front_nodes = Front_Paroi.getNodes();
	std::vector<TCellID> front_edges = Front_Paroi.getEdges();

	Variable<int>* couche_id = m_meshQ->getVariable<int, GMDS_NODE>("GMDS_Couche_Id");

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