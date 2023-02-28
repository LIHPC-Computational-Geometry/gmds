//
// Created by rochec on 23/03/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AeroBoundaries_2D.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AeroBoundaries_2D::AeroBoundaries_2D(Mesh *AMesh) :
  AbstractAeroBoundaries(AMesh)
{

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroBoundaries_2D::MarkBoundariesNodes(){

	//Get the boundary node ids
	BoundaryOperator2D bnd_op(m_mesh);
	//std::vector<TCellID> bnd_node_ids;
	bnd_op.getBoundaryNodes(m_bnd_nodes_ids);

	// Initialise une marque sur les noeuds du bord
	for (auto n_id:m_bnd_nodes_ids){
		Node n = m_mesh->get<Node>(n_id);
		m_mesh->mark(n, m_markBoundaryNodes);
	}

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroBoundaries_2D::WhichColorIsAmont(){

	// Calcul des boites englobantes
	std::vector<double> x_min(m_nbrBords,0.0);
	std::vector<double> y_min(m_nbrBords,0.0);
	std::vector<double> x_max(m_nbrBords,0.0);
	std::vector<double> y_max(m_nbrBords,0.0);

	for (auto n_id:m_bnd_nodes_ids ) {
		Node n = m_mesh->get<Node>(n_id);
		math::Point p = n.point();
		// A quel bord appartient le noeud n_id
		int couleur = m_var_color_bords->value(n_id);
		if(p.X() < x_min[couleur-1]){
			x_min[couleur-1] = p.X();
		}
		if(p.X() > x_max[couleur-1]){
			x_max[couleur-1] = p.X();
		}
		if(p.Y() < y_min[couleur-1]){
			y_min[couleur-1] = p.Y();
		}
		if(p.Y() > y_max[couleur-1]){
			y_max[couleur-1] = p.Y();
		}
	}

	// On recherche la plus grosse boite englobante
	m_color_Amont = 1;
	for(int i=2;i<=m_nbrBords;i++){
		// On compte le nombre de boites que la couleur i englobe
		int nbr_boite_englobe = 0;
		for (int j=1;j<=m_nbrBords;j++){
			if (i != j){
				if ( (x_min[i-1] < x_min[j-1]) &&
				    (y_min[i-1] < y_min[j-1]) &&
				    (x_max[j-1] < x_max[i-1]) &&
				    (y_max[j-1] < y_max[i-1]) ) {
					nbr_boite_englobe += 1;
				}
			}
		}
		if(nbr_boite_englobe == m_nbrBords-1){
			m_color_Amont = i;
		}
	}

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
std::vector<TCellID> AeroBoundaries_2D::BndNodesOrdered(int color){

	std::vector<TCellID> bnd_nodes_id_ordered;

	//Get the boundary node ids
	BoundaryOperator2D bnd_op(m_mesh);
	std::vector<TCellID> bnd_node_ids;
	bnd_op.getBoundaryNodes(bnd_node_ids);

	int markCurveEdges = m_mesh->newMark<Edge>();
	int markCurveNodes = m_mesh->newMark<Node>();
	int markPointNodes = m_mesh->newMark<Node>();
	int markAloneNodes = m_mesh->newMark<Node>();

	bnd_op.markCellOnGeometry(markCurveEdges, markCurveNodes,
	                          markPointNodes, markAloneNodes);

	m_mesh->unmarkAll<Node>(markAloneNodes);
	m_mesh->freeMark<Node>(markAloneNodes);
	m_mesh->unmarkAll<Edge>(markCurveEdges);
	m_mesh->freeMark<Edge>(markCurveEdges);
	m_mesh->unmarkAll<Node>(markCurveNodes);
	m_mesh->freeMark<Node>(markCurveNodes);

	TCellID n0_id = NullID ; // Choix du noeud de départ

	// STRATEGIE 1 : On récuère un noeud de ce bord qui est sur un sommet
	for(auto n_it = m_mesh->nodes_begin(); n_it!= m_mesh->nodes_end() && (n0_id == NullID);++n_it){
		TCellID n_id = *n_it;
		//if ( m_mesh->isMarked<Node>(n_id, m_markNodesParoi) &&
		if ( m_var_color_bords->value(n_id) == color &&
		    m_mesh->isMarked<Node>(n_id, markPointNodes) ){
			n0_id = n_id;
		}
	}

	m_mesh->unmarkAll<Node>(markPointNodes);
	m_mesh->freeMark<Node>(markPointNodes);

	// STRATEGIE 2 : On prend le point d'arrêt, positionné au x_min.
	// Attention, plusieurs noeuds peuvent être au x_min. Cette méthode
	// en retourne un au hasard.
	n0_id = PointArret(color);

	// STRATEGIE 3 : On sélectionne un noeud au hasard
	/*
	for (auto n_it = m_mesh->nodes_begin(); n_it != m_mesh->nodes_end() && (n0_id == NullID); ++n_it) {
	   TCellID n_id = *n_it;
	   if (m_Bnd->isParoi(n_id) &&
	      m_Bnd->getNodeColor(n_id) == color) {
	      n0_id = n_id;
	   }
	}
	 */

	// Initialisation de la première valeur du vecteur
	bnd_nodes_id_ordered.push_back(n0_id);

	// Initialisation des valeurs pour le parcours du bord
	TCellID n1_id;
	TCellID n2_id = n0_id;
	TCellID n3_id = n0_id;

	Node n0 = m_mesh->get<Node>(n0_id);
	std::vector<Edge> adj_edges = n0.get<Edge>();
	for (auto const &e:adj_edges){
		Node ne = e.getOppositeNode(n0);
		if ( m_mesh->isMarked<Node>(ne.id(), m_markNodesParoi) ){
			n3_id = ne.id();
		}
	}

	while (n3_id != n0_id){

		bnd_nodes_id_ordered.push_back(n3_id);		// On ajout l'id du noeud n3 au vecteur

		// On cherche le prochain noeud
		n1_id = n2_id;
		n2_id = n3_id;
		Node n2 = m_mesh->get<Node>(n2_id);
		adj_edges = n2.get<Edge>();
		for (auto const &e:adj_edges){
			Node ne = e.getOppositeNode(n2);
			if ( m_mesh->isMarked<Node>(ne.id(), m_markNodesParoi) &&
			    ne.id() != n1_id ){
				n3_id = ne.id();
			}
		}

	}

	return bnd_nodes_id_ordered;

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
double AeroBoundaries_2D::ComputeBoundaryLength(int color){

	double length(0.0);
	std::vector<TCellID> bnd_nodes_id_ordered = BndNodesOrdered(color);

	for(int i=0;i<bnd_nodes_id_ordered.size()-1;i++){
		Node n0 = m_mesh->get<Node>(bnd_nodes_id_ordered[i]);
		Node n1 = m_mesh->get<Node>(bnd_nodes_id_ordered[i+1]);
		math::Point p0 = n0.point();
		math::Point p1 = n1.point();
		math::Vector3d v = p1-p0;
		length += v.norm();
	}

	// Ajout de la dernière longueur pour "fermer la boucle"
	Node n0 = m_mesh->get<Node>(bnd_nodes_id_ordered[bnd_nodes_id_ordered.size()-1]);
	Node n1 = m_mesh->get<Node>(bnd_nodes_id_ordered[0]);
	math::Point p0 = n0.point();
	math::Point p1 = n1.point();
	math::Vector3d v = p1-p0;
	length += v.norm();

	return length;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
std::vector<TCellID> AeroBoundaries_2D::BndEdges(int color)
{
	std::vector<TCellID> edges_id;

	for (auto n_id:m_mesh->nodes())
	{
		m_map_color_bords[n_id] = m_var_color_bords->value(n_id);
	}

	for (auto e_id:m_mesh->edges())
	{
		Edge e = m_mesh->get<Edge>(e_id);
		std::vector<Node> nodes = e.get<Node>();
		if (m_map_color_bords[nodes[0].id()] == color
		       && m_map_color_bords[nodes[1].id()] == color)
		{
			edges_id.push_back(e_id);
		}
	}

	return edges_id;
}
/*------------------------------------------------------------------------*/