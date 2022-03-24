//
// Created by rochec on 23/03/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/Utils.h>
#include <gmds/claire/AeroBoundaries_2D.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <iostream>
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