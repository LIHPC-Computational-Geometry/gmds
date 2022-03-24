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
AbstractAeroBoundaries::STATUS AeroBoundaries_2D::execute(){

	MarkBoundariesNodes();

	return AbstractAeroBoundaries::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroBoundaries_2D::MarkBoundariesNodes(){

	//Get the boundary node ids
	BoundaryOperator2D bnd_op(m_mesh);
	std::vector<TCellID> bnd_node_ids;
	bnd_op.getBoundaryNodes(bnd_node_ids);

	// Initialise une marque sur les noeuds du bord
	for (auto n_id:bnd_node_ids){
		Node n = m_mesh->get<Node>(n_id);
		m_mesh->mark(n, m_markBoundaryNodes);
	}

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroBoundaries_2D::ColoriageBordsConnexes(){

	int color = 0; //Default value is 0
	int markTreated = m_mesh->newMark<Node>();

	for (auto n_id:m_mesh->nodes())
	{
		// Si un noeud est marqué sur le bord et qu'il n'est pas encore traité
		if (m_mesh->isMarked<Node>(n_id, m_markBoundaryNodes) &&
		    !m_mesh->isMarked<Node>(n_id, markTreated)){
			Node n = m_mesh->get<Node>(n_id);

			// Nouveau bord, nouvelle couleur
			color++;
			m_mesh->mark(n, markTreated);
			(*m_var_color_bords)[n_id] = color;

			// Propagation à tous les noeuds connexes qui sont sur le bord
			std::vector<Node> next;
			next.push_back(n);

			while (!next.empty()) {
				// On récupère un noeud de la liste de noeuds next à traiter
				Node current_node = next.back();
				next.pop_back();

				// On récupère les noeuds adjacents au noeud traité
				std::vector<Edge> adjacent_edges = current_node.get<Edge>() ;
				std::vector<Node> adjacent_nodes;
				for (auto e:adjacent_edges){
					TCellID ne_id = e.getOppositeNodeId(current_node);
					Node ne = m_mesh->get<Node>(ne_id);
					adjacent_nodes.push_back(ne);
				}

				for (auto n_adj: adjacent_nodes) {
					TCellID n_adj_id = n_adj.id();
					if(m_mesh->isMarked<Node>(n_adj_id, m_markBoundaryNodes) &&
					    !m_mesh->isMarked<Node>(n_adj_id, markTreated)){
						// Si le noeud est sur le bord et qu'il n'a pas été traité
						// On met à jour sa couleur et on le marque comme traité
						m_mesh->mark(n_adj, markTreated);
						(*m_var_color_bords)[n_adj_id] = color;

						// Ajout du noeud dans la liste des noeuds à proprager
						next.push_back(n_adj);
					}
				}

			}

		}
	}

	m_nbrBords = color;
	m_nbrBordsParoi = color-1;

	if(color < 2){
		m_isImmerged = false;
	}

	m_mesh->unmarkAll<Node>(markTreated);
	m_mesh->freeMark<Node>(markTreated);

}
/*------------------------------------------------------------------------*/