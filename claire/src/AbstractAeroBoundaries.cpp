//
// Created by rochec on 23/03/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AbstractAeroBoundaries.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
AbstractAeroBoundaries::AbstractAeroBoundaries(Mesh *AMesh) : m_mesh(AMesh), m_isImmerged(true)
{
	m_markNodesParoi = m_mesh->newMark<gmds::Node>();
	m_markNodesAmont = m_mesh->newMark<gmds::Node>();
	m_markBoundaryNodes = m_mesh->newMark<Node>();
	m_var_color_bords = m_mesh->getOrCreateVariable<int, GMDS_NODE>("COLOR_BORDS");
}
/*------------------------------------------------------------------------*/

AbstractAeroBoundaries::~AbstractAeroBoundaries()
{
	m_mesh->unmarkAll<Node>(m_markBoundaryNodes);
	m_mesh->freeMark<Node>(m_markBoundaryNodes);
	m_mesh->unmarkAll<Node>(m_markNodesParoi);
	m_mesh->freeMark<Node>(m_markNodesParoi);
	m_mesh->unmarkAll<Node>(m_markNodesAmont);
	m_mesh->freeMark<Node>(m_markNodesAmont);
}
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
AbstractAeroBoundaries::STATUS
AbstractAeroBoundaries::execute()
{

	MarkBoundariesNodes();
	ColoriageBordsConnexes();
	WhichColorIsAmont();
	MarkAmontAndParoiNodes();

	return AbstractAeroBoundaries::SUCCESS;
}
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
bool
AbstractAeroBoundaries::isBnd(TCellID n_id)
{
	return m_mesh->isMarked<Node>(n_id, m_markBoundaryNodes);
}
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
bool
AbstractAeroBoundaries::isAmont(TCellID n_id)
{
	return m_mesh->isMarked<Node>(n_id, m_markNodesAmont);
}
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
bool
AbstractAeroBoundaries::isParoi(TCellID n_id)
{
	return m_mesh->isMarked<Node>(n_id, m_markNodesParoi);
}
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
int
AbstractAeroBoundaries::getMarkBnd()
{
	return m_markBoundaryNodes;
}
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
int
AbstractAeroBoundaries::getMarkAmont()
{
	return m_markNodesAmont;
}
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
int
AbstractAeroBoundaries::getMarkParoi()
{
	return m_markNodesParoi;
}
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
bool
AbstractAeroBoundaries::isImmerged()
{
	return m_isImmerged;
}
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
int
AbstractAeroBoundaries::getNbrBords()
{
	return m_nbrBords;
}
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
int
AbstractAeroBoundaries::getColorAmont()
{
	return m_color_Amont;
}
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
int
AbstractAeroBoundaries::getNodeColor(TCellID n_id)
{
	return m_var_color_bords->value(n_id);
}
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
TCellID
AbstractAeroBoundaries::PointArret(int color)
{

	TCellID n_arret_id(NullID);
	double x_min = std::numeric_limits<double>::max();

	for (auto n_it = m_mesh->nodes_begin(); n_it != m_mesh->nodes_end(); ++n_it) {
		TCellID n_id = *n_it;
		Node n = m_mesh->get<Node>(n_id);
		math::Point p = n.point();
		if (m_mesh->isMarked<Node>(n_id, m_markNodesParoi) && m_var_color_bords->value(n_id) == color && p.X() < x_min) {
			n_arret_id = n_id;
			x_min = p.X();
		}
	}
	return n_arret_id;
}
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
void
AbstractAeroBoundaries::ColoriageBordsConnexes()
{

	int color = 0;     // Default value is 0
	int markTreated = m_mesh->newMark<Node>();

	for (auto n_id : m_mesh->nodes()) {
		// Si un noeud est marqué sur le bord et qu'il n'est pas encore traité
		if (m_mesh->isMarked<Node>(n_id, m_markBoundaryNodes) && !m_mesh->isMarked<Node>(n_id, markTreated)) {
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
				std::vector<Edge> adjacent_edges = current_node.get<Edge>();
				std::vector<Node> adjacent_nodes;
				for (auto e : adjacent_edges) {
					TCellID ne_id = e.getOppositeNodeId(current_node);
					Node ne = m_mesh->get<Node>(ne_id);
					adjacent_nodes.push_back(ne);
				}

				for (auto n_adj : adjacent_nodes) {
					TCellID n_adj_id = n_adj.id();
					if (m_mesh->isMarked<Node>(n_adj_id, m_markBoundaryNodes) && !m_mesh->isMarked<Node>(n_adj_id, markTreated)) {
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

	if (color < 2) {
		m_isImmerged = false;
	}

	m_mesh->unmarkAll<Node>(markTreated);
	m_mesh->freeMark<Node>(markTreated);
}
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
void
AbstractAeroBoundaries::MarkAmontAndParoiNodes()
{
	// On marque les fronts paroi et extérieur
	for (auto n_id : m_bnd_nodes_ids) {
		int couleur = m_var_color_bords->value(n_id);
		if (couleur == m_color_Amont) {
			m_mesh->mark<Node>(n_id, m_markNodesAmont);
		}
		else {
			m_mesh->mark<Node>(n_id, m_markNodesParoi);
		}
	}
}
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
TCellID
AbstractAeroBoundaries::ClosestNodeOnBnd(int color, math::Point p)
{
	TCellID closest_node_id;
	double min(std::numeric_limits<double>::max());

	for (auto n_id : m_mesh->nodes()) {
		Node n = m_mesh->get<Node>(n_id);
		math::Vector3d vec = p - n.point();
		if (m_var_color_bords->value(n_id) == color && vec.norm() < min) {
			min = vec.norm();
			closest_node_id = n_id;
		}
	}
	return closest_node_id;
}
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
TCellID
AbstractAeroBoundaries::RandomNodeOnBnd(int color)
{
	for (auto n_id : m_mesh->nodes()) {
		if (m_var_color_bords->value(n_id) == color) {
			return n_id;
		}
	}
	return NullID;
}
/*------------------------------------------------------------------------*/