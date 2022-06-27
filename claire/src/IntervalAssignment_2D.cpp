//
// Created by rochec on 27/06/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/IntervalAssignment_2D.h>
#include <gmds/ig/Mesh.h>
#include <gmds/claire/AeroExtrusion_2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

IntervalAssignment_2D::IntervalAssignment_2D(Blocking2D* ABlocking2D, ParamsAero Aparams_aero) {
	m_blocking = ABlocking2D;
	m_params_aero = Aparams_aero;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
IntervalAssignment_2D::STATUS
IntervalAssignment_2D::execute()
{

	std::map<int, std::vector<TCellID>> map_chords = ComputeChords();
	std::cout << "NOMBRE DE CORDES : " << map_chords.size() << std::endl;

	return IntervalAssignment_2D::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
std::map<int, std::vector<TCellID>>
IntervalAssignment_2D::ComputeChords(){

	std::map<int, std::vector<TCellID>> map_chords;

	int mark_isTreated = m_blocking->newMark<Edge>();
	int color_chord(0);

	for (auto e_id:m_blocking->edges()) {
		// If an edge is marked not treated
		if (!m_blocking->isMarked<Edge>(e_id, mark_isTreated))
		{
			std::vector<TCellID> chord_edges_ids;

			// Propagation à tous les noeuds connexes qui sont sur le bord
			std::vector<TCellID> next_edges;
			next_edges.push_back(e_id);

			while (!next_edges.empty()) {
				// On récupère un noeud de la liste de noeuds next à traiter
				TCellID current_e_id = next_edges.back();
				next_edges.pop_back();
				chord_edges_ids.push_back(current_e_id);

				Edge current_e = m_blocking->get<Edge>(current_e_id);
				m_blocking->mark(current_e, mark_isTreated);

				// Get the opposite edges of the edge current_e_id
				std::vector<TCellID> opp_edges = ComputeOppositeEdges(current_e_id);
				for (auto opp_e:opp_edges)
				{
					Edge opp_edge = m_blocking->get<Edge>(opp_e);
					// Add to the heap the edges not treated
					if (!m_blocking->isMarked(opp_edge, mark_isTreated)) {
						next_edges.push_back(opp_e);
					}
				}

			}

			map_chords[color_chord] = chord_edges_ids;
			color_chord++;

		}
	}

	m_blocking->unmarkAll<Node>(mark_isTreated);
	m_blocking->freeMark<Node>(mark_isTreated);

	return map_chords;

}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
std::vector<TCellID>
IntervalAssignment_2D::ComputeOppositeEdges(TCellID e_id)
{
	std::vector<TCellID> opposite_edges_id;
	Edge e = m_blocking->get<Edge>(e_id);
	std::vector<Node> e_nodes = e.get<Node>();

	std::vector<Face> e_blocks = e.get<Face>();

	for (auto f:e_blocks)
	{
		std::vector<Edge> f_edges = f.get<Edge>();
		for (auto f_edge:f_edges)
		{
			std::vector<Node> f_edge_nodes = f_edge.get<Node>();
			if (e_nodes[0].id() != f_edge_nodes[0].id() && e_nodes[0].id() != f_edge_nodes[1].id()
			    && e_nodes[1].id() != f_edge_nodes[0].id() && e_nodes[1].id() != f_edge_nodes[1].id())
			{
				opposite_edges_id.push_back(f_edge.id());
			}

		}
	}

	return opposite_edges_id;
}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void
IntervalAssignment_2D::EdgeConstraint(TCellID e_id, int N_ideal, bool hardConstraint)
{
	Edge e = m_blocking->get<Edge>(e_id) ;
	Variable<int>* var_layer_id = m_blocking->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche");

	// Default parameters
	hardConstraint = false;
	N_ideal = int(e.length()/m_params_aero.cell_size_default);

	std::vector<Face> e_faces = e.get<Face>() ;
	std::vector<Node> face_nodes = e_faces[0].get<Node>();

	// If the edge is on a border
	if (e_faces.size() == 1)
	{
		if (var_layer_id->value(face_nodes[0].id()) == 0
		    && var_layer_id->value(face_nodes[1].id()) == 0)
		{
			N_ideal = int(e.length()/m_params_aero.cell_size_default) ;
			hardConstraint = true;
		}
	}

	// If the edge is ortho to the wall, in the boundary layer
	if (var_layer_id->value(face_nodes[0].id()) == 0
	    xor var_layer_id->value(face_nodes[1].id()) == 0)
	{
		N_ideal = m_params_aero.nbrCellsInCL ;
		hardConstraint = true;
	}



	if (N_ideal <= 0)
	{
		N_ideal = 0;
	}

}
/*-------------------------------------------------------------------*/

