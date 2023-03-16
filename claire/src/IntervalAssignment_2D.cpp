//
// Created by rochec on 27/06/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/IntervalAssignment_2D.h>
#include <gmds/ig/Mesh.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

IntervalAssignment_2D::IntervalAssignment_2D(Blocking2D* ABlocking2D, ParamsAero& Aparams_aero) {
	m_blocking = ABlocking2D;
	m_params_aero = Aparams_aero;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
IntervalAssignment_2D::STATUS
IntervalAssignment_2D::execute()
{

	std::map<int, std::vector<TCellID>> map_chords = ComputeChords();
	//std::cout << "NOMBRE DE CORDES : " << map_chords.size() << std::endl;

	Variable<int>* var_NbrCells = m_blocking->getOrCreateVariable<int, GMDS_EDGE>("NbrCells");

	for (auto & chord:map_chords)
	{
		int Nb_cells = ComputeChordDiscretization(chord.second);
		//std::cout << "Nombre de cellules dans la corde " << chord.first << " : " << Nb_cells << std::endl;
		for (auto e_id:chord.second)
		{
			var_NbrCells->set(e_id, Nb_cells);
		}
	}

	/*
	for (auto bloc_id:m_blocking->faces())
	{
		Face bloc = m_blocking->get<Face>(bloc_id);
		Blocking2D::Block B = m_blocking->block(bloc_id);
		std::vector<Edge> bloc_edges = bloc.get<Edge>() ; // e0, e2 are subdivsion I, e1, e3 are subdivision J (according to Blocking2D class)

		B.setNbDiscretizationI(var_NbrCells->value(bloc_edges[0].id()));
		B.setNbDiscretizationJ(var_NbrCells->value(bloc_edges[1].id()));

	}
	 */
	for (auto bloc:m_blocking->allBlocks())
	{
		Edge e_i = bloc.getEdgeI();
		Edge e_j = bloc.getEdgeJ();
		bloc.setNbDiscretizationI(var_NbrCells->value(e_i.id())+1);
		bloc.setNbDiscretizationJ(var_NbrCells->value(e_j.id())+1);
	}

	m_blocking->deleteVariable(GMDS_EDGE, "NbrCells");

	return IntervalAssignment_2D::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
std::map<int, std::vector<TCellID>>
IntervalAssignment_2D::ComputeChords(){

	std::map<int, std::vector<TCellID>> map_chords;

	TInt mark_isTreated = m_blocking->newMark<Edge>();
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

	for (auto const& f:e_blocks)
	{
		std::vector<Edge> f_edges = f.get<Edge>();
		for (auto const& f_edge:f_edges)
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
IntervalAssignment_2D::EdgeConstraint(TCellID e_id, int &N_ideal, bool &hardConstraint)
{
	Edge e = m_blocking->get<Edge>(e_id) ;
	Variable<int>* var_layer_id = m_blocking->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche");

	// Default parameters
	hardConstraint = false;
	N_ideal = int(e.length()/m_params_aero.edge_size_default);

	std::vector<Face> e_faces = e.get<Face>() ;
	std::vector<Node> e_nodes = e.get<Node>() ;

	// If the edge is on a border
	if (e_faces.size() == 1)
	{
		if (var_layer_id->value(e_nodes[0].id()) == 0
		    && var_layer_id->value(e_nodes[1].id()) == 0)
		{
			N_ideal = int(e.length()/m_params_aero.edge_size_wall) ;
			hardConstraint = true;
		}
	}

	// If the edge is ortho to the wall, in the boundary layer
	if ((var_layer_id->value(e_nodes[0].id()) == 0)
	    ^ (var_layer_id->value(e_nodes[1].id()) == 0))
	{
		bool inserted_edge(false);
		std::vector<Node> f0_nodes = e_faces[0].get<Node>();
		std::vector<Node> f1_nodes = e_faces[1].get<Node>();
		if ( ( (var_layer_id->value(f0_nodes[0].id()) == 0)
		    ^ (var_layer_id->value(f0_nodes[1].id()) == 0)
		    ^ (var_layer_id->value(f0_nodes[2].id()) == 0)
		    ^ (var_layer_id->value(f0_nodes[3].id()) == 0) )
		    &&
		    ( (var_layer_id->value(f1_nodes[0].id()) == 0)
		     ^ (var_layer_id->value(f1_nodes[1].id()) == 0)
		     ^ (var_layer_id->value(f1_nodes[2].id()) == 0)
		     ^ (var_layer_id->value(f1_nodes[3].id()) == 0) ) )
		{
			inserted_edge = true;
		}

		if (!inserted_edge) {
			N_ideal = m_params_aero.nbrCellsInCL;
			hardConstraint = true;
		}
	}



	if (N_ideal <= 0)
	{
		N_ideal = 1;
	}

}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
int
IntervalAssignment_2D::ComputeChordDiscretization(std::vector<TCellID>& chord)
{
	int Nbr_cells;

	bool chord_hardConstrained(false);
	double sum_num(0.0);

	for(auto e_id:chord)
	{
		int N_ideal;
		bool hardConstraint;
		EdgeConstraint(e_id, N_ideal, hardConstraint);
		if (hardConstraint)
		{
			Nbr_cells = N_ideal;
			chord_hardConstrained = true;
		}
		else
		{
			sum_num += N_ideal;
		}
	}

	if (!chord_hardConstrained)
	{
		Nbr_cells = int(1.0*sum_num/double(chord.size()));
	}

	return Nbr_cells;
}
/*-------------------------------------------------------------------*/