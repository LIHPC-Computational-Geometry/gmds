//
// Created by rochec on 23/03/2023.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/IntervalAssignment_3D.h>
#include <gmds/claire/Utils.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

IntervalAssignment_3D::IntervalAssignment_3D(Blocking3D* ABlocking3D,
                                             Blocking3D* ACtrlPnts3D,
                                             ParamsAero& Aparams_aero,
                                             Variable<int>* AedgesDiscretization) {
	m_Blocking3D = ABlocking3D;
	m_CtrlPts3D = ACtrlPnts3D;
	m_params_aero = Aparams_aero;
	m_edgesDiscretization = AedgesDiscretization;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
IntervalAssignment_3D::STATUS
IntervalAssignment_3D::execute()
{
	std::map<int, std::vector<TCellID>> map_sheets = ComputeSheets();

	for (auto const &sheet:map_sheets)
	{
		int Nbr_cells = ComputeSheetDiscretization(sheet.second);
		for (auto e_id:sheet.second)
		{
			m_edgesDiscretization->set(e_id, Nbr_cells);
		}
	}

	// Set the discretization on each edge of each block
	for (auto bloc: m_Blocking3D->allBlocks())
	{
		Edge e_i = bloc.getEdgeI();
		Edge e_j = bloc.getEdgeJ();
		Edge e_k = bloc.getEdgeK();
		bloc.setNbDiscretizationI(m_edgesDiscretization->value(e_i.id())+1);
		bloc.setNbDiscretizationJ(m_edgesDiscretization->value(e_j.id())+1);
		bloc.setNbDiscretizationK(m_edgesDiscretization->value(e_k.id())+1);
	}

	int node_compteur(0);
	for (auto bloc:m_Blocking3D->allBlocks())
	{
		node_compteur += (bloc.getNbDiscretizationI()-2)*(bloc.getNbDiscretizationJ()-2)*(bloc.getNbDiscretizationK()-2) ;
	}
	std::cout << "Underestimated predicted approximative number of final nodes: " << node_compteur << std::endl;

	return IntervalAssignment_3D::SUCCESS;
}
/*------------------------------------------------------------------------*/
std::vector<TCellID>
IntervalAssignment_3D::ComputeSingleSheet(TCellID e_id)
{
	std::vector<TCellID> sheet;

	TInt mark_isTreated = m_Blocking3D->newMark<Edge>();

	// Propagation à tous les noeuds connexes qui sont sur le bord
	std::vector<TCellID> next_edges;
	next_edges.push_back(e_id);

	while (!next_edges.empty())
	{
		// On récupère un noeud de la liste de noeuds next à traiter
		TCellID current_e_id = next_edges.back();
		next_edges.pop_back();
		sheet.push_back(current_e_id);

		Edge current_e = m_Blocking3D->get<Edge>(current_e_id);
		m_Blocking3D->mark(current_e, mark_isTreated);	// Mark the current edge as treated

		// Get the opposite edges of the edge current_e
		std::vector<Face> current_e_faces = current_e.get<Face>();
		for (auto const &f:current_e_faces)
		{
			Edge e_opp = math::Utils::oppositeEdgeInFace(m_Blocking3D, current_e_id, f.id());
			if (!m_Blocking3D->isMarked(e_opp, mark_isTreated))
			{
				next_edges.push_back(e_opp.id());
			}
		}

	}

	m_Blocking3D->unmarkAll<Edge>(mark_isTreated);
	m_Blocking3D->freeMark<Edge>(mark_isTreated);

	return sheet;
}
/*------------------------------------------------------------------------*/
std::map<int, std::vector<TCellID>>
IntervalAssignment_3D::ComputeSheets(){

	std::map<int, std::vector<TCellID>> map_chords;

	TInt mark_isTreated = m_Blocking3D->newMark<Edge>();
	int color_sheet(0);

	for (auto e_id: m_Blocking3D->edges()) {
		// If an edge is marked not treated
		if (!m_Blocking3D->isMarked<Edge>(e_id, mark_isTreated))
		{
			std::vector<TCellID> sheet_edges_ids = ComputeSingleSheet(e_id);
			for (auto sheet_edge_id:sheet_edges_ids)
			{
				Edge sheet_edge = m_Blocking3D->get<Edge>(sheet_edge_id);
				m_Blocking3D->mark(sheet_edge, mark_isTreated);
			}

			map_chords[color_sheet] = sheet_edges_ids;
			color_sheet++;

		}
	}

	m_Blocking3D->unmarkAll<Edge>(mark_isTreated);
	m_Blocking3D->freeMark<Edge>(mark_isTreated);

	return map_chords;

}
/*------------------------------------------------------------------------*/
void
IntervalAssignment_3D::EdgeConstraint(TCellID e_id, int &N_ideal, bool &hardConstraint)
{

	Edge e = m_Blocking3D->get<Edge>(e_id) ;
	Variable<int>* var_layer_id = m_Blocking3D->getOrCreateVariable<int, GMDS_NODE>("GMDS_Couche");

	// Default parameters
	hardConstraint = false;
	double edge_length = BezierEdgeLength(e);
	N_ideal = int(edge_length/m_params_aero.edge_size_default)+1;
	//N_ideal = int(e.length()/m_params_aero.edge_size_default)+1;

	std::vector<Face> e_faces = e.get<Face>() ;
	std::vector<Node> e_nodes = e.get<Node>() ;

	// If the edge is on the geometry surface (on the layer 0)
	if (var_layer_id->value(e_nodes[0].id()) == 0
	    && var_layer_id->value(e_nodes[1].id()) == 0)
	{
		//N_ideal = int(e.length()/m_params_aero.edge_size_wall)+1 ;
		N_ideal = int(edge_length/m_params_aero.edge_size_wall)+1 ;
		//hardConstraint = true;
	}

	// If the edge is ortho to the wall, in the boundary layer
	if ((var_layer_id->value(e_nodes[0].id()) == 0)
	    ^ (var_layer_id->value(e_nodes[1].id()) == 0))
	{
		N_ideal = m_params_aero.nbrCellsInCL;
		hardConstraint = true;
	}

	// To avoid to have a negative discretization on a sheet
	if (N_ideal <= 0)
	{
		N_ideal = 1;
	}

}
/*------------------------------------------------------------------------*/
int
IntervalAssignment_3D::ComputeSheetDiscretization(const std::vector<TCellID>& sheet)
{
	int Nbr_cells;

	bool sheet_hardConstrained(false);
	double sum_num(0.0);

	for(auto e_id:sheet)
	{
		int N_ideal;
		bool hardConstraint;
		EdgeConstraint(e_id, N_ideal, hardConstraint);
		if (hardConstraint)
		{
			Nbr_cells = N_ideal;
			sheet_hardConstrained = true;
		}
		else
		{
			sum_num += N_ideal;
		}
	}

	if (!sheet_hardConstrained)
	{
		Nbr_cells = int(1.0*sum_num/double(sheet.size()))+1;
	}

	return Nbr_cells;
}
/*------------------------------------------------------------------------*/
double
IntervalAssignment_3D::BezierEdgeLength(const Edge Ae)
{
	/*
	std::vector<Region> e_blocks = Ae.get<Region>();

	Blocking3D::Block b = m_Blocking3D->block(e_blocks[0].id());
	Blocking3D::Block b_ctrlpts = m_CtrlPts3D->block(e_blocks[0].id());

	std::vector<math::Point> ctrlpts(b_ctrlpts.getNbDiscretizationI());
	int nb_I = b_ctrlpts.getNbDiscretizationI();
	int nb_J = b_ctrlpts.getNbDiscretizationJ();
	int nb_K = b_ctrlpts.getNbDiscretizationK();

	if ( (b_ctrlpts.getNode(0).id() == e_nodes[0].id() && b_ctrlpts.getNode(1).id() == e_nodes[1].id())
	    || (b_ctrlpts.getNode(0).id() == e_nodes[1].id() && b_ctrlpts.getNode(1).id() == e_nodes[0].id()))
	{
		for (int i=0;i<nb_I;i++)
		{
			ctrlpts[i] = b_ctrlpts(i,0,0).point();
		}
	}
	else if ( (b_ctrlpts.getNode(2).id() == e_nodes[0].id() && b_ctrlpts.getNode(3).id() == e_nodes[1].id())
	    || (b_ctrlpts.getNode(2).id() == e_nodes[1].id() && b_ctrlpts.getNode(3).id() == e_nodes[0].id()))
	{
		for (int i=0;i<nb_I;i++)
		{
			ctrlpts[i] = b_ctrlpts(i,nb_J-1,0).point();
		}
	}
	else if ( (b_ctrlpts.getNode(7).id() == e_nodes[0].id() && b_ctrlpts.getNode(6).id() == e_nodes[1].id())
	    || (b_ctrlpts.getNode(7).id() == e_nodes[1].id() && b_ctrlpts.getNode(6).id() == e_nodes[0].id()))
	{
		for (int i=0;i<nb_I;i++)
		{
			ctrlpts[i] = b_ctrlpts(i,nb_J-1,nb_K-1).point();
		}
	}
	else if ( (b_ctrlpts.getNode(4).id() == e_nodes[0].id() && b_ctrlpts.getNode(5).id() == e_nodes[1].id())
	         || (b_ctrlpts.getNode(4).id() == e_nodes[1].id() && b_ctrlpts.getNode(5).id() == e_nodes[0].id()))
	{
		for (int i=0;i<nb_I;i++)
		{
			ctrlpts[i] = b_ctrlpts(i,0,nb_K-1).point();
		}
	}
	else if ( (b_ctrlpts.getNode(0).id() == e_nodes[0].id() && b_ctrlpts.getNode(3).id() == e_nodes[1].id())
	         || (b_ctrlpts.getNode(0).id() == e_nodes[1].id() && b_ctrlpts.getNode(3).id() == e_nodes[0].id()))
	{
		for (int j=0;j<nb_J;j++)
		{
			ctrlpts[j] = b_ctrlpts(0,j,0).point();
		}
	}
	else if ( (b_ctrlpts.getNode(1).id() == e_nodes[0].id() && b_ctrlpts.getNode(2).id() == e_nodes[1].id())
	         || (b_ctrlpts.getNode(1).id() == e_nodes[1].id() && b_ctrlpts.getNode(2).id() == e_nodes[0].id()))
	{
		for (int j=0;j<nb_J;j++)
		{
			ctrlpts[j] = b_ctrlpts(nb_I-1,j,0).point();
		}
	}
	else if ( (b_ctrlpts.getNode(5).id() == e_nodes[0].id() && b_ctrlpts.getNode(6).id() == e_nodes[1].id())
	         || (b_ctrlpts.getNode(5).id() == e_nodes[1].id() && b_ctrlpts.getNode(6).id() == e_nodes[0].id()))
	{
		for (int j=0;j<nb_J;j++)
		{
			ctrlpts[j] = b_ctrlpts(nb_I-1,j,nb_K-1).point();
		}
	}
	else if ( (b_ctrlpts.getNode(4).id() == e_nodes[0].id() && b_ctrlpts.getNode(7).id() == e_nodes[1].id())
	         || (b_ctrlpts.getNode(4).id() == e_nodes[1].id() && b_ctrlpts.getNode(7).id() == e_nodes[0].id()))
	{
		for (int j=0;j<nb_J;j++)
		{
			ctrlpts[j] = b_ctrlpts(0,j,nb_K-1).point();
		}
	}
	else if ( (b_ctrlpts.getNode(0).id() == e_nodes[0].id() && b_ctrlpts.getNode(4).id() == e_nodes[1].id())
	         || (b_ctrlpts.getNode(0).id() == e_nodes[1].id() && b_ctrlpts.getNode(4).id() == e_nodes[0].id()))
	{
		for (int k=0;k<nb_K;k++)
		{
			ctrlpts[k] = b_ctrlpts(0,0,k).point();
		}
	}
	else if ( (b_ctrlpts.getNode(1).id() == e_nodes[0].id() && b_ctrlpts.getNode(5).id() == e_nodes[1].id())
	         || (b_ctrlpts.getNode(1).id() == e_nodes[1].id() && b_ctrlpts.getNode(5).id() == e_nodes[0].id()))
	{
		for (int k=0;k<nb_K;k++)
		{
			ctrlpts[k] = b_ctrlpts(nb_I-1,0,k).point();
		}
	}
	else if ( (b_ctrlpts.getNode(2).id() == e_nodes[0].id() && b_ctrlpts.getNode(6).id() == e_nodes[1].id())
	         || (b_ctrlpts.getNode(2).id() == e_nodes[1].id() && b_ctrlpts.getNode(6).id() == e_nodes[0].id()))
	{
		for (int k=0;k<nb_K;k++)
		{
			ctrlpts[k] = b_ctrlpts(nb_I-1,nb_J-1,k).point();
		}
	}
	else if ( (b_ctrlpts.getNode(3).id() == e_nodes[0].id() && b_ctrlpts.getNode(7).id() == e_nodes[1].id())
	         || (b_ctrlpts.getNode(3).id() == e_nodes[1].id() && b_ctrlpts.getNode(7).id() == e_nodes[0].id()))
	{
		for (int k=0;k<nb_K;k++)
		{
			ctrlpts[k] = b_ctrlpts(0,nb_J-1,k).point();
		}
	}
	else
	{
		throw GMDSException("ERROR in IntervalAssignment_3D: the edge does not exist");
	}
	 */

	std::vector<math::Point> ctrlpts = m_CtrlPts3D->getEdgeNodesPoints(Ae.id());
	math::BezierCurve bc = math::BezierCurve(ctrlpts);
	return math::Utils::lengthBezierCurve(&bc);

}
/*------------------------------------------------------------------------*/