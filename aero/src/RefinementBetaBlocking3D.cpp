//
// Created by rochec on 08/12/23.
//
/*------------------------------------------------------------------------*/
#include <gmds/aero/RefinementBetaBlocking3D.h>
#include <gmds/aero/RefinementBetaBlock3D.h>
#include <gmds/aero/RefinementBeta.h>
#include <gmds/aero/Utils.h>
#include <gmds/ig/Mesh.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

RefinementBetaBlocking3D::RefinementBetaBlocking3D(Blocking3D *ABlocking3D,
                                                   Variable<int>* Avar_LayerID,
                                                   double Asize_first_edge)
{
	m_blocking3D = ABlocking3D;
	m_var_LayerID = Avar_LayerID;
	m_size_first_edge = Asize_first_edge;
}
/*------------------------------------------------------------------------*/
RefinementBetaBlocking3D::STATUS
RefinementBetaBlocking3D::execute()
{
	//std::cout << "1." << std::endl;
	computeSheetsToRefine();
	//std::cout << "2." << std::endl;
	for (auto sheet: m_sheets)
	{
		//std::cout << "3." << std::endl;
		sheetRefinement(&sheet);
	}

	/*
	m_sheets = m_blocking3D->computeAllSheets();

	std::vector<std::vector<TCellID>> sheetsNeedToBeRefined;
	Variable<int>* var_couche = m_blocking3D->getVariable<int, GMDS_NODE>("GMDS_Couche");

	// Compute the sheets to refine and store them in the vector sheetsNeedToBeRefined
	for (auto &sheet:m_sheets)
	{
		bool sheet_to_refine(false);
		int compteur_blocks_layer_0 = 0;

		std::vector<TCellID> sheet_blocks = m_blocking3D->computeSheetBlocks(sheet[0]) ;

		for (auto bf_id:sheet_blocks)
		{
			Blocking3D::BlockFace bf = m_blocking3D->blockFace(bf_id);

			Node n0 = bf.getNode(0);
			Node n1 = bf.getNode(1);
			Node n2 = bf.getNode(2);
			Node n3 = bf.getNode(3);

			if (var_couche->value(n0.id()) == 0
			    || var_couche->value(n1.id()) == 0
			    || var_couche->value(n2.id()) == 0
			    || var_couche->value(n3.id()) == 0)
			{
				compteur_blocks_layer_0++;
			}

			if ( (var_couche->value(n0.id()) == 0 && var_couche->value(n1.id()) == 0 )
			    || (var_couche->value(n1.id()) == 0 && var_couche->value(n2.id()) == 0 )
			    || (var_couche->value(n2.id()) == 0 && var_couche->value(n3.id()) == 0 )
			    || (var_couche->value(n3.id()) == 0 && var_couche->value(n0.id()) == 0 ))
			{
				sheet_to_refine = true;
			}

		}

		if (sheet_to_refine && compteur_blocks_layer_0 > 1){
			sheetsNeedToBeRefined.push_back(sheet);
		}

	}

	// Sheet refinement
	for (auto sheet:sheetsNeedToBeRefined)
	{
		sheetRefinement(&sheet);
	}
	 */

	return RefinementBetaBlocking3D::SUCCESS;
}
/*------------------------------------------------------------------------*/
void
RefinementBetaBlocking3D::sheetRefinement(std::vector<std::pair<TCellID,TCellID>>* sheet)
{
	for (auto pair:(*sheet))
	{
		Blocking3D::Block b = m_blocking3D->block(pair.first);
		Blocking3D::BlockFace bf = m_blocking3D->blockFace(pair.second);
		RefinementBetaBlock3D algo_ref = RefinementBetaBlock3D(&b, &bf, m_size_first_edge);
		RefinementBetaBlock3D::STATUS algo_res = algo_ref.execute();
	}

	/*
	Variable<int>* var_couche = m_blocking3D->getVariable<int, GMDS_NODE>("GMDS_Couche");

	int compteur_treated_bloc = 0;

	std::map<TCellID, bool> map_isNodeTreated;
	std::map<TCellID, bool> map_refinementNeededOnNode;

	// Initialize the maps
	TInt mark_isFaceOnSheet = m_blocking3D->newMark<Face>();

	//===============================
	// Store the faces of the sheet
	//===============================
	std::vector<TCellID> sheet_faces = m_blocking3D->computeSheetBlocks((*sheet)[0]);


	for (auto f_id:sheet_faces)
	{
		Blocking3D::BlockFace bf = m_blocking3D->blockFace(f_id);
		Node n0 = bf.getNode(0);
		Node n1 = bf.getNode(1);
		Node n2 = bf.getNode(2);
		Node n3 = bf.getNode(3);

		bool block_treated(false);

		if (var_couche->value(n0.id()) == 0
		    && var_couche->value(n1.id()) == 0)
		{
			map_refinementNeededOnNode[n0.id()] = true;
			map_refinementNeededOnNode[n1.id()] = true;
			block_treated = true;
		}
		else if (var_couche->value(n1.id()) == 0
		         && var_couche->value(n2.id()) == 0)
		{
			map_refinementNeededOnNode[n1.id()] = true;
			map_refinementNeededOnNode[n2.id()] = true;
			block_treated = true;
		}
		else if (var_couche->value(n2.id()) == 0
		         && var_couche->value(n3.id()) == 0)
		{
			map_refinementNeededOnNode[n2.id()] = true;
			map_refinementNeededOnNode[n3.id()] = true;
			block_treated = true;
		}
		else if (var_couche->value(n3.id()) == 0
		         && var_couche->value(n0.id()) == 0)
		{
			map_refinementNeededOnNode[n3.id()] = true;
			map_refinementNeededOnNode[n0.id()] = true;
			block_treated = true;
		}

		if (block_treated)
		{
			map_isNodeTreated[n0.id()] = true;
			map_isNodeTreated[n1.id()] = true;
			map_isNodeTreated[n2.id()] = true;
			map_isNodeTreated[n3.id()] = true;
			compteur_treated_bloc++;
		}

	}


	while(compteur_treated_bloc != sheet_faces.size())
	{
		for (auto f_id:sheet_faces)
		{
			Blocking3D::BlockFace bf = m_blocking3D->blockFace(f_id);
			Node n0 = bf.getNode(0);
			Node n1 = bf.getNode(1);
			Node n2 = bf.getNode(2);
			Node n3 = bf.getNode(3);
			bool block_treated(false);

			if (!map_isNodeTreated[n0.id()]
			    || !map_isNodeTreated[n1.id()]
			    || !map_isNodeTreated[n2.id()]
			    || !map_isNodeTreated[n3.id()])
			{
				if (map_refinementNeededOnNode[n0.id()]
				    && !map_refinementNeededOnNode[n1.id()]
				    && map_isNodeTreated[n1.id()])
				{
					map_refinementNeededOnNode[n3.id()] = true;
					block_treated = true;
				}
				else if (map_refinementNeededOnNode[n0.id()]
				         && !map_refinementNeededOnNode[n3.id()]
				         && map_isNodeTreated[n3.id()])
				{
					map_refinementNeededOnNode[n1.id()] = true;
					block_treated = true;
				}

				else if (map_refinementNeededOnNode[n1.id()]
				         && !map_refinementNeededOnNode[n0.id()]
				         && map_isNodeTreated[n0.id()])
				{
					map_refinementNeededOnNode[n2.id()] = true;
					block_treated = true;
				}
				else if (map_refinementNeededOnNode[n1.id()]
				         && !map_refinementNeededOnNode[n2.id()]
				         && map_isNodeTreated[n2.id()])
				{
					map_refinementNeededOnNode[n0.id()] = true;
					block_treated = true;
				}
				else if (map_refinementNeededOnNode[n2.id()]
				         && !map_refinementNeededOnNode[n1.id()]
				         && map_isNodeTreated[n1.id()])
				{
					map_refinementNeededOnNode[n3.id()] = true;
					block_treated = true;
				}
				else if (map_refinementNeededOnNode[n2.id()]
				         && !map_refinementNeededOnNode[n3.id()]
				         && map_isNodeTreated[n3.id()])
				{
					map_refinementNeededOnNode[n1.id()] = true;
					block_treated = true;
				}
				else if (map_refinementNeededOnNode[n3.id()]
				         && !map_refinementNeededOnNode[n0.id()]
				         && map_isNodeTreated[n0.id()])
				{
					map_refinementNeededOnNode[n2.id()] = true;
					block_treated = true;
				}
				else if (map_refinementNeededOnNode[n3.id()]
				         && !map_refinementNeededOnNode[n2.id()]
				         && map_isNodeTreated[n2.id()])
				{
					map_refinementNeededOnNode[n0.id()] = true;
					block_treated = true;
				}


				if (block_treated)
				{
					map_isNodeTreated[n0.id()] = true;
					map_isNodeTreated[n1.id()] = true;
					map_isNodeTreated[n2.id()] = true;
					map_isNodeTreated[n3.id()] = true;
					compteur_treated_bloc++;
				}

			}
		}
		//std::cout << "------------" << std::endl;
		//std::cout << "Compteur : " << compteur_treated_bloc << " sur " << m_map_chords[ind_chord].size() << std::endl;
	}

	 */

	// Refinement
	//for (auto b:m_blocking->allBlocks())
	/*
	for (auto f_id:sheet_faces)
	{
		Blocking3D::BlockFace bf = m_blocking3D->blockFace(f_id) ;
		int Nx = bf.getNbDiscretizationI();
		int Ny = bf.getNbDiscretizationJ();

		Node n0 = bf.getNode(0);
		Node n1 = bf.getNode(1);
		Node n2 = bf.getNode(2);
		Node n3 = bf.getNode(3);

		if (map_refinementNeededOnNode[n0.id()]
		    && map_refinementNeededOnNode[n1.id()])
		{
			for (int i = 0; i < Nx; i++) {
				std::vector<math::Point> Points;
				for (int j = 0; j < Ny; j++) {
					Points.push_back(bf(i, j).point());
				}
				RefinementBeta ref(Points, m_params_aero.edge_size_first_ortho_wall);
				ref.execute();
				Points = ref.GetNewPositions();

				for (int j = 1; j < Ny - 1; j++) {
					bf(i, j).setPoint({Points[j]});
				}
			}
		}

		if (map_refinementNeededOnNode[n2.id()]
		    && map_refinementNeededOnNode[n3.id()])
		{
			for (int i = 0; i < Nx; i++) {
				std::vector<math::Point> Points;
				for (int j = Ny-1; j >= 0; j--) {
					Points.push_back(bf(i, j).point());
				}
				RefinementBeta ref(Points, m_params_aero.edge_size_first_ortho_wall);
				ref.execute();
				Points = ref.GetNewPositions();

				for (int j = 1; j < Ny - 1; j++) {
					bf(i, j).setPoint({Points[Ny-1-j]});
				}
			}
		}

		if (map_refinementNeededOnNode[n0.id()]
		    && map_refinementNeededOnNode[n3.id()])
		{
			for (int j = 0; j < Ny; j++) {
				std::vector<math::Point> Points;
				for (int i = 0; i < Nx; i++) {
					Points.push_back(bf(i, j).point());
				}
				RefinementBeta ref(Points, m_params_aero.edge_size_first_ortho_wall);
				ref.execute();
				Points = ref.GetNewPositions();

				for (int i = 1; i < Nx - 1; i++) {
					bf(i, j).setPoint({Points[i]});
				}
			}
		}

		if (map_refinementNeededOnNode[n1.id()]
		    && map_refinementNeededOnNode[n2.id()])
		{
			for (int j = 0; j < Ny; j++) {
				std::vector<math::Point> Points;
				for (int i = Nx-1; i >= 0; i--) {
					Points.push_back(bf(i, j).point());
				}
				RefinementBeta ref(Points, m_params_aero.edge_size_first_ortho_wall);
				ref.execute();
				Points = ref.GetNewPositions();

				for (int i = 1; i < Nx - 1; i++) {
					bf(i, j).setPoint({Points[Nx-1-i]});
				}
			}
		}

	}
	 */

	// Free marks
	/*
	m_blocking3D->unmarkAll<Face>(mark_isFaceOnSheet);
	m_blocking3D->freeMark<Face>(mark_isFaceOnSheet);
	 */

}
/*------------------------------------------------------------------------*/
void
RefinementBetaBlocking3D::computeSheetsToRefine()
{
	m_sheets.clear();

	TInt mark_isTreated = m_blocking3D->newMark<Face>();
	//std::cout << "1.1." << std::endl;
	for (auto bf_id:m_blocking3D->faces())
	{
		//std::cout << "1.2." << std::endl;
		Face bf = m_blocking3D->get<Face>(bf_id);
		std::vector<Node> bf_nodes = bf.get<Node>() ;
		//std::cout << "1.3." << std::endl;
		//std::cout << m_var_LayerID->value(bf_nodes[0].id()) << std::endl;
		if (m_var_LayerID->value(bf_nodes[0].id()) == 0
		    && m_var_LayerID->value(bf_nodes[1].id()) == 0
		    && m_var_LayerID->value(bf_nodes[2].id()) == 0
		    && m_var_LayerID->value(bf_nodes[3].id()) == 0
		    && !m_blocking3D->isMarked(bf, mark_isTreated))	// If the block face is on the geometry and not already treated
		{
			//std::cout << "1.5." << std::endl;
			std::vector<Region> bf_regions = bf.get<Region>() ;
			std::vector<std::pair<TCellID,TCellID>> sheet = computeOneSheet(bf_regions[0].id(), bf_id);
			std::cout << "Sheet size: " << sheet.size() << std::endl;
			m_sheets.push_back(sheet);
			// Mark all the faces of the sheet as treated
			for (auto couple:sheet)
			{
				m_blocking3D->mark(m_blocking3D->get<Face>(couple.second), mark_isTreated);
			}
		}
	}
	m_blocking3D->unmarkAll<Face>(mark_isTreated);
	m_blocking3D->freeMark<Face>(mark_isTreated);

}
/*------------------------------------------------------------------------*/
std::vector<std::pair<TCellID,TCellID>>
RefinementBetaBlocking3D::computeOneSheet(TCellID b_id, TCellID bf_id)
{
	std::vector<std::pair<TCellID,TCellID>> sheet;
	std::pair<TCellID,TCellID> p_init(b_id,bf_id);
	//sheet.push_back(p_init);
	std::vector<std::pair<TCellID,TCellID>> next_blocks;
	next_blocks.push_back(p_init);

	TInt mark_isTreated = m_blocking3D->newMark<Face>();

	while (!next_blocks.empty())
	{
		// Get the last block id of the list of blocks to treat
		std::pair<TCellID,TCellID> current_b_id = next_blocks.back();
		next_blocks.pop_back();		// Remove the block from the list
		sheet.push_back(current_b_id);
		m_blocking3D->mark(m_blocking3D->get<Face>(current_b_id.second), mark_isTreated);

		// 1. Get the 4 faces of the blocks to propagate, and the edge connected to the face of both blocks.
		std::vector<std::pair<TCellID,TCellID>> adj_faces = compute4adjacentFacesandEdgestoFaceinBlock(current_b_id.first, current_b_id.second);

		// 2. Get the two adjacent blocks to each face. If only one: do nothing. If two: select the next one.
		for (auto f_id:adj_faces)
		{
			// We put adjacent cells in the heap
			Face f = m_blocking3D->get<Face>(f_id.first);
			Edge e = m_blocking3D->get<Edge>(f_id.second);
			std::vector<Region> f_regions = f.get<Region>();
			std::pair<TCellID,TCellID> next_pair;
			if (f_regions.size()==2)	// Then we propagate to the adjacent block if not already treated
			{
				if (f_regions[0].id() == current_b_id.first)
				{
					Blocking3D::Block b = m_blocking3D->block(f_regions[1].id());
					next_pair.first = b.id();
					std::vector<TCellID> candidate_faces = b.getEdgeFaces(f_id.second);
					if (candidate_faces[0] == f_id.first)
					{
						next_pair.second = candidate_faces[1];
					}
					else
					{
						next_pair.second = candidate_faces[0];
					}
					if ( !m_blocking3D->isMarked(m_blocking3D->get<Face>(next_pair.second), mark_isTreated) )
					{
						next_blocks.push_back(next_pair);
					}
				}
				else
				{
					Blocking3D::Block b = m_blocking3D->block(f_regions[0].id());
					next_pair.first = b.id();
					std::vector<TCellID> candidate_faces = b.getEdgeFaces(f_id.second);
					if (candidate_faces[0] == f_id.first)
					{
						next_pair.second = candidate_faces[1];
					}
					else
					{
						next_pair.second = candidate_faces[0];
					}
					if ( !m_blocking3D->isMarked(m_blocking3D->get<Face>(next_pair.second), mark_isTreated) )
					{
						next_blocks.push_back(next_pair);
					}

				}
			}
		}

	}

	m_blocking3D->unmarkAll<Face>(mark_isTreated);
	m_blocking3D->freeMark<Face>(mark_isTreated);

	return sheet;
}
/*------------------------------------------------------------------------*/
std::vector<std::pair<TCellID,TCellID>>
RefinementBetaBlocking3D::compute4adjacentFacesandEdgestoFaceinBlock(TCellID b_id, TCellID bf_id)
{
	std::vector<std::pair<TCellID,TCellID>> adj_faces_by_edges;

	Blocking3D::Block b = m_blocking3D->block(b_id);
	Face f_i0 = b.getFace(0,3,4,7);
	Face f_imax = b.getFace(1,2,5,6);
	Face f_j0 = b.getFace(0,1,4,5);
	Face f_jmax = b.getFace(2,3,6,7);
	Face f_k0 = b.getFace(0,1,2,3);
	Face f_kmax = b.getFace(4,5,6,7);

	if (bf_id == f_i0.id())
	{
		adj_faces_by_edges.push_back(std::pair(f_j0.id(), b.getEdge(0,4).id()));
		adj_faces_by_edges.push_back(std::pair(f_jmax.id(), b.getEdge(3,7).id()));
		adj_faces_by_edges.push_back(std::pair(f_k0.id(), b.getEdge(0,3).id()));
		adj_faces_by_edges.push_back(std::pair(f_kmax.id(), b.getEdge(4,7).id()));
	}
	else if (bf_id == f_imax.id())
	{
		adj_faces_by_edges.push_back(std::pair(f_j0.id(), b.getEdge(1,5).id()));
		adj_faces_by_edges.push_back(std::pair(f_jmax.id(), b.getEdge(2,6).id()));
		adj_faces_by_edges.push_back(std::pair(f_k0.id(), b.getEdge(1,2).id()));
		adj_faces_by_edges.push_back(std::pair(f_kmax.id(), b.getEdge(5,6).id()));
	}
	else if (bf_id == f_j0.id())
	{
		adj_faces_by_edges.push_back(std::pair(f_i0.id(), b.getEdge(0,4).id()));
		adj_faces_by_edges.push_back(std::pair(f_imax.id(), b.getEdge(1,5).id()));
		adj_faces_by_edges.push_back(std::pair(f_k0.id(), b.getEdge(0,1).id()));
		adj_faces_by_edges.push_back(std::pair(f_kmax.id(), b.getEdge(4,5).id()));
	}
	else if (bf_id == f_jmax.id())
	{
		adj_faces_by_edges.push_back(std::pair(f_i0.id(), b.getEdge(3,7).id()));
		adj_faces_by_edges.push_back(std::pair(f_imax.id(), b.getEdge(2,6).id()));
		adj_faces_by_edges.push_back(std::pair(f_k0.id(), b.getEdge(2,3).id()));
		adj_faces_by_edges.push_back(std::pair(f_kmax.id(), b.getEdge(6,7).id()));
	}
	else if (bf_id == f_k0.id())
	{
		adj_faces_by_edges.push_back(std::pair(f_i0.id(), b.getEdge(0,3).id()));
		adj_faces_by_edges.push_back(std::pair(f_imax.id(), b.getEdge(1,2).id()));
		adj_faces_by_edges.push_back(std::pair(f_j0.id(), b.getEdge(0,1).id()));
		adj_faces_by_edges.push_back(std::pair(f_jmax.id(), b.getEdge(2,3).id()));
	}
	else if (bf_id == f_kmax.id())
	{
		adj_faces_by_edges.push_back(std::pair(f_i0.id(), b.getEdge(4,7).id()));
		adj_faces_by_edges.push_back(std::pair(f_imax.id(), b.getEdge(5,6).id()));
		adj_faces_by_edges.push_back(std::pair(f_j0.id(), b.getEdge(4,5).id()));
		adj_faces_by_edges.push_back(std::pair(f_jmax.id(), b.getEdge(6,7).id()));
	}

	return adj_faces_by_edges;
}
/*------------------------------------------------------------------------*/