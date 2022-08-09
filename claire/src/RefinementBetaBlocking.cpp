//
// Created by rochec on 08/08/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/RefinementBetaBlocking.h>
#include <gmds/claire/RefinementBeta.h>
#include <gmds/ig/Mesh.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

RefinementBetaBlocking::RefinementBetaBlocking(Blocking2D *ABlocking2D, ParamsAero Aparams_aero) {
	m_blocking = ABlocking2D;
	m_params_aero = Aparams_aero;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
RefinementBetaBlocking::STATUS
RefinementBetaBlocking::execute()
{
	m_map_chords = ComputeChords();
	std::vector<int> chords_to_refine ;
	Variable<int>* var_couche = m_blocking->getVariable<int, GMDS_NODE>("GMDS_Couche");

	// Compute the chords to refine and store them in the vector chords_to_refined
	for (auto chord:m_map_chords)
	{
		/*
		for (auto b_id:chord.second)
		{
			std::cout << "Corde " << chord.first << ", bloc d'id : " << b_id << std::endl;
			Blocking2D::Block b = m_blocking->block(b_id);
			Node n0 = b.getNode(0);
			Node n1 = b.getNode(1);
			Node n2 = b.getNode(2);
			Node n3 = b.getNode(3);

			std::cout << "Noeud " << n0.id() << ", sur la couche : " << var_couche->value(n0.id()) << std::endl;
			std::cout << "Noeud " << n1.id() << ", sur la couche : " << var_couche->value(n1.id()) << std::endl;
			std::cout << "Noeud " << n2.id() << ", sur la couche : " << var_couche->value(n2.id()) << std::endl;
			std::cout << "Noeud " << n3.id() << ", sur la couche : " << var_couche->value(n3.id()) << std::endl;
		}
		 */

		bool chord_to_refine(false);
		int compteur_blocks_layer_0 = 0;

		for (auto b_id:chord.second)
		{
			Blocking2D::Block b = m_blocking->block(b_id);

			Node n0 = b.getNode(0);
			Node n1 = b.getNode(1);
			Node n2 = b.getNode(2);
			Node n3 = b.getNode(3);
			/*
			std::cout << " ------------------------- " << std::endl;
			std::cout << "BLOCK " << b_id << std::endl;
			std::cout << "Node " << n0.id() << " , layer : " << var_couche->value(n0.id()) << std::endl;
			std::cout << "Node " << n1.id() << " , layer : " << var_couche->value(n1.id()) << std::endl;
			std::cout << "Node " << n2.id() << " , layer : " << var_couche->value(n2.id()) << std::endl;
			std::cout << "Node " << n3.id() << " , layer : " << var_couche->value(n3.id()) << std::endl;
			 */

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
				chord_to_refine = true;
			}

			/*
			if (!chord_to_refine
			    && var_couche->value(n0.id()) == 0
			    && var_couche->value(n1.id()) == 0)
			{
				chord_to_refine = true;
			}
			else if (!chord_to_refine
			    && var_couche->value(n1.id()) == 0
			    && var_couche->value(n2.id()) == 0)
			{
				chord_to_refine = true;
			}
			else if (!chord_to_refine
			         && var_couche->value(n2.id()) == 0
			         && var_couche->value(n3.id()) == 0)
			{
				chord_to_refine = true;
			}
			else if (!chord_to_refine
			         && var_couche->value(n3.id()) == 0
			         && var_couche->value(n0.id()) == 0)
			{
				chord_to_refine = true;
			}
			 */

		}

		if (chord_to_refine && compteur_blocks_layer_0 > 1){
			chords_to_refine.push_back(chord.first);
			std::cout << "Indice de corde à raffiner : " << chord.first << std::endl;
		}

	}

	// Chord refinement
	for (auto ind_chord:chords_to_refine)
	{
		ChordRefinement(ind_chord);
	}

	return RefinementBetaBlocking::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
void
RefinementBetaBlocking::ChordRefinement(int ind_chord){

	Variable<int>* var_couche = m_blocking->getVariable<int, GMDS_NODE>("GMDS_Couche");

	/*
	// Mark the bloc nodes for a chord to refine
	//int mark_isTreated = m_blocking->newMark<Node>();
	//int mark_refinementNeeded = m_blocking->newMark<Node>();
	//std::cout << "marque 1 : " << mark_isTreated << " marque 2 : " << mark_isTreated << std::endl;
	int compteur_treated_bloc = 0;

	for (auto b_id:m_map_chords[ind_chord])
	{
		Blocking2D::Block b = m_blocking->block(b_id);
		Node n0 = b.getNode(0);
		Node n1 = b.getNode(1);
		Node n2 = b.getNode(2);
		Node n3 = b.getNode(3);

		bool block_treated(false);

		if (var_couche->value(n0.id()) == 0
		    && var_couche->value(n1.id()) == 0)
		{
			m_blocking->mark(n0, mark_refinementNeeded);
			m_blocking->mark(n1, mark_refinementNeeded);
			block_treated = true;
		}
		else if (var_couche->value(n1.id()) == 0
		         && var_couche->value(n2.id()) == 0)
		{
			m_blocking->mark(n1, mark_refinementNeeded);
			m_blocking->mark(n2, mark_refinementNeeded);
			block_treated = true;
		}
		else if (var_couche->value(n2.id()) == 0
		         && var_couche->value(n3.id()) == 0)
		{
			m_blocking->mark(n2, mark_refinementNeeded);
			m_blocking->mark(n3, mark_refinementNeeded);
			block_treated = true;
		}
		else if (var_couche->value(n3.id()) == 0
		         && var_couche->value(n0.id()) == 0)
		{
			m_blocking->mark(n3, mark_refinementNeeded);
			m_blocking->mark(n0, mark_refinementNeeded);
			block_treated = true;
		}

		if (block_treated)
		{
			m_blocking->mark(n0, mark_isTreated);
			m_blocking->mark(n1, mark_isTreated);
			m_blocking->mark(n2, mark_isTreated);
			m_blocking->mark(n3, mark_isTreated);
			compteur_treated_bloc++;
		}

		std::cout << "------------------------" << std::endl;
		std::cout << "Bloc " << b.id() << std::endl;
		std::cout << "Node " << n0.id() << ", treated : " << m_blocking->isMarked(n0, mark_isTreated) << ", refined : " << m_blocking->isMarked(n0, mark_refinementNeeded) << ", layer : " << var_couche->value(n0.id()) << std::endl;
		std::cout << "Node " << n1.id() << ", treated : " << m_blocking->isMarked(n1, mark_isTreated) << ", refined : " << m_blocking->isMarked(n1, mark_refinementNeeded) << ", layer : " << var_couche->value(n1.id()) << std::endl;
		std::cout << "Node " << n2.id() << ", treated : " << m_blocking->isMarked(n2, mark_isTreated) << ", refined : " << m_blocking->isMarked(n2, mark_refinementNeeded) << ", layer : " << var_couche->value(n2.id()) << std::endl;
		std::cout << "Node " << n3.id() << ", treated : " << m_blocking->isMarked(n3, mark_isTreated) << ", refined : " << m_blocking->isMarked(n3, mark_refinementNeeded) << ", layer : " << var_couche->value(n3.id()) << std::endl;


	}

	std::cout << "===============================" << std::endl;

	while(compteur_treated_bloc != m_map_chords[ind_chord].size())
	{
		for (auto b_id:m_map_chords[ind_chord])
		{
			Blocking2D::Block b = m_blocking->block(b_id);
			Node n0 = b.getNode(0);
			Node n1 = b.getNode(1);
			Node n2 = b.getNode(2);
			Node n3 = b.getNode(3);
			bool block_treated(false);

			std::cout << "------------------------" << std::endl;
			std::cout << "Bloc " << b.id() << std::endl;
			std::cout << "Node " << n0.id() << ", treated : " << m_blocking->isMarked(n0, mark_isTreated) << ", refined : " << m_blocking->isMarked(n0, mark_refinementNeeded) << ", layer : " << var_couche->value(n0.id()) << std::endl;
			std::cout << "Node " << n1.id() << ", treated : " << m_blocking->isMarked(n1, mark_isTreated) << ", refined : " << m_blocking->isMarked(n1, mark_refinementNeeded) << ", layer : " << var_couche->value(n1.id()) << std::endl;
			std::cout << "Node " << n2.id() << ", treated : " << m_blocking->isMarked(n2, mark_isTreated) << ", refined : " << m_blocking->isMarked(n2, mark_refinementNeeded) << ", layer : " << var_couche->value(n2.id()) << std::endl;
			std::cout << "Node " << n3.id() << ", treated : " << m_blocking->isMarked(n3, mark_isTreated) << ", refined : " << m_blocking->isMarked(n3, mark_refinementNeeded) << ", layer : " << var_couche->value(n3.id()) << std::endl;


			if (!m_blocking->isMarked(n0, mark_isTreated)
			    || !m_blocking->isMarked(n1, mark_isTreated)
			    || !m_blocking->isMarked(n2, mark_isTreated)
			    || !m_blocking->isMarked(n3, mark_isTreated))
			{
				if (m_blocking->isMarked(n0, mark_refinementNeeded)
				    && !m_blocking->isMarked(n1, mark_refinementNeeded)
				    && m_blocking->isMarked(n1, mark_isTreated))
				{
					m_blocking->mark(n3, mark_refinementNeeded);
					block_treated = true;
				}
				else if (m_blocking->isMarked(n0, mark_refinementNeeded)
				         && !m_blocking->isMarked(n3, mark_refinementNeeded)
				         && m_blocking->isMarked(n3, mark_isTreated))
				{
					m_blocking->mark(n1, mark_refinementNeeded);
					block_treated = true;
				}

				else if (m_blocking->isMarked(n1, mark_refinementNeeded)
				         && !m_blocking->isMarked(n0, mark_refinementNeeded)
				         && m_blocking->isMarked(n0, mark_isTreated))
				{
					m_blocking->mark(n2, mark_refinementNeeded);
					block_treated = true;
				}
				else if (m_blocking->isMarked(n1, mark_refinementNeeded)
				         && !m_blocking->isMarked(n2, mark_refinementNeeded)
				         && m_blocking->isMarked(n2, mark_isTreated))
				{
					m_blocking->mark(n0, mark_refinementNeeded);
					block_treated = true;
				}

				else if (m_blocking->isMarked(n2, mark_refinementNeeded)
				         && !m_blocking->isMarked(n1, mark_refinementNeeded)
				         && m_blocking->isMarked(n1, mark_isTreated))
				{
					m_blocking->mark(n3, mark_refinementNeeded);
					block_treated = true;
				}
				else if (m_blocking->isMarked(n2, mark_refinementNeeded)
				         && !m_blocking->isMarked(n3, mark_refinementNeeded)
				         && m_blocking->isMarked(n3, mark_isTreated))
				{
					m_blocking->mark(n1, mark_refinementNeeded);
					block_treated = true;
				}

				else if (m_blocking->isMarked(n3, mark_refinementNeeded)
				         && !m_blocking->isMarked(n0, mark_refinementNeeded)
				         && m_blocking->isMarked(n0, mark_isTreated))
				{
					m_blocking->mark(n2, mark_refinementNeeded);
					block_treated = true;
				}
				else if (m_blocking->isMarked(n3, mark_refinementNeeded)
				         && !m_blocking->isMarked(n2, mark_refinementNeeded)
				         && m_blocking->isMarked(n2, mark_isTreated))
				{
					m_blocking->mark(n0, mark_refinementNeeded);
					block_treated = true;
				}


				if (block_treated)
				{
					m_blocking->mark(n0, mark_isTreated);
					m_blocking->mark(n1, mark_isTreated);
					m_blocking->mark(n2, mark_isTreated);
					m_blocking->mark(n3, mark_isTreated);
					compteur_treated_bloc++;
				}

			}
		}
		std::cout << "------------" << std::endl;
		std::cout << "Compteur : " << compteur_treated_bloc << " sur " << m_map_chords[ind_chord].size() << std::endl;
	}


	std::cout << "Nbr of blocs in the chord : " << m_map_chords[ind_chord].size() << std::endl;
	std::cout << "Compteur : " << compteur_treated_bloc << std::endl;
	
	
	
	// Refinement
	for (auto b:m_blocking->allBlocks())
	{
		int Nx = b.getNbDiscretizationI();
		int Ny = b.getNbDiscretizationJ();
		
		Variable<int>* var_couche = m_blocking->getVariable<int, GMDS_NODE>("GMDS_Couche");
		Node n0 = b.getNode(0);
		Node n1 = b.getNode(1);
		Node n2 = b.getNode(2);
		Node n3 = b.getNode(3);

		if (m_blocking->isMarked(n0, mark_refinementNeeded)
		    && m_blocking->isMarked(n1, mark_refinementNeeded))
		{
			for (int i = 0; i < Nx; i++) {
				std::vector<math::Point> Points;
				for (int j = 0; j < Ny; j++) {
					Points.push_back(b(i, j).point());
				}
				RefinementBeta ref(Points, m_params_aero.edge_size_first_ortho_wall);
				ref.execute();
				Points = ref.GetNewPositions();

				for (int j = 1; j < Ny - 1; j++) {
					b(i, j).setPoint({Points[j]});
				}
			}
		}

		if (m_blocking->isMarked(n2, mark_refinementNeeded)
		    && m_blocking->isMarked(n3, mark_refinementNeeded))
		{
			for (int i = 0; i < Nx; i++) {
				std::vector<math::Point> Points;
				for (int j = Ny-1; j >= 0; j--) {
					Points.push_back(b(i, j).point());
				}
				RefinementBeta ref(Points, m_params_aero.edge_size_first_ortho_wall);
				ref.execute();
				Points = ref.GetNewPositions();

				for (int j = 1; j < Ny - 1; j++) {
					b(i, j).setPoint({Points[Ny-1-j]});
				}
			}
		}

		if (m_blocking->isMarked(n0, mark_refinementNeeded)
		    && m_blocking->isMarked(n3, mark_refinementNeeded))
		{
			for (int j = 0; j < Ny; j++) {
				std::vector<math::Point> Points;
				for (int i = 0; i < Nx; i++) {
					Points.push_back(b(i, j).point());
				}
				RefinementBeta ref(Points, m_params_aero.edge_size_first_ortho_wall);
				ref.execute();
				Points = ref.GetNewPositions();

				for (int i = 1; i < Nx - 1; i++) {
					b(i, j).setPoint({Points[i]});
				}
			}
		}

		if (m_blocking->isMarked(n1, mark_refinementNeeded)
		    && m_blocking->isMarked(n2, mark_refinementNeeded))
		{
			for (int j = 0; j < Ny; j++) {
				std::vector<math::Point> Points;
				for (int i = Nx-1; i >= 0; i--) {
					Points.push_back(b(i, j).point());
				}
				RefinementBeta ref(Points, m_params_aero.edge_size_first_ortho_wall);
				ref.execute();
				Points = ref.GetNewPositions();

				for (int i = 1; i < Nx - 1; i++) {
					b(i, j).setPoint({Points[Nx-1-i]});
				}
			}
		}

	}
	
	
	// Free marks
	m_blocking->unmarkAll<Node>(mark_isTreated);
	m_blocking->freeMark<Node>(mark_isTreated);
	m_blocking->unmarkAll<Node>(mark_refinementNeeded);
	m_blocking->freeMark<Node>(mark_refinementNeeded);

	 */


	int compteur_treated_bloc = 0;
	std::map<TCellID, bool> map_isTreated;
	std::map<TCellID, bool> map_refinementNeeded;

	for (auto n_id:m_blocking->nodes())
	{
		map_isTreated[n_id] = false;
		map_refinementNeeded[n_id] = false;
	}



	for (auto b_id:m_map_chords[ind_chord])
	{
		Blocking2D::Block b = m_blocking->block(b_id);
		Node n0 = b.getNode(0);
		Node n1 = b.getNode(1);
		Node n2 = b.getNode(2);
		Node n3 = b.getNode(3);

		bool block_treated(false);

		if (var_couche->value(n0.id()) == 0
		    && var_couche->value(n1.id()) == 0)
		{
			map_refinementNeeded[n0.id()] = true;
			map_refinementNeeded[n1.id()] = true;
			block_treated = true;
		}
		else if (var_couche->value(n1.id()) == 0
		         && var_couche->value(n2.id()) == 0)
		{
			map_refinementNeeded[n1.id()] = true;
			map_refinementNeeded[n2.id()] = true;
			block_treated = true;
		}
		else if (var_couche->value(n2.id()) == 0
		         && var_couche->value(n3.id()) == 0)
		{
			map_refinementNeeded[n2.id()] = true;
			map_refinementNeeded[n3.id()] = true;
			block_treated = true;
		}
		else if (var_couche->value(n3.id()) == 0
		         && var_couche->value(n0.id()) == 0)
		{
			map_refinementNeeded[n3.id()] = true;
			map_refinementNeeded[n0.id()] = true;
			block_treated = true;
		}

		if (block_treated)
		{
			map_isTreated[n0.id()] = true;
			map_isTreated[n1.id()] = true;
			map_isTreated[n2.id()] = true;
			map_isTreated[n3.id()] = true;
			compteur_treated_bloc++;
		}

		/*
		std::cout << "------------------------" << std::endl;
		std::cout << "Bloc " << b.id() << std::endl;
		std::cout << "Node " << n0.id() << ", treated : " << m_blocking->isMarked(n0, mark_isTreated) << ", refined : " << m_blocking->isMarked(n0, mark_refinementNeeded) << ", layer : " << var_couche->value(n0.id()) << std::endl;
		std::cout << "Node " << n1.id() << ", treated : " << m_blocking->isMarked(n1, mark_isTreated) << ", refined : " << m_blocking->isMarked(n1, mark_refinementNeeded) << ", layer : " << var_couche->value(n1.id()) << std::endl;
		std::cout << "Node " << n2.id() << ", treated : " << m_blocking->isMarked(n2, mark_isTreated) << ", refined : " << m_blocking->isMarked(n2, mark_refinementNeeded) << ", layer : " << var_couche->value(n2.id()) << std::endl;
		std::cout << "Node " << n3.id() << ", treated : " << m_blocking->isMarked(n3, mark_isTreated) << ", refined : " << m_blocking->isMarked(n3, mark_refinementNeeded) << ", layer : " << var_couche->value(n3.id()) << std::endl;
		*/

	}

	std::cout << "===============================" << std::endl;

	while(compteur_treated_bloc != m_map_chords[ind_chord].size())
	{
		for (auto b_id:m_map_chords[ind_chord])
		{
			Blocking2D::Block b = m_blocking->block(b_id);
			Node n0 = b.getNode(0);
			Node n1 = b.getNode(1);
			Node n2 = b.getNode(2);
			Node n3 = b.getNode(3);
			bool block_treated(false);

			/*
			int compteur_refinementneeded = 0;
			if (map_refinementNeeded[n0.id()])
			{
				compteur_refinementneeded++;
			}
			if (map_refinementNeeded[n1.id()])
			{
				compteur_refinementneeded++;
			}
			if (map_refinementNeeded[n2.id()])
			{
				compteur_refinementneeded++;
			}
			if (map_refinementNeeded[n3.id()])
			{
				compteur_refinementneeded++;
			}
			if (compteur_refinementneeded==2)
			{
				block_treated = true;
			}
			 */

			/*
			std::cout << "------------------------" << std::endl;
			std::cout << "Bloc " << b.id() << std::endl;
			std::cout << "Node " << n0.id() << ", treated : " << m_blocking->isMarked(n0, mark_isTreated) << ", refined : " << m_blocking->isMarked(n0, mark_refinementNeeded) << ", layer : " << var_couche->value(n0.id()) << std::endl;
			std::cout << "Node " << n1.id() << ", treated : " << m_blocking->isMarked(n1, mark_isTreated) << ", refined : " << m_blocking->isMarked(n1, mark_refinementNeeded) << ", layer : " << var_couche->value(n1.id()) << std::endl;
			std::cout << "Node " << n2.id() << ", treated : " << m_blocking->isMarked(n2, mark_isTreated) << ", refined : " << m_blocking->isMarked(n2, mark_refinementNeeded) << ", layer : " << var_couche->value(n2.id()) << std::endl;
			std::cout << "Node " << n3.id() << ", treated : " << m_blocking->isMarked(n3, mark_isTreated) << ", refined : " << m_blocking->isMarked(n3, mark_refinementNeeded) << ", layer : " << var_couche->value(n3.id()) << std::endl;
			*/

			if (!map_isTreated[n0.id()]
			    || !map_isTreated[n1.id()]
			    || !map_isTreated[n2.id()]
			    || !map_isTreated[n3.id()])
			{
				if (map_refinementNeeded[n0.id()]
				    && !map_refinementNeeded[n1.id()]
				    && map_isTreated[n1.id()])
				{
					map_refinementNeeded[n3.id()] = true;
					block_treated = true;
				}
				else if (map_refinementNeeded[n0.id()]
				         && !map_refinementNeeded[n3.id()]
				         && map_isTreated[n3.id()])
				{
					map_refinementNeeded[n1.id()] = true;
					block_treated = true;
				}

				else if (map_refinementNeeded[n1.id()]
				         && !map_refinementNeeded[n0.id()]
				         && map_isTreated[n0.id()])
				{
					map_refinementNeeded[n2.id()] = true;
					block_treated = true;
				}
				else if (map_refinementNeeded[n1.id()]
				         && !map_refinementNeeded[n2.id()]
				         && map_isTreated[n2.id()])
				{
					map_refinementNeeded[n0.id()] = true;
					block_treated = true;
				}
				else if (map_refinementNeeded[n2.id()]
				         && !map_refinementNeeded[n1.id()]
				         && map_isTreated[n1.id()])
				{
					map_refinementNeeded[n3.id()] = true;
					block_treated = true;
				}
				else if (map_refinementNeeded[n2.id()]
				         && !map_refinementNeeded[n3.id()]
				         && map_isTreated[n3.id()])
				{
					map_refinementNeeded[n1.id()] = true;
					block_treated = true;
				}
				else if (map_refinementNeeded[n3.id()]
				         && !map_refinementNeeded[n0.id()]
				         && map_isTreated[n0.id()])
				{
					map_refinementNeeded[n2.id()] = true;
					block_treated = true;
				}
				else if (map_refinementNeeded[n3.id()]
				         && !map_refinementNeeded[n2.id()]
				         && map_isTreated[n2.id()])
				{
					map_refinementNeeded[n0.id()] = true;
					block_treated = true;
				}


				if (block_treated)
				{
					map_isTreated[n0.id()] = true;
					map_isTreated[n1.id()] = true;
					map_isTreated[n2.id()] = true;
					map_isTreated[n3.id()] = true;
					compteur_treated_bloc++;
				}

			}
		}
		std::cout << "------------" << std::endl;
		std::cout << "Compteur : " << compteur_treated_bloc << " sur " << m_map_chords[ind_chord].size() << std::endl;
	}


	std::cout << "Nbr of blocs in the chord : " << m_map_chords[ind_chord].size() << std::endl;
	std::cout << "Compteur : " << compteur_treated_bloc << std::endl;



	// Refinement
	for (auto b:m_blocking->allBlocks())
	{
		int Nx = b.getNbDiscretizationI();
		int Ny = b.getNbDiscretizationJ();

		Variable<int>* var_couche = m_blocking->getVariable<int, GMDS_NODE>("GMDS_Couche");
		Node n0 = b.getNode(0);
		Node n1 = b.getNode(1);
		Node n2 = b.getNode(2);
		Node n3 = b.getNode(3);

		if (map_refinementNeeded[n0.id()]
		    && map_refinementNeeded[n1.id()])
		{
			for (int i = 0; i < Nx; i++) {
				std::vector<math::Point> Points;
				for (int j = 0; j < Ny; j++) {
					Points.push_back(b(i, j).point());
				}
				RefinementBeta ref(Points, m_params_aero.edge_size_first_ortho_wall);
				ref.execute();
				Points = ref.GetNewPositions();

				for (int j = 1; j < Ny - 1; j++) {
					b(i, j).setPoint({Points[j]});
				}
			}
		}

		if (map_refinementNeeded[n2.id()]
		    && map_refinementNeeded[n3.id()])
		{
			for (int i = 0; i < Nx; i++) {
				std::vector<math::Point> Points;
				for (int j = Ny-1; j >= 0; j--) {
					Points.push_back(b(i, j).point());
				}
				RefinementBeta ref(Points, m_params_aero.edge_size_first_ortho_wall);
				ref.execute();
				Points = ref.GetNewPositions();

				for (int j = 1; j < Ny - 1; j++) {
					b(i, j).setPoint({Points[Ny-1-j]});
				}
			}
		}

		if (map_refinementNeeded[n0.id()]
		    && map_refinementNeeded[n3.id()])
		{
			for (int j = 0; j < Ny; j++) {
				std::vector<math::Point> Points;
				for (int i = 0; i < Nx; i++) {
					Points.push_back(b(i, j).point());
				}
				RefinementBeta ref(Points, m_params_aero.edge_size_first_ortho_wall);
				ref.execute();
				Points = ref.GetNewPositions();

				for (int i = 1; i < Nx - 1; i++) {
					b(i, j).setPoint({Points[i]});
				}
			}
		}

		if (map_refinementNeeded[n1.id()]
		    && map_refinementNeeded[n2.id()])
		{
			for (int j = 0; j < Ny; j++) {
				std::vector<math::Point> Points;
				for (int i = Nx-1; i >= 0; i--) {
					Points.push_back(b(i, j).point());
				}
				RefinementBeta ref(Points, m_params_aero.edge_size_first_ortho_wall);
				ref.execute();
				Points = ref.GetNewPositions();

				for (int i = 1; i < Nx - 1; i++) {
					b(i, j).setPoint({Points[Nx-1-i]});
				}
			}
		}

	}

}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
std::map<int, std::vector<TCellID>>
RefinementBetaBlocking::ComputeChords(){

	std::map<int, std::vector<TCellID>> map_chords;
	std::map<int, std::vector<TCellID>> map_chords_blocs;

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

	m_blocking->unmarkAll<Edge>(mark_isTreated);
	m_blocking->freeMark<Edge>(mark_isTreated);

	// Créé les cordes sous forme de blocs
	for (auto chord:map_chords)
	{
		int mark_isOnChord = m_blocking->newMark<Face>();

		std::vector<TCellID> chord_blocks;

		for (auto e_id:chord.second)
		{
			Edge e = m_blocking->get<Edge>(e_id);
			std::vector<Face> blocks = e.get<Face>();
			for (auto b:blocks)
			{
				m_blocking->mark(b, mark_isOnChord);
			}
		}

		for (auto b_id:m_blocking->faces())
		{
			Face b = m_blocking->get<Face>(b_id) ;
			if (m_blocking->isMarked(b, mark_isOnChord))
			{
				chord_blocks.push_back(b_id);
			}

		}

		map_chords_blocs[chord.first] = chord_blocks ;

		m_blocking->unmarkAll<Face>(mark_isOnChord);
		m_blocking->freeMark<Face>(mark_isOnChord);
	}



	return map_chords_blocs;

}
/*-------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
std::vector<TCellID>
RefinementBetaBlocking::ComputeOppositeEdges(TCellID e_id)
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