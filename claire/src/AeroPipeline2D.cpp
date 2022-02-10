//
// Created by rochec on 09/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AeroPipeline2D.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <unit_test_config.h>
#include <iostream>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AeroPipeline2D::AeroPipeline2D(ParamsAero Aparams) {
	m_mesh = NULL;
	m_params = Aparams;
	m_isOver = false;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline2D::execute(){

	// Je n'ai pas trouvé d'autre façon de faire pour initialiser le maillage
	// pour l'instant...
	gmds::Mesh m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
	                             gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));
	m_mesh = &m;

	LectureMaillage();
	InitialisationFronts();
	EcritureMaillage();

	m_isOver = true;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline2D::LectureMaillage(){

	std::cout << "-> Lecture du maillage ..." << std::endl;

	// Lecture du maillage

	std::string dir(TEST_SAMPLES_DIR);
	std::string vtk_file = dir+"/Carre.vtk";

	gmds::IGMeshIOService ioService(m_mesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(vtk_file);

	gmds::MeshDoctor doc(m_mesh);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline2D::EcritureMaillage(){

	std::cout << "-> Ecriture du maillage ..." << std::endl;

	gmds::IGMeshIOService ioService(m_mesh);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	std::string dir(TEST_SAMPLES_DIR);
	vtkWriter.write(m_params.output_file);

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline2D::InitialisationFronts(){

	std::cout << "-> Initialisation des fronts" << std::endl;

	// Marques sur les deux fronts d'intérêt pour l'aéro :
	// Paroi -> Noeuds sur la paroi
	// Ext -> Noeuds sur la frontière extérieur
	m_markFrontNodesParoi = m_mesh->newMark<gmds::Node>();
	m_markFrontNodesExt = m_mesh->newMark<gmds::Node>();

	//Get the boundary node ids
	BoundaryOperator2D bnd_op(m_mesh);
	std::vector<TCellID> bnd_node_ids;
	bnd_op.getBoundaryNodes(bnd_node_ids);

	int markBoundaryNodes = m_mesh->newMark<Node>();

	// Initialise une marque sur les noeuds du bord
	for (auto n_id:bnd_node_ids){
		Node n = m_mesh->get<Node>(n_id);
		m_mesh->mark(n, markBoundaryNodes);
	}

	// Variable qui contient la couleur du bord
	Variable<int>* var_color_bords ;
	var_color_bords = m_mesh->newVariable<int, GMDS_NODE>("COLOR_BORDS");

	int color = 0; //Default value is 0
	int markTreated = m_mesh->newMark<Node>();

	for (auto n_id:m_mesh->nodes())
	{
		// Si un noeud est marqué sur le bord et qu'il n'est pas encore traité
		if (m_mesh->isMarked<Node>(n_id, markBoundaryNodes) &&
		    !m_mesh->isMarked<Node>(n_id, markTreated)){
			Node n = m_mesh->get<Node>(n_id);

			// Nouveau bord, nouvelle couleur
			color++;
			m_mesh->mark(n, markTreated);
			(*var_color_bords)[n_id] = color;

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
					if(m_mesh->isMarked<Node>(n_adj_id, markBoundaryNodes) &&
					    !m_mesh->isMarked<Node>(n_adj_id, markTreated)){
						// Si le noeud est sur le bord et qu'il n'a pas été traité
						// On met à jour sa couleur et on le marque comme traité
						m_mesh->mark(n_adj, markTreated);
						(*var_color_bords)[n_adj_id] = color;

						// Ajout du noeud dans la liste des noeuds à proprager
						next.push_back(n_adj);
					}
				}

			}

		}

	}

	m_mesh->unmarkAll<Node>(markBoundaryNodes);
	m_mesh->freeMark<Node>(markBoundaryNodes);

	m_mesh->unmarkAll<Node>(markTreated);
	m_mesh->freeMark<Node>(markTreated);




	/*
	int markFrontNodes = m.newMark<gmds::Node>();
	for(auto id:bnd_node_ids){
		Node n = m.get<Node>(id);
		double coord_y = n.Y() ;
		if (coord_y == 0) {
			// For the square test case, the front to advance is the boundary where y=0
			m.mark<Node>(id,markFrontNodes);
		}
	}
	 */
}
/*------------------------------------------------------------------------*/