//
// Created by rochec on 09/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AeroPipeline3D.h>
#include <gmds/claire/LevelSetCombined.h>
#include <gmds/claire/LeastSquaresGradientComputation.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <unit_test_config.h>
#include <iostream>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AeroPipeline3D::AeroPipeline3D(ParamsAero Aparams) :
	AbstractAeroPipeline(Aparams),
  m_mTetra(gmds::MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
                                        F2E | E2F | R2E | N2R | N2F | N2E)),
  m_mHexa( gmds::MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
                         F2E | E2F | R2E | N2R | N2F | N2E) )
{
	m_couche_id = m_mHexa.newVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	m_mesh = &m_mTetra;
	m_meshGen = &m_mHexa;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline3D::execute(){
	LectureMaillage();
	EcritureMaillage();
	m_isOver = true;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline3D::LectureMaillage(){

	// Lecture du maillage
	std::cout << "-> Lecture du maillage ..." << std::endl;

	gmds::IGMeshIOService ioService(m_mesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::R);
	vtkReader.read(m_params.input_file);

	gmds::MeshDoctor doctor(m_mesh);
	doctor.buildFacesAndR2F();
	doctor.buildEdgesAndX2E();
	doctor.updateUpwardConnectivity();

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline3D::EcritureMaillage(){

	std::cout << "-> Ecriture du maillage ..." << std::endl;

	// Ecriture du maillage généré
	gmds::IGMeshIOService ioService(m_meshGen);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	std::string dir(TEST_SAMPLES_DIR);
	vtkWriter.write(m_params.output_file);

	// Ecriture du maillage initial (tetra)
	ioService = m_mesh;
	gmds::VTKWriter vtkWriter2(&ioService);
	vtkWriter2.setCellOptions(gmds::N|gmds::F);
	vtkWriter2.setDataOptions(gmds::N|gmds::F);
	vtkWriter2.write("AeroPipeline3D_Test1_Tetra.vtk");

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline3D::InitialisationFronts(){

	std::cout << "-> Initialisation des fronts" << std::endl;

	// Marques sur les deux fronts d'intérêt pour l'aéro :
	// Paroi -> Noeuds sur la paroi
	// Ext -> Noeuds sur la frontière extérieure
	m_markFrontNodesParoi = m_mesh->newMark<gmds::Node>();
	m_markFrontNodesExt = m_mesh->newMark<gmds::Node>();

	//Get the boundary node ids
	BoundaryOperator bnd_op(m_mesh);
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
				std::vector<Node> adjacent_nodes = AdjacentNodes(current_node) ;

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

}
/*------------------------------------------------------------------------*/