//
// Created by rochec on 09/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AeroPipeline3D.h>
#include <gmds/claire/LevelSetCombined.h>
#include <gmds/claire/LeastSquaresGradientComputation.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
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
AbstractAeroPipeline::STATUS AeroPipeline3D::execute(){
	LectureMaillage();
	InitialisationFronts();

	// Calcul du level set
	std::cout << "-> Calcul des Level Sets" << std::endl;
	m_mesh->newVariable<double,GMDS_NODE>("GMDS_Distance");
	m_mesh->newVariable<double,GMDS_NODE>("GMDS_Distance_Int");
	m_mesh->newVariable<double,GMDS_NODE>("GMDS_Distance_Out");
	LevelSetCombined lsCombined(m_mesh, m_markFrontNodesParoi, m_markFrontNodesExt,
	                            m_mesh->getVariable<double,GMDS_NODE>("GMDS_Distance"),
	                            m_mesh->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
	                            m_mesh->getVariable<double,GMDS_NODE>("GMDS_Distance_Out"));
	lsCombined.execute();

	// Calcul du gradient du champ de Level Set
	std::cout << "-> Calcul du gradient du champ des Level Sets" << std::endl;
	m_mesh->newVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient");
	LeastSquaresGradientComputation grad2D(m_mesh, m_mesh->getVariable<double,GMDS_NODE>("GMDS_Distance"),
	                                       m_mesh->getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	grad2D.execute();

	EcritureMaillage();

	return AbstractAeroPipeline::SUCCESS;
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
	std::vector<TCellID> bnd_node_ids;
	bnd_node_ids = getBndNodes();

	int markBoundaryNodes = m_mesh->newMark<Node>();

	// Initialise une marque sur les noeuds du bord
	for (auto n_id:bnd_node_ids){
		Node n = m_mesh->get<Node>(n_id);
		m_mesh->mark(n, markBoundaryNodes);
	}

	// Variable qui contient la couleur du bord
	Variable<int>* var_color_bords ;
	var_color_bords = m_mesh->newVariable<int, GMDS_NODE>("COLOR_BORDS");

	// ----- Marque chaque bord avec une couleur différente -----
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

	// ----- Fin de marquage des différents bords -----

	if (color < 2){
		std::cout << "Attention : il n'y a pas deux bords minimum." << std::endl;
	}

	m_mesh->unmarkAll<Node>(markBoundaryNodes);
	m_mesh->freeMark<Node>(markBoundaryNodes);

	m_mesh->unmarkAll<Node>(markTreated);
	m_mesh->freeMark<Node>(markTreated);


	// ----- Initialisation quel bord est le front ext, les autres bords sont paroi -----
	// Calcul des boites englobantes
	std::vector<double> x_min(color,0.0);
	std::vector<double> y_min(color,0.0);
	std::vector<double> z_min(color,0.0);
	std::vector<double> x_max(color,0.0);
	std::vector<double> y_max(color,0.0);
	std::vector<double> z_max(color,0.0);

	for (auto n_id:bnd_node_ids ) {
		Node n = m_mesh->get<Node>(n_id);
		math::Point p = n.point();
		// A quel bord appartient le noeud n_id
		int couleur = var_color_bords->value(n_id);
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
		if(p.Z() < z_min[couleur-1]){
			z_min[couleur-1] = p.Z();
		}
		if(p.Z() > z_max[couleur-1]){
			z_max[couleur-1] = p.Z();
		}
	}

	// On recherche la plus grosse boite englobante
	int bigest_color = 1;
	for(int i=2;i<=color;i++){
		// On compte le nombre de boites que la couleur i englobe
		int nbr_boite_englobe = 0;
		for (int j=1;j<=color;j++){
			if (i != j){
				if ( (x_min[i-1] < x_min[j-1]) &&
				    (y_min[i-1] < y_min[j-1]) &&
				    (z_min[i-1] < z_min[j-1]) &&
				    (x_max[j-1] < x_max[i-1]) &&
				    (y_max[j-1] < y_max[i-1]) &&
				    (z_max[j-1] < z_max[i-1]) ) {
					nbr_boite_englobe += 1;
				}
			}
		}
		if(nbr_boite_englobe == color-1){
			bigest_color = i;
		}
	}
	//std::cout << "Boite englobante : " << bigest_color << std::endl;

	// On marque les fronts paroi et extérieur
	// Variable qui contient la couleur du bord
	Variable<int>* var_color_paroi ;
	var_color_paroi = m_mesh->newVariable<int, GMDS_NODE>("COLOR_PAROI");
	for (auto n_id:bnd_node_ids){
		int couleur = var_color_bords->value(n_id);
		if(couleur == bigest_color){
			m_mesh->mark<Node>(n_id,m_markFrontNodesExt);
			(*var_color_paroi)[n_id] = 1;
		}
		else{
			m_mesh->mark<Node>(n_id,m_markFrontNodesParoi);
			(*var_color_paroi)[n_id] = 2;
		}
	}

	// ----- Fin initialisation des bords -----

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
std::vector<TCellID> AeroPipeline3D::getBndNodes(){
	std::vector<TCellID> bnd_nodes;

	for (auto f_id:m_mesh->faces())
	{
		Face f = m_mesh->get<Face>(f_id);

		if (f.get<Region>().size() == 1)
		{
			std::vector<Node> face_nodes = f.get<Node>() ;
			for (auto n:face_nodes){
				bnd_nodes.push_back(n.id());
			}
		}
	}

	return bnd_nodes;
}
/*------------------------------------------------------------------------*/
