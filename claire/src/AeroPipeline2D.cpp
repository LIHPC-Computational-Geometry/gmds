//
// Created by rochec on 09/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AeroPipeline2D.h>
#include <gmds/claire/LevelSetCombined.h>
#include <gmds/claire/LeastSquaresGradientComputation.h>
#include <gmds/claire/AdvectedPointRK4_2D.h>
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

AeroPipeline2D::AeroPipeline2D(ParamsAero Aparams) :
  AbstractAeroPipeline(Aparams),
  m_m(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
                                       gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F)),
  m_mGen(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
                                       gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F))
{
	m_couche_id = m_mGen.newVariable<int, GMDS_NODE>("GMDS_Couche_Id");
	m_mesh = &m_m;
	m_meshGen = &m_mGen;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline2D::execute(){

	LectureMaillage();
	InitialisationFronts();

	// Calcul du level set
	m_mesh->newVariable<double,GMDS_NODE>("GMDS_Distance");
	m_mesh->newVariable<double,GMDS_NODE>("GMDS_Distance_Int");
	m_mesh->newVariable<double,GMDS_NODE>("GMDS_Distance_Out");
	LevelSetCombined lsCombined(m_mesh, m_markFrontNodesParoi, m_markFrontNodesExt,
	                            m_mesh->getVariable<double,GMDS_NODE>("GMDS_Distance"),
	                            m_mesh->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
	                            m_mesh->getVariable<double,GMDS_NODE>("GMDS_Distance_Out"));
	lsCombined.execute();

	// Calcul du gradient du champ de Level Set
	m_mesh->newVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient");
	LeastSquaresGradientComputation grad2D(m_mesh, m_mesh->getVariable<double,GMDS_NODE>("GMDS_Distance"),
	                                       m_mesh->getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	grad2D.execute();

	InitialisationMeshGen();

	GenerationCouche(1, 0.05);
	GenerationCouche(2, 1);

	EcritureMaillage(m_meshGen);

	m_isOver = true;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline2D::LectureMaillage(){

	// Lecture du maillage
	std::cout << "-> Lecture du maillage ..." << std::endl;

	gmds::IGMeshIOService ioService(m_mesh);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(m_params.input_file);

	gmds::MeshDoctor doc(m_mesh);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	MeshCleaner();

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline2D::MeshCleaner(){
	for (auto n_id:m_mesh->nodes())
	{
		Node n = m_mesh->get<Node>(n_id);
		if (n.get<Face>().empty()) {
			//std::cout << "Noeud isolé : " << n_id << std::endl;
			m_mesh->deleteNode(n_id);
		}
	}
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline2D::EcritureMaillage(Mesh* p_mesh){

	std::cout << "-> Ecriture du maillage ..." << std::endl;

	gmds::IGMeshIOService ioService(p_mesh);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	std::string dir(TEST_SAMPLES_DIR);
	vtkWriter.write(m_params.output_file);

	// Ecriture du maillage en triangles initial pour visualisation et débug
	ioService = m_mesh;
	gmds::VTKWriter vtkWriter2(&ioService);
	vtkWriter2.setCellOptions(gmds::N|gmds::F);
	vtkWriter2.setDataOptions(gmds::N|gmds::F);
	vtkWriter2.write("AeroPipeline2D_Test1_Triangles.vtk");

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline2D::InitialisationFronts(){

	std::cout << "-> Initialisation des fronts" << std::endl;

	// Marques sur les deux fronts d'intérêt pour l'aéro :
	// Paroi -> Noeuds sur la paroi
	// Ext -> Noeuds sur la frontière extérieure
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

	if (color < 2){
		std::cout << "Attention : il n'y a pas deux bords minimum." << std::endl;
	}

	m_mesh->unmarkAll<Node>(markBoundaryNodes);
	m_mesh->freeMark<Node>(markBoundaryNodes);

	m_mesh->unmarkAll<Node>(markTreated);
	m_mesh->freeMark<Node>(markTreated);


	// Calcul des boites englobantes
	std::vector<double> x_min(color,0.0);
	std::vector<double> y_min(color,0.0);
	std::vector<double> x_max(color,0.0);
	std::vector<double> y_max(color,0.0);

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
				    (x_max[j-1] < x_max[i-1]) &&
				    (y_max[j-1] < y_max[i-1]) ) {
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

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline2D::InitialisationMeshGen(){

	std::cout << "-> Discrétisation du front int du maillage quad" << std::endl;

	//Get the boundary node ids
	BoundaryOperator2D bnd_op(m_mesh);
	std::vector<TCellID> bnd_node_ids;
	bnd_op.getBoundaryNodes(bnd_node_ids);

	Variable<double> *var_newID = m_mesh->newVariable<double, GMDS_NODE>("GMDS_NewID");
	for (auto n_id:m_mesh->nodes()){
		var_newID->set(n_id, NullID);
	}

	double err=pow(10,-6);

	// Ajout des noeuds de la frontière intérieure au maillage à généré
	for (auto n_id:bnd_node_ids){
		Node n = m_mesh->get<Node>(n_id);
		math::Point p = n.point();
		//std::cout << "node : " << n << std::endl;
		if ( m_mesh->isMarked<Node>(n_id, m_markFrontNodesParoi) ){
			Node n_new = m_meshGen->newNode(p);
			m_couche_id->set(n_new.id(), 0);
			var_newID->set(n_id, n_new.id());	// On stocke l'indice du noeud n dans le nouveau maillage sur l'ancien maillage
		}
	}

	// Ajout des arêtes
	for (auto e_id:m_mesh->edges()){
		Edge e = m_mesh->get<Edge>(e_id);
		std::vector<Node> adj_nodes = e.get<Node>();
		TCellID n0_id = adj_nodes[0].id();
		TCellID n1_id = adj_nodes[1].id();
		if ( m_mesh->isMarked<Node>(n0_id, m_markFrontNodesParoi) &&
		   m_mesh->isMarked<Node>(n1_id, m_markFrontNodesParoi)){
			TCellID n0_id_new = var_newID->value(n0_id) ;
			TCellID n1_id_new = var_newID->value(n1_id) ;
			Edge e_new = m_meshGen->newEdge(n0_id_new, n1_id_new);
			// Ajout des connectivités Node -> Edge
			Node n0 = m_meshGen->get<Node>(n0_id_new);
			Node n1 = m_meshGen->get<Node>(n1_id_new);
			n0.add<Edge>(e_new);
			n1.add<Edge>(e_new);
		}
	}

	gmds::MeshDoctor doc(m_meshGen);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	m_mesh->deleteVariable(GMDS_NODE, "GMDS_NewID");

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline2D::GenerationCouche(int couche_id, double dist){

	std::cout << "-> Génération de la couche " << couche_id << std::endl ;

	std::vector<Node> nodes_couche_id;	// Noeuds sur la douche couche_id

	// Création des nouveaux noeuds de la couche couche_id naïvement
	// + aretes entre la couche couche_id-1 et couche_id
	for (auto n_id:m_meshGen->nodes()){
		if (m_couche_id->value(n_id) == couche_id-1){
			// Placement du point P à la distance souhaitée suivant le champ de gradient
			Node n = m_meshGen->get<Node>(n_id);
			math::Point M = n.point();
			AdvectedPointRK4_2D advpoint(m_mesh, M, dist, m_m.getVariable<double, GMDS_NODE>("GMDS_Distance"),
			   m_m.getVariable<math::Vector3d ,GMDS_NODE>("GMDS_Gradient"));
			advpoint.execute();
			math::Point P = advpoint.getPend();

			// Nouveau noeud dans le maillage
			Node n_new = m_mGen.newNode(P);
			m_couche_id->set(n_new.id(),couche_id);
			nodes_couche_id.push_back(n_new);

			// Création de l'arête pour relier le noeud à la couche n-1
			Edge e_new = m_mGen.newEdge(n, n_new);
			n_new.add<Edge>(e_new);	// Ajout de la connectivité N -> E
			n.add<Edge>(e_new);	// Ajout de la connectivité N -> E

		}
	}

	// Création des faces de la couche
	for (auto n0:nodes_couche_id){

		Node n1 = AnteriorNode(n0);
		std::vector<Node> adj_nodes_n1_in_anterior_layer = AdjNodesInLayer(n1, couche_id-1);
		// Première face, on regarde si elle existe, si elle n'existe pas, on la créé
		bool exist = isQuadCreated(n0, n1, adj_nodes_n1_in_anterior_layer[0]);
		if (!exist){
			CreateQuadAndConnectivities(n0, n1, adj_nodes_n1_in_anterior_layer[0]);
		}
		// Deuxième face, on regarde si elle existe, si elle n'existe pas, on la créé
		exist = isQuadCreated(n0, n1, adj_nodes_n1_in_anterior_layer[1]);
		if (!exist){
			CreateQuadAndConnectivities(n0, n1, adj_nodes_n1_in_anterior_layer[1]);
		}

	}

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
Node AeroPipeline2D::SuccessorNode(Node n0){
	Node n;
	int couche_i = m_couche_id->value(n0.id());
	std::vector<Edge> adj_edges = n0.get<Edge>() ;
	for (auto e:adj_edges){
		Node n_opp = e.getOppositeNode(n0);
		if (m_couche_id->value(n_opp.id()) == couche_i+1){
			n = n_opp;
		}
	}
	return n;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
Node AeroPipeline2D::AnteriorNode(Node n0){
	Node n;
	int couche_i = m_couche_id->value(n0.id());
	std::vector<Edge> adj_edges = n0.get<Edge>() ;
	for (auto e:adj_edges){
		Node n_opp = e.getOppositeNode(n0);
		if (m_couche_id->value(n_opp.id()) == couche_i-1){
			n = n_opp;
		}
	}
	return n;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
std::vector<Node> AeroPipeline2D::AdjNodesInLayer(Node n0, int couche_i){
	std::vector<Node> adj_nodes_in_layer ;
	std::vector<Edge> adj_edges = n0.get<Edge>() ; // Ensemble des arêtes adjacentes au noeud n0
	for (auto e:adj_edges){
		Node n_opp = e.getOppositeNode(n0);
		if (m_couche_id->value(n_opp.id()) == couche_i){
			adj_nodes_in_layer.push_back(n_opp);
		}
	}
	return adj_nodes_in_layer;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline2D::CreateQuadAndConnectivities(Node n0, Node n1, Node n2){

	Node n3 = SuccessorNode(n2);	// On récupère le noeud n3, 4ème noeud du quad, dans la couche i

	// Création de l'arête pour relier les doeux noeuds de la couche i
	Edge e3 = m_mGen.newEdge(n0, n3);
	n0.add<Edge>(e3);	// Connectivités N->E
	n3.add<Edge>(e3);

	// Création de la face de type quad
	Face f = m_mGen.newQuad(n0, n1, n2, n3) ;
	n0.add<Face>(f);	// Connectivités N -> F
	n1.add<Face>(f);
	n2.add<Face>(f);
	n3.add<Face>(f);
	f.add<Edge>(e3); // Connectivité F -> E
	e3.add<Face>(f); // Connectivité E -> F

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
bool AeroPipeline2D::isQuadCreated(Node n0, Node n1, Node n2){
	bool exist(false);

	// On regarde si l'arête entre les noeuds n0 et n3 est construite
	Node n3 = SuccessorNode(n2) ;
	TCellID n3_id = n3.id();
	std::vector<Edge> adj_edges = n0.get<Edge>() ;

	for (auto e:adj_edges){
		Node n_opp = e.getOppositeNode(n0);
		if (n_opp.id() == n3_id){
			exist = true;
		}
	}

	return exist;
}
/*------------------------------------------------------------------------*/