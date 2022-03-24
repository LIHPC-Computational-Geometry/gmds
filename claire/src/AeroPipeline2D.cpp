//
// Created by rochec on 09/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/Utils.h>
#include <gmds/claire/AeroPipeline2D.h>
#include <gmds/claire/AeroBoundaries_2D.h>
#include <gmds/claire/LevelSetCombined.h>
#include <gmds/claire/LeastSquaresGradientComputation.h>
#include <gmds/claire/AdvectedPointRK4_2D.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
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
	//m_var_color_bords = m_m.newVariable<int, GMDS_NODE>("COLOR_BORDS");
	m_mesh = &m_m;
	m_meshGen = &m_mGen;
	m_Bnd = new AeroBoundaries_2D(m_mesh) ;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
AbstractAeroPipeline::STATUS AeroPipeline2D::execute(){

	LectureMaillage();
	m_Bnd->execute();
	//InitialisationFronts();

	// Calcul du level set
	m_mesh->newVariable<double,GMDS_NODE>("GMDS_Distance");
	m_mesh->newVariable<double,GMDS_NODE>("GMDS_Distance_Int");
	m_mesh->newVariable<double,GMDS_NODE>("GMDS_Distance_Out");
	LevelSetCombined lsCombined(m_mesh, m_Bnd->getMarkParoi(), m_Bnd->getMarkAmont(),
	                            m_mesh->getVariable<double,GMDS_NODE>("GMDS_Distance"),
	                            m_mesh->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
	                            m_mesh->getVariable<double,GMDS_NODE>("GMDS_Distance_Out"));
	lsCombined.execute();

	// Calcul du gradient du champ de Level Set
	m_mesh->newVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient");
	LeastSquaresGradientComputation grad2D(m_mesh, m_mesh->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
	                                       m_mesh->getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	grad2D.execute();

	//InitialisationMeshGen();
	//DiscretisationBlocParoi();

	for (int i=1;i<=m_Bnd->getNbrBords();i++){
		if(i != m_Bnd->getColorAmont()){
			DiscretisationBlocsBord(i);
		}
	}

	GenerationCouche(1, 0.25);
	GenerationCouche(2, 1);


	EcritureMaillage();

	return AbstractAeroPipeline::SUCCESS;

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

	math::Utils::MeshCleaner(m_mesh);

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline2D::EcritureMaillage(){

	std::cout << "-> Ecriture du maillage ..." << std::endl;

	gmds::IGMeshIOService ioService(m_meshGen);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	std::string dir(".");//TEST_SAMPLES_DIR);
	vtkWriter.write(m_params.output_file);

	// Ecriture du maillage en triangles initial pour visualisation et débug
	ioService = m_mesh;
	gmds::VTKWriter vtkWriter2(&ioService);
	vtkWriter2.setCellOptions(gmds::N|gmds::F);
	vtkWriter2.setDataOptions(gmds::N|gmds::F);
	vtkWriter2.write("AeroPipeline2D_Triangles.vtk");

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
/*
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
	//Variable<int>* var_color_bords ;
	//var_color_bords = m_mesh->newVariable<int, GMDS_NODE>("COLOR_BORDS");

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
			//(*m_var_color_bords)[n_id] = color;

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
						//(*m_var_color_bords)[n_adj_id] = color;

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
	else{
		m_nbrBordsParoi = color-1 ;
	}

	m_nbrBords = color;

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
		//int couleur = m_var_color_bords->value(n_id);
		int couleur = m_Bnd->getNodeColor(n_id);
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
	m_bigest_color = 1;
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
			m_bigest_color = i;
		}
	}
	//std::cout << "Boite englobante : " << bigest_color << std::endl;

	// On marque les fronts paroi et extérieur
	// Variable qui contient la couleur du bord
	Variable<int>* var_color_paroi ;
	var_color_paroi = m_mesh->newVariable<int, GMDS_NODE>("COLOR_PAROI");
	for (auto n_id:bnd_node_ids){
		//int couleur = m_var_color_bords->value(n_id);
		int couleur = m_Bnd->getNodeColor(n_id);
		if(couleur == m_bigest_color){
			m_mesh->mark<Node>(n_id,m_markFrontNodesExt);
			(*var_color_paroi)[n_id] = 1;
		}
		else{
			m_mesh->mark<Node>(n_id,m_markFrontNodesParoi);
			(*var_color_paroi)[n_id] = 2;
		}
	}


}
*/
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline2D::DiscretisationBlocParoi(){

	std::cout << "-> Discrétisation du front int du maillage quad" << std::endl;

	//Get the boundary node ids
	BoundaryOperator2D bnd_op(m_mesh);
	std::vector<TCellID> bnd_node_ids;
	bnd_op.getBoundaryNodes(bnd_node_ids);

	int markCurveEdges = m_mesh->newMark<Edge>();
	int markCurveNodes = m_mesh->newMark<Node>();
	int markPointNodes = m_mesh->newMark<Node>();
	int markAloneNodes = m_mesh->newMark<Node>();

	bnd_op.markCellOnGeometry(markCurveEdges, markCurveNodes,
	                 markPointNodes, markAloneNodes);

	m_mesh->unmarkAll<Node>(markAloneNodes);
	m_mesh->freeMark<Node>(markAloneNodes);



	// Discrétisation pour un bord de paroi
	int color = 1;

	// On récuère un noeud de ce bord qui est sur un sommet
	TCellID n0_id = NullID ;
	for(auto n_it = m_mesh->nodes_begin(); n_it!= m_mesh->nodes_end() && (n0_id == NullID);++n_it){
		TCellID n_id = *n_it;
		if ( m_mesh->isMarked<Node>(n_id, m_markFrontNodesParoi) && m_mesh->isMarked<Node>(n_id, markPointNodes) ){
			n0_id = n_id;
		}
	}

	// Si le maillage ne comporte pas de sommet, alors, on prend un noeud au hasard
	double x_min = std::numeric_limits<double>::max();
	/*
	if (n0_id==NullID) {
		for (auto n_it = m_mesh->nodes_begin(); n_it != m_mesh->nodes_end() && (n0_id == NullID); ++n_it) {
			TCellID n_id = *n_it;
			Node n = m_mesh->get<Node>(n_id);
			math::Point p = n.point();
			if (m_mesh->isMarked<Node>(n_id, m_markFrontNodesParoi) && p.X() < x_min ) {
				n0_id = n_id;
				x_min = p.X();
			}
		}
	}
	*/

	n0_id = PointArret(color);

	// Initialisation des longueurs pour définir les sommets de blocs
	Node n0 = m_mesh->get<Node>(n0_id);
	x_min = n0.X();
	double perim = ComputeBoundaryLength(color);		// Calcul du périmètre du bord regardé
	double Lmax = perim/ m_params.nbrMinBloc ;
	double l(0);
	std::cout << "Périmètre : " << perim << std::endl ;
	//std::cout << "Limite taille bloc : " << Lmax << std::endl;


	// Parcours du bord de proche en proche en passant par l'opposé
	TCellID n1_id = n0_id;
	TCellID n2_id = NullID;
	std::vector<Edge> adj_edges = n0.get<Edge>();
	for (auto e:adj_edges){
		Node ne = e.getOppositeNode(n0);
		if ( m_mesh->isMarked<Node>(ne.id(), m_markFrontNodesParoi) ){
			n2_id = ne.id();
			l = e.length();
		}
	}

	Node nstart_quad = m_meshGen->newNode(n0.point()); // Premier noeud du nouveau maillage
	Node n0_quad = nstart_quad;
	Node n1_quad;
	TCellID nstart_quad_id = nstart_quad.id() ;
	TCellID n0_quad_id = nstart_quad_id;
	TCellID n1_quad_id = NullID;

	//std::cout << "First node : " << n0_quad.point() << std::endl ;

	while (n2_id != n0_id){

		Node n_test = m_mesh->get<Node>(n2_id);
		math::Point p2 = n_test.point();

		// Si le noeud est sur un sommet
		if ( m_mesh->isMarked<Node>(n2_id, markPointNodes) || l >= Lmax || abs(l-Lmax) <= pow(10,-6) || abs(p2.X() - x_min) <= pow(10,-6) ){
			Node n2 = m_mesh->get<Node>(n2_id);
			n1_quad = m_meshGen->newNode(n2.point());
			Edge e = m_meshGen->newEdge(n0_quad, n1_quad);
			// Ajout des connectivités Node -> Edge
			n0_quad.add<Edge>(e);
			n1_quad.add<Edge>(e);
			// Mise à jour du noeud n0
			n0_quad = n1_quad;
			n0_quad_id = n1_quad.id();
			//std::cout << "l = " << l << std::endl;
			//std::cout << "pos : " << n2.point() << std::endl;
			//std::cout << "pos quad : " << n1_quad.point() << std::endl;
			l = 0;	// On remet l à 0
		}

		// On cherche le prochain noeud
		TCellID n_old_id = n1_id;
		n1_id = n2_id;
		Node n2 = m_mesh->get<Node>(n2_id);
		adj_edges = n2.get<Edge>();
		for (auto e:adj_edges){
			Node ne = e.getOppositeNode(n2);
			if ( m_mesh->isMarked<Node>(ne.id(), m_markFrontNodesParoi) && ne.id() != n_old_id ){
				n2_id = ne.id();
				l = l + e.length();
			}
		}
		//std::cout << "l = " << l << std::endl;
	}

	// Ajout de la dernière arête pour fermer le blocking du bord
	Edge e = m_meshGen->newEdge(n1_quad, nstart_quad);
	// Ajout des connectivités Node -> Edge
	n1_quad.add<Edge>(e);
	nstart_quad.add<Edge>(e);

	m_mesh->unmarkAll<Node>(markCurveEdges);
	m_mesh->freeMark<Node>(markCurveEdges);

	m_mesh->unmarkAll<Node>(markCurveNodes);
	m_mesh->freeMark<Node>(markCurveNodes);

	m_mesh->unmarkAll<Node>(markPointNodes);
	m_mesh->freeMark<Node>(markPointNodes);

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


/*------------------------------------------------------------------------*/
double AeroPipeline2D::ComputeBoundaryLength(int color){

	double length(0.0);
	std::vector<TCellID> bnd_nodes_id_ordered = BndNodesOrdered(color);

	for(int i=0;i<bnd_nodes_id_ordered.size()-1;i++){
		Node n0 = m_mesh->get<Node>(bnd_nodes_id_ordered[i]);
		Node n1 = m_mesh->get<Node>(bnd_nodes_id_ordered[i+1]);
		math::Point p0 = n0.point();
		math::Point p1 = n1.point();
		math::Vector3d v = p1-p0;
		length += v.norm();
	}

	// Ajout de la dernière longueur pour "fermer la boucle"
	Node n0 = m_mesh->get<Node>(bnd_nodes_id_ordered[bnd_nodes_id_ordered.size()-1]);
	Node n1 = m_mesh->get<Node>(bnd_nodes_id_ordered[0]);
	math::Point p0 = n0.point();
	math::Point p1 = n1.point();
	math::Vector3d v = p1-p0;
	length += v.norm();

	return length;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
TCellID AeroPipeline2D::PointArret(int color){

	TCellID n_arret_id(NullID);
	double x_min = std::numeric_limits<double>::max();

	for (auto n_it = m_mesh->nodes_begin(); n_it != m_mesh->nodes_end() ; ++n_it) {
		TCellID n_id = *n_it;
		Node n = m_mesh->get<Node>(n_id);
		math::Point p = n.point();
		//if (m_mesh->isMarked<Node>(n_id, m_markFrontNodesParoi) &&
		//    m_var_color_bords->value(n_id) == color &&
		//    p.X() < x_min ) {
		if (m_mesh->isMarked<Node>(n_id, m_markFrontNodesParoi) &&
		    m_Bnd->getNodeColor(n_id) == color &&
		    p.X() < x_min ) {
			n_arret_id = n_id;
			x_min = p.X();
		}
	}
	return n_arret_id;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
std::vector<TCellID> AeroPipeline2D::BndNodesOrdered(int color){

	std::vector<TCellID> bnd_nodes_id_ordered;

	//Get the boundary node ids
	BoundaryOperator2D bnd_op(m_mesh);
	std::vector<TCellID> bnd_node_ids;
	bnd_op.getBoundaryNodes(bnd_node_ids);

	int markCurveEdges = m_mesh->newMark<Edge>();
	int markCurveNodes = m_mesh->newMark<Node>();
	int markPointNodes = m_mesh->newMark<Node>();
	int markAloneNodes = m_mesh->newMark<Node>();

	bnd_op.markCellOnGeometry(markCurveEdges, markCurveNodes,
	                          markPointNodes, markAloneNodes);

	m_mesh->unmarkAll<Node>(markAloneNodes);
	m_mesh->freeMark<Node>(markAloneNodes);

	TCellID n0_id = NullID ; // Choix du noeud de départ

	// STRATEGIE 1 : On récuère un noeud de ce bord qui est sur un sommet
	for(auto n_it = m_mesh->nodes_begin(); n_it!= m_mesh->nodes_end() && (n0_id == NullID);++n_it){
		TCellID n_id = *n_it;
		//std::cout << "Test 1 : " << n_id << " " << m_mesh->isMarked<Node>(n_id, m_markFrontNodesParoi) << std::endl;
		if ( m_mesh->isMarked<Node>(n_id, m_markFrontNodesParoi) &&
		    m_Bnd->getNodeColor(n_id) == color &&
		   m_mesh->isMarked<Node>(n_id, markPointNodes) ){
			n0_id = n_id;
		}
	}

	m_mesh->unmarkAll<Edge>(markCurveEdges);
	m_mesh->freeMark<Edge>(markCurveEdges);

	m_mesh->unmarkAll<Node>(markCurveNodes);
	m_mesh->freeMark<Node>(markCurveNodes);

	m_mesh->unmarkAll<Node>(markPointNodes);
	m_mesh->freeMark<Node>(markPointNodes);

	// STRATEGIE 2 : On prend le point d'arrêt, positionné au x_min.
	// Attention, plusieurs noeuds peuvent être au x_min. Cette méthode
	// en retourne un au hasard.
	n0_id = PointArret(color);

	// STRATEGIE 3 : On sélectionne un noeud au hasard
	/*
	for (auto n_it = m_mesh->nodes_begin(); n_it != m_mesh->nodes_end() && (n0_id == NullID); ++n_it) {
	   TCellID n_id = *n_it;
	   if (m_mesh->isMarked<Node>(n_id, m_markFrontNodesParoi) &&
	      m_var_color_bords->value(n_id) == color) {
	      n0_id = n_id;
	   }
	}
	 */

	// Initialisation de la première valeur du vecteur
	bnd_nodes_id_ordered.push_back(n0_id);

	// Initialisation des valeurs pour le parcours du bord
	TCellID n1_id = n0_id;
	TCellID n2_id = n0_id;
	TCellID n3_id = n0_id;

	Node n0 = m_mesh->get<Node>(n0_id);
	std::vector<Edge> adj_edges = n0.get<Edge>();
	for (auto e:adj_edges){
		Node ne = e.getOppositeNode(n0);
		if ( m_mesh->isMarked<Node>(ne.id(), m_markFrontNodesParoi) ){
			n3_id = ne.id();
		}
	}

	while (n3_id != n0_id){

		bnd_nodes_id_ordered.push_back(n3_id);		// On ajout l'id du noeud n3 au vecteur

		// On cherche le prochain noeud
		n1_id = n2_id;
		n2_id = n3_id;
		Node n2 = m_mesh->get<Node>(n2_id);
		adj_edges = n2.get<Edge>();
		for (auto e:adj_edges){
			Node ne = e.getOppositeNode(n2);
			if ( m_mesh->isMarked<Node>(ne.id(), m_markFrontNodesParoi) &&
			   ne.id() != n1_id ){
				n3_id = ne.id();
			}
		}

	}

	return bnd_nodes_id_ordered;

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline2D::DiscretisationBlocsBord(int color){

	std::cout << "-> Discrétisation du bord de couleur " << color << std::endl;

	std::vector<TCellID> bnd_nodes_id_ordered = BndNodesOrdered(color);

	//Get the boundary node ids
	BoundaryOperator2D bnd_op(m_mesh);
	std::vector<TCellID> bnd_node_ids;
	bnd_op.getBoundaryNodes(bnd_node_ids);

	int markCurveEdges = m_mesh->newMark<Edge>();
	int markCurveNodes = m_mesh->newMark<Node>();
	int markPointNodes = m_mesh->newMark<Node>();
	int markAloneNodes = m_mesh->newMark<Node>();

	bnd_op.markCellOnGeometry(markCurveEdges, markCurveNodes,
	                          markPointNodes, markAloneNodes);

	m_mesh->unmarkAll<Node>(markAloneNodes);
	m_mesh->freeMark<Node>(markAloneNodes);

	// Calcul de la longueur max des arêtes de bloc
	double perim = ComputeBoundaryLength(color);		// Calcul du périmètre du bord regardé
	double Lmax = perim/ m_params.nbrMinBloc ;
	double l(0);

	// Calcul de la valeur min de x sur les noeuds du bord
	TCellID n_arret_tri_id = PointArret(color);
	Node n_arret_tri = m_mesh->get<Node>(n_arret_tri_id) ;
	math::Point point_arret = n_arret_tri.point();
	double x_min = point_arret.X() ;

	Node n0_tri = m_mesh->get<Node>(bnd_nodes_id_ordered[0]) ;
	Node n0_quad = m_meshGen->newNode(n0_tri.point()); // Premier noeud du nouveau maillage

	Node n1_quad = n0_quad;
	Node n2_quad = n0_quad;

	for(int i=1;i<bnd_nodes_id_ordered.size();i++){

		Node n = m_mesh->get<Node>(bnd_nodes_id_ordered[i]);
		math::Point p = n.point() ;
		l += math::Utils::distFromNodeIds(m_mesh, bnd_nodes_id_ordered[i], bnd_nodes_id_ordered[i-1]);

		if ( m_mesh->isMarked<Node>(bnd_nodes_id_ordered[i], markPointNodes)
		    || l >= Lmax
		    || abs(l-Lmax) <= pow(10,-6)
		    || abs(p.X() - x_min) <= pow(10,-6) ){

			n1_quad = n2_quad;
			n2_quad = m_meshGen->newNode(n.point());
			Edge e = m_meshGen->newEdge(n1_quad, n2_quad);

			// Ajout des connectivités Node -> Edge
			n1_quad.add<Edge>(e);
			n2_quad.add<Edge>(e);

			l = 0;	// On remet l à 0

		}

	}

	// Ajout de la dernière arête pour fermer le blocking du bord
	Edge e = m_meshGen->newEdge(n2_quad, n0_quad);
	// Ajout des connectivités Node -> Edge
	n2_quad.add<Edge>(e);
	n0_quad.add<Edge>(e);

	m_mesh->unmarkAll<Edge>(markCurveEdges);
	m_mesh->freeMark<Edge>(markCurveEdges);

	m_mesh->unmarkAll<Node>(markCurveNodes);
	m_mesh->freeMark<Node>(markCurveNodes);

	m_mesh->unmarkAll<Node>(markPointNodes);
	m_mesh->freeMark<Node>(markPointNodes);

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
