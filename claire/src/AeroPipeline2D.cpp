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
  AbstractAeroPipeline(Aparams)
{
	m_meshTet = new Mesh(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
                                       gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));
	m_meshHex = new Mesh(gmds::MeshModel(gmds::DIM3|gmds::F|gmds::N|gmds::E| gmds::N2E|
                                       gmds::N2F|gmds::F2N|gmds::E2N|gmds::F2E|gmds::E2F));
	m_Bnd = new AeroBoundaries_2D(m_meshTet) ;
	m_couche_id = m_meshHex->newVariable<int, GMDS_NODE>("GMDS_Couche_Id");
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
AbstractAeroPipeline::STATUS AeroPipeline2D::execute(){

	LectureMaillage();
	m_Bnd->execute();

	m_manager->initAndLinkFrom2DMesh(m_meshTet, m_linker_TG);

	// Calcul du level set
	m_meshTet->newVariable<double,GMDS_NODE>("GMDS_Distance");
	m_meshTet->newVariable<double,GMDS_NODE>("GMDS_Distance_Int");
	m_meshTet->newVariable<double,GMDS_NODE>("GMDS_Distance_Out");
	LevelSetCombined lsCombined(m_meshTet, m_Bnd->getMarkParoi(), m_Bnd->getMarkAmont(),
	                            m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance"),
	                            m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
	                            m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance_Out"));
	lsCombined.execute();

	// Calcul du gradient du champ de Level Set
	m_meshTet->newVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient");
	LeastSquaresGradientComputation grad2D(m_meshTet, m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
	                                       m_meshTet->getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	grad2D.execute();

	for (int i=1;i<=m_Bnd->getNbrBords();i++){
		if(i != m_Bnd->getColorAmont()){
			DiscretisationParoi(i);
		}
	}

	//m_linker_QHG.setGeometry(&m_manager);
	//m_linker_QHG.setMesh(&m_mQuad);

	GenerationCouche(1, 0.25);
	GenerationCouche(2, 1);

	ConvertisseurMeshToBlocking();

	EcritureMaillage();

	//std::cout << "Nbr faces : " << m_mQuad.getNbFaces() << std::endl;
	//std::cout << "Nbr arêtes : " << m_mQuad.getNbEdges() << std::endl;

	return AbstractAeroPipeline::SUCCESS;

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline2D::LectureMaillage(){

	// Lecture du maillage
	std::cout << "-> Lecture du maillage ..." << std::endl;

	gmds::IGMeshIOService ioService(m_meshTet);
	gmds::VTKReader vtkReader(&ioService);
	vtkReader.setCellOptions(gmds::N|gmds::F);
	vtkReader.read(m_params.input_file);

	gmds::MeshDoctor doc(m_meshTet);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	math::Utils::MeshCleaner(m_meshTet);

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline2D::EcritureMaillage(){

	std::cout << "-> Ecriture du maillage ..." << std::endl;

	//gmds::IGMeshIOService ioService(m_meshHex);
	gmds::IGMeshIOService ioService(&m_Blocking2D);
	gmds::VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(gmds::N|gmds::F);
	vtkWriter.setDataOptions(gmds::N|gmds::F);
	std::string dir(".");//TEST_SAMPLES_DIR);
	vtkWriter.write(m_params.output_file);

	// Ecriture du maillage en triangles initial pour visualisation et débug
	ioService = m_meshTet;
	gmds::VTKWriter vtkWriter2(&ioService);
	vtkWriter2.setCellOptions(gmds::N|gmds::F);
	vtkWriter2.setDataOptions(gmds::N|gmds::F);
	vtkWriter2.write("AeroPipeline2D_Triangles.vtk");



}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline2D::GenerationCouche(int couche_id, double dist){

	std::cout << "-> Génération de la couche " << couche_id << std::endl ;

	std::vector<Node> nodes_couche_id;	// Noeuds sur la douche couche_id

	// Création des nouveaux noeuds de la couche couche_id naïvement
	// + aretes entre la couche couche_id-1 et couche_id
	for (auto n_id:m_meshHex->nodes()){
		if (m_couche_id->value(n_id) == couche_id-1){
			// Placement du point P à la distance souhaitée suivant le champ de gradient
			Node n = m_meshHex->get<Node>(n_id);
			math::Point M = n.point();
			AdvectedPointRK4_2D advpoint(m_meshTet, M, dist, m_meshTet->getVariable<double, GMDS_NODE>("GMDS_Distance"),
			   m_meshTet->getVariable<math::Vector3d ,GMDS_NODE>("GMDS_Gradient"));
			advpoint.execute();
			math::Point P = advpoint.getPend();

			// Nouveau noeud dans le maillage
			Node n_new = m_meshHex->newNode(P);
			m_couche_id->set(n_new.id(),couche_id);
			nodes_couche_id.push_back(n_new);

			// Création de l'arête pour relier le noeud à la couche n-1
			Edge e_new = m_meshHex->newEdge(n, n_new);
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
	Edge e3 = m_meshHex->newEdge(n0, n3);
	n0.add<Edge>(e3);	// Connectivités N->E
	n3.add<Edge>(e3);

	// Création de la face de type quad
	Face f = m_meshHex->newQuad(n0, n1, n2, n3) ;
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
void AeroPipeline2D::DiscretisationParoi(int color){

	std::cout << "-> Discrétisation du bord de couleur " << color << std::endl;

	std::vector<TCellID> bnd_nodes_id_ordered = m_Bnd->BndNodesOrdered(color);	// Tri les noeuds du bord concerné dans l'ordre

	int markCurveEdges = m_meshTet->newMark<Edge>();
	int markCurveNodes = m_meshTet->newMark<Node>();
	int markPointNodes = m_meshTet->newMark<Node>();
	int markAloneNodes = m_meshTet->newMark<Node>();

	BoundaryOperator2D bnd_op(m_meshTet);
	bnd_op.markCellOnGeometry(markCurveEdges, markCurveNodes,
	                          markPointNodes, markAloneNodes);

	m_meshTet->unmarkAll<Node>(markAloneNodes);
	m_meshTet->freeMark<Node>(markAloneNodes);
	m_meshTet->unmarkAll<Edge>(markCurveEdges);
	m_meshTet->freeMark<Edge>(markCurveEdges);
	m_meshTet->unmarkAll<Node>(markCurveNodes);
	m_meshTet->freeMark<Node>(markCurveNodes);

	// Calcul de la longueur max des arêtes de bloc
	double perim = m_Bnd->ComputeBoundaryLength(color);		// Calcul du périmètre du bord regardé
	double Lmax = perim/ m_params.nbrMinBloc ;
	double l(0);

	// Calcul de la valeur min de x sur les noeuds du bord
	TCellID n_arret_tri_id = m_Bnd->PointArret(color);
	Node n_arret_tri = m_meshTet->get<Node>(n_arret_tri_id) ;
	math::Point point_arret = n_arret_tri.point();
	double x_min = point_arret.X() ;

	Node n0_tri = m_meshTet->get<Node>(bnd_nodes_id_ordered[0]) ;
	Node n0_quad = m_meshHex->newNode(n0_tri.point()); // Premier noeud du nouveau maillage

	Node n1_quad = n0_quad;
	Node n2_quad = n0_quad;

	for(int i=1;i<bnd_nodes_id_ordered.size();i++){

		Node n = m_meshTet->get<Node>(bnd_nodes_id_ordered[i]);
		math::Point p = n.point() ;
		l += math::Utils::distFromNodeIds(m_meshTet, bnd_nodes_id_ordered[i], bnd_nodes_id_ordered[i-1]);

		Node n_gauche = m_meshTet->get<Node>(bnd_nodes_id_ordered[i-1]);
		Node n_droit;
		if(i == bnd_nodes_id_ordered.size()-1) {
			n_droit = m_meshTet->get<Node>(bnd_nodes_id_ordered[0]);
		}
		else{
			n_droit = m_meshTet->get<Node>(bnd_nodes_id_ordered[i+1]);
		}

		bool isExtremum(false);
		if ( (n.X() < n_gauche.X() && n.X() < n_droit.X() )
		    || (n.X() > n_gauche.X() && n.X() > n_droit.X() )
		    || (n.Y() < n_gauche.Y() && n.Y() < n_droit.Y() )
		    || (n.Y() > n_gauche.Y() && n.Y() > n_droit.Y() )) {
			isExtremum = true;
		}



		if ( m_meshTet->isMarked<Node>(bnd_nodes_id_ordered[i], markPointNodes)
		    || l >= Lmax
		    || abs(l-Lmax) <= pow(10,-6)
		    || abs(p.X() - x_min) <= pow(10,-6)
		    || isExtremum ){

			n1_quad = n2_quad;
			n2_quad = m_meshHex->newNode(n.point());
			Edge e = m_meshHex->newEdge(n1_quad, n2_quad);

			// Ajout des connectivités Node -> Edge
			n1_quad.add<Edge>(e);
			n2_quad.add<Edge>(e);

			l = 0;	// On remet l à 0

		}

	}

	// Ajout de la dernière arête pour fermer le blocking du bord
	Edge e = m_meshHex->newEdge(n2_quad, n0_quad);
	// Ajout des connectivités Node -> Edge
	n2_quad.add<Edge>(e);
	n0_quad.add<Edge>(e);

	m_meshTet->unmarkAll<Node>(markPointNodes);
	m_meshTet->freeMark<Node>(markPointNodes);

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void AeroPipeline2D::ConvertisseurMeshToBlocking(){

	Variable<TCellID>* var_new_id = m_meshHex->newVariable<TCellID,GMDS_NODE>("New_ID");

	for (auto n_id:m_meshHex->nodes()){
		Node n = m_meshHex->get<Node>(n_id);
		Node n_blocking = m_Blocking2D.newBlockCorner(n.point());
		var_new_id->set(n_id, n_blocking.id());
	}

	for (auto f_id:m_meshHex->faces()){
		Face f = m_meshHex->get<Face>(f_id);
		std::vector<Node> quad_nodes = f.get<Node>() ;
		Blocking2D::Block B0 = m_Blocking2D.newBlock( var_new_id->value(quad_nodes[0].id()), var_new_id->value(quad_nodes[1].id()),
		                      var_new_id->value(quad_nodes[2].id()), var_new_id->value(quad_nodes[3].id()));
		B0.seNbDiscretizationI(10);
		B0.seNbDiscretizationJ(10);
	}

	m_Blocking2D.initializeGridPoints();

	m_meshHex->deleteVariable(GMDS_NODE, "New_ID");

}
/*------------------------------------------------------------------------*/


