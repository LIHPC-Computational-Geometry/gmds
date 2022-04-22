//
// Created by rochec on 09/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/Utils.h>
#include <gmds/claire/AeroPipeline_2D.h>
#include <gmds/claire/AeroBoundaries_2D.h>
#include <gmds/claire/LevelSetCombined.h>
#include <gmds/claire/LeastSquaresGradientComputation.h>
#include <gmds/claire/AdvectedPointRK4_2D.h>
#include <gmds/claire/Grid_Smooth2D.h>
#include <gmds/claire/AeroExtrusion_2D.h>

#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/igalgo/BoundaryOperator2D.h>
#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <iostream>
#include <chrono>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AeroPipeline_2D::AeroPipeline_2D(ParamsAero Aparams) :
  AbstractAeroPipeline(Aparams),
  m_linker_BG(new cad::GeomMeshLinker())
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
AbstractAeroPipeline::STATUS
AeroPipeline_2D::execute(){

	clock_t t_start, t_end;

	LectureMaillage();

	std::cout << "-> Analyse des bords" << std::endl;
	t_start = clock();
	m_Bnd->execute();
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;
	std::cout << " " << std::endl;

	m_manager->initAndLinkFrom2DMesh(m_meshTet, m_linker_TG);

	// Calcul du level set
	std::cout << "-> Calcul Level Sets" << std::endl;
	t_start = clock();
	m_meshTet->newVariable<double,GMDS_NODE>("GMDS_Distance");
	m_meshTet->newVariable<double,GMDS_NODE>("GMDS_Distance_Int");
	m_meshTet->newVariable<double,GMDS_NODE>("GMDS_Distance_Out");
	LevelSetCombined lsCombined(m_meshTet, m_Bnd->getMarkParoi(), m_Bnd->getMarkAmont(),
	                            m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance"),
	                            m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
	                            m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance_Out"));
	lsCombined.execute();
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;
	std::cout << " " << std::endl;

	// Calcul du gradient du champ de Level Set
	std::cout << "-> Calcul Gradient" << std::endl;
	t_start = clock();
	m_meshTet->newVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient");
	LeastSquaresGradientComputation grad2D(m_meshTet, m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
	                                       m_meshTet->getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	grad2D.execute();
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;
	std::cout << " " << std::endl;

	m_linker_HG->setGeometry(m_manager);
	m_linker_HG->setMesh(m_meshHex);

	std::cout << "-> Discrétisation du bord" << std::endl;
	t_start = clock();
	for (int i=1;i<=m_Bnd->getNbrBords();i++){
		if(i != m_Bnd->getColorAmont()){
			DiscretisationParoi(i);
		}
	}
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;
	std::cout << " " << std::endl;

	// Extrusion par couches
	/*
	m_nbr_couches=1;
	for(int couche=1;couche<=m_nbr_couches;couche++){
		t_start = clock();
		GenerationCouche(couche, 1.0/m_nbr_couches);
		t_end = clock();
		std::cout << "........................................ temps : " << 1.0*(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;
	}
	 */

	std::cout << "-> Extrusion" << std::endl;
	t_start = clock();
	AeroExtrusion_2D aero_extrusion(m_meshTet, m_meshHex);
	aero_extrusion.execute();
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;

	// Mise à jour du linker pour la dernière couche de noeuds
	UpdateLinkerLastLayer(1);

	// Conversion structure maillage à blocking + maillage des blocs par transfinies
	ConvertisseurMeshToBlocking();

	// Lissage
	/*
	std::cout << "-> Lissage final" << std::endl;
	t_start = clock();
	Grid_Smooth2D smoother(&m_Blocking2D, 100);
	smoother.execute();
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;
	 */


	// Ecriture finale des maillages
	EcritureMaillage();

	std::cout << " " << std::endl;
	std::cout << "INFORMATIONS COMPLEMENTAIRES :" << std::endl;
	std::cout << "Nbr faces : " << m_meshHex->getNbFaces() << std::endl;
	std::cout << "Nbr arêtes : " << m_meshHex->getNbEdges() << std::endl;

	return AbstractAeroPipeline::SUCCESS;

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void
AeroPipeline_2D::LectureMaillage(){

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
void
AeroPipeline_2D::EcritureMaillage(){

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
void
AeroPipeline_2D::DiscretisationParoi(int color){

	//std::cout << "-> Discrétisation du bord de couleur " << color << std::endl;

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
	//TCellID n_arret_tri_id = m_Bnd->PointArret(color);
	//Node n_arret_tri = m_meshTet->get<Node>(n_arret_tri_id) ;
	//math::Point point_arret = n_arret_tri.point();
	//double x_min = point_arret.X() ;

	Node n0_tri = m_meshTet->get<Node>(bnd_nodes_id_ordered[0]) ;
	Node n0_quad = m_meshHex->newNode(n0_tri.point()); // Premier noeud du nouveau maillage

	UpdateLinker(m_linker_TG, n0_tri, m_linker_HG, n0_quad);

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
		    //|| abs(p.X() - x_min) <= pow(10,-6)
		    || isExtremum ){

			n1_quad = n2_quad;
			n2_quad = m_meshHex->newNode(n.point());
			UpdateLinker(m_linker_TG, n, m_linker_HG, n2_quad);
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
void
AeroPipeline_2D::GenerationCouche(int couche_id, double dist){

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
Node
AeroPipeline_2D::SuccessorNode(Node n0){
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
Node
AeroPipeline_2D::AnteriorNode(Node n0){
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
std::vector<Node>
AeroPipeline_2D::AdjNodesInLayer(Node n0, int couche_i){
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
void
AeroPipeline_2D::CreateQuadAndConnectivities(Node n0, Node n1, Node n2){

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
bool
AeroPipeline_2D::isQuadCreated(Node n0, Node n1, Node n2){
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
void
AeroPipeline_2D::ConvertisseurMeshToBlocking(){

	std::cout << "-> Conversion Maillage en Blocking2D et connectivités" << std::endl;

	Variable<TCellID>* var_new_id = m_meshHex->newVariable<TCellID,GMDS_NODE>("New_ID");

	m_linker_BG->setGeometry(m_manager);
	m_linker_BG->setMesh(&m_Blocking2D);

	Variable<int>* var_couche = m_Blocking2D.newVariable<int, GMDS_NODE>("GMDS_Couche");

	for (auto n_id:m_meshHex->nodes()){
		Node n = m_meshHex->get<Node>(n_id);
		Node n_blocking = m_Blocking2D.newBlockCorner(n.point());
		var_new_id->set(n_id, n_blocking.id());
		var_couche->set(n_blocking.id(), m_couche_id->value(n_id));

		UpdateLinker(m_linker_HG, n, m_linker_BG, n_blocking); // Init linker_BG
	}

	for (auto f_id:m_meshHex->faces()){
		Face f = m_meshHex->get<Face>(f_id);
		std::vector<Node> quad_nodes = f.get<Node>() ;
		Blocking2D::Block B0 = m_Blocking2D.newBlock( var_new_id->value(quad_nodes[0].id()), var_new_id->value(quad_nodes[1].id()),
		                      var_new_id->value(quad_nodes[2].id()), var_new_id->value(quad_nodes[3].id()));
		B0.setNbDiscretizationI(10);
		B0.setNbDiscretizationJ(10);
	}

	m_Blocking2D.initializeGridPoints();	// Maillage des blocs par transfinies

	m_meshHex->deleteVariable(GMDS_NODE, "New_ID");

	BlockingClassification();

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void
AeroPipeline_2D::UpdateLinker(cad::GeomMeshLinker* linker_1, Node n_1, cad::GeomMeshLinker* linker_2, Node n_2){
	int geom_dim = linker_1->getGeomDim<Node>(n_1.id());
	if(geom_dim == 1){
		linker_2->linkNodeToPoint(n_2.id(), linker_1->getGeomId<Node>(n_1.id()));
	}
	else if(geom_dim==2){
		linker_2->linkNodeToCurve(n_2.id(), linker_1->getGeomId<Node>(n_1.id()));
	}
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void
AeroPipeline_2D::UpdateLinkerLastLayer(int layer_id){

	for (auto n_id:m_meshHex->nodes()){
		if( m_couche_id->value(n_id) == layer_id ) {
			Node n = m_meshHex->get<Node>(n_id);
			math::Point p = n.point();
			TCellID closest_n_id = m_Bnd->ClosestNodeOnBnd(m_Bnd->getColorAmont(), p);
			UpdateLinker(m_linker_TG, m_meshTet->get<Node>(closest_n_id), m_linker_HG, n);
			Node n_test = m_meshTet->get<Node>(closest_n_id);
		}
	}

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void
AeroPipeline_2D::BlockingClassification(){

	Variable<int>* var_couche = m_Blocking2D.getVariable<int, GMDS_NODE>("GMDS_Couche");

	for (auto B0:m_Blocking2D.allBlocks()){
		int Nx = B0.getNbDiscretizationI()-1;
		int Ny = B0.getNbDiscretizationJ()-1;

		if ( var_couche->value( B0(0,0).id() ) == var_couche->value( B0(Nx,0).id() ) ) {
			int geom_id_corner_1 = m_linker_BG->getGeomId<Node>(B0(0,0).id()) ;
			int geom_id_corner_2 = m_linker_BG->getGeomId<Node>(B0(Nx,0).id()) ;
			int dim_corner_1 = m_linker_BG->getGeomId<Node>(B0(0,0).id()) ;
			int dim_corner_2 = m_linker_BG->getGeomId<Node>(B0(Nx,0).id());
			if( dim_corner_1 == dim_corner_2 && dim_corner_1 == 2 )
			{
				for (int i=1;i<Nx;i++){
					Node n = B0(i,0);
					m_linker_BG->linkNodeToCurve(n.id(), geom_id_corner_1);
					cad::GeomCurve* curve = m_manager->getCurve(geom_id_corner_1);
					math::Point p = n.point();
					curve->project(p);
					n.setPoint(p);
				}
			}
		}

		if ( var_couche->value( B0(0,0).id() ) == var_couche->value( B0(0,Ny).id() ) ) {
			int geom_id_corner_1 = m_linker_BG->getGeomId<Node>(B0(0,0).id()) ;
			int geom_id_corner_2 = m_linker_BG->getGeomId<Node>(B0(0,Ny).id()) ;
			int dim_corner_1 = m_linker_BG->getGeomDim<Node>(B0(0,0).id()) ;
			int dim_corner_2 = m_linker_BG->getGeomDim<Node>(B0(0,Ny).id());
			//std::cout << "dim 1 : " << dim_corner_1 << std::endl;
			//std::cout << "dim 2 : " << dim_corner_2 << std::endl;
			if( dim_corner_1 == dim_corner_2 && dim_corner_1 == 2 )
			{
				for (int j=1;j<Ny;j++){
					Node n = B0(0,j);
					m_linker_BG->linkNodeToCurve(n.id(), geom_id_corner_1);
					cad::GeomCurve* curve = m_manager->getCurve(geom_id_corner_1);
					math::Point p = n.point();
					curve->project(p);
					n.setPoint(p);
				}
			}
		}

		if ( var_couche->value( B0(0,Ny).id() ) == var_couche->value( B0(Nx,Ny).id() ) ) {
			int geom_id_corner_1 = m_linker_BG->getGeomId<Node>(B0(0,Ny).id()) ;
			int geom_id_corner_2 = m_linker_BG->getGeomId<Node>(B0(Nx,Ny).id()) ;
			int dim_corner_1 = m_linker_BG->getGeomId<Node>(B0(0,Ny).id()) ;
			int dim_corner_2 = m_linker_BG->getGeomId<Node>(B0(Nx,Ny).id());
			if( dim_corner_1 == dim_corner_2 && dim_corner_1 == 2 )
			{
				for (int i=1;i<Nx;i++){
					Node n = B0(i,Ny);
					m_linker_BG->linkNodeToCurve(n.id(), geom_id_corner_1);
					cad::GeomCurve* curve = m_manager->getCurve(geom_id_corner_1);
					math::Point p = n.point();
					curve->project(p);
					n.setPoint(p);
				}
			}
		}

		if ( var_couche->value( B0(Nx,0).id() ) == var_couche->value( B0(Nx,Ny).id() ) ) {
			int geom_id_corner_1 = m_linker_BG->getGeomId<Node>(B0(Nx,0).id()) ;
			int geom_id_corner_2 = m_linker_BG->getGeomId<Node>(B0(Nx,Ny).id()) ;
			int dim_corner_1 = m_linker_BG->getGeomId<Node>(B0(Nx,0).id()) ;
			int dim_corner_2 = m_linker_BG->getGeomId<Node>(B0(Nx,Ny).id());
			if( dim_corner_1 == dim_corner_2 && dim_corner_1 == 2 )
			{
				for (int j=1;j<Ny;j++){
					Node n = B0(Nx,j);
					m_linker_BG->linkNodeToCurve(n.id(), geom_id_corner_1);
					cad::GeomCurve* curve = m_manager->getCurve(geom_id_corner_1);
					math::Point p = n.point();
					curve->project(p);
					n.setPoint(p);
				}
			}
		}

	}

}
/*------------------------------------------------------------------------*/
