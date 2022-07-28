//
// Created by rochec on 09/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/Utils.h>
#include <gmds/claire/AeroPipeline_2D.h>
#include <gmds/claire/AeroBoundaries_2D.h>
#include <gmds/claire/LevelSetExtended.h>
#include <gmds/claire/LevelSetCombined.h>
#include <gmds/claire/LeastSquaresGradientComputation.h>
#include <gmds/claire/AdvectedPointRK4_2D.h>
#include <gmds/claire/Grid_Smooth2D.h>
#include <gmds/claire/SmoothLineSweepingYao.h>
#include <gmds/claire/SmoothLineSweepingOrtho.h>
#include <gmds/claire/AeroExtrusion_2D.h>
#include <gmds/claire/SU2Writer.h>
#include <gmds/claire/IntervalAssignment_2D.h>
#include <gmds/math/TransfiniteInterpolation.h>
#include <gmds/claire/RefinementBeta.h>

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
	std::cout << "Nombre de bords : " << m_Bnd->getNbrBords()-1 << std::endl;
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

	/*

   Variable<double> * var_d_int = m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance_Int");
	Variable<double> * var_d_comb = m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance");
	double d_max_amont(0.0);
	for (auto n_id:m_meshTet->nodes()){
		Node n = m_meshTet->get<Node>(n_id);
		math::Point p = n.point();
		if ( p.X() < m_params.x_lim && d_max_amont < var_d_int->value(n_id) ) {
         d_max_amont = var_d_int->value(n_id);
		}
	}

	int m_mark_z1 = m_meshTet->newMark<Node>();
	for (auto n_id:m_meshTet->nodes()){
		if ( var_d_int->value(n_id) > d_max_amont || abs(var_d_comb->value(n_id)-1.0) < pow(10,-6) )
		{
			Node n = m_meshTet->get<Node>(n_id);
         m_meshTet->mark(n, m_mark_z1);
		}
	}

	//LevelSetExtended ls_test(m_meshTet, m_mark_z1,m_meshTet->newVariable<double,GMDS_NODE>("GMDS_Distance_Z1"));
	//ls_test.execute();

	LevelSetCombined lsCombined2(m_meshTet, m_Bnd->getMarkParoi(), m_mark_z1,
	                            m_meshTet->newVariable<double,GMDS_NODE>("GMDS_Distance_2"),
	                            m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
	                             m_meshTet->newVariable<double,GMDS_NODE>("GMDS_Distance_Z1"));
	lsCombined2.execute();

	Variable<double> * var_d_comb_2 = m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance_2");
	Variable<double> * var_d_comb_3 = m_meshTet->newVariable<double,GMDS_NODE>("GMDS_Distance_3");

	for (auto n_id:m_meshTet->nodes())
	{
		var_d_comb_3->set(n_id, var_d_comb_2->value(n_id));
	}

	for (auto f_id:m_meshTet->faces())
	{
		Face f = m_meshTet->get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();

		if ( abs( var_d_comb_2->value(nodes[0].id()) - 1.0) < pow(10,-6)
		    && abs( var_d_comb_2->value(nodes[1].id()) - 1.0) < pow(10,-6)
		    && abs( var_d_comb_2->value(nodes[2].id()) - 1.0) < pow(10,-6) )
		{
			var_d_comb_3->set(nodes[0].id(), 1.1);
			var_d_comb_3->set(nodes[1].id(), 1.1);
			var_d_comb_3->set(nodes[2].id(), 1.1);
		}

	}

	m_meshTet->unmarkAll<Node>(m_mark_z1);
	m_meshTet->freeMark<Node>(m_mark_z1);


	*/


	// Calcul du gradient du champ de Level Set
	std::cout << "-> Calcul Gradient" << std::endl;
	t_start = clock();
	//m_meshTet->newVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient");
	//LeastSquaresGradientComputation grad2D(m_meshTet, m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
	//                                       m_meshTet->getVariable<math::Vector3d, GMDS_NODE>("GMDS_Gradient"));
	//grad2D.execute();
	Variable<math::Vector3d>* var_VectorField = m_meshTet->newVariable<math::Vector3d, GMDS_NODE>("VectorField_Extrusion");
	ComputeVectorFieldForExtrusion();
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
	AeroExtrusion_2D aero_extrusion(m_meshTet, m_meshHex, m_params, var_VectorField);
	aero_extrusion.execute();
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;
	std::cout << " " << std::endl;

	std::cout << "-> Classification géométrique" << std::endl;
	t_start = clock();
	// Mise à jour du linker pour la dernière couche de noeuds
	UpdateLinkerLastLayer();
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;
	std::cout << " " << std::endl;

	std::cout << "-> Conversion maillage en blocking" << std::endl;
	t_start = clock();
	// Conversion structure maillage à blocking + maillage des blocs par transfinies
	ConvertisseurMeshToBlocking();
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;
	std::cout << " " << std::endl;

	// Lissage
	std::cout << "-> Lissage final" << std::endl;
	t_start = clock();
	//Grid_Smooth2D smoother(&m_Blocking2D, 2000);
	//smoother.execute();

	for (auto b:m_Blocking2D.allBlocks())
	{
		SmoothLineSweepingYao smoother( &b, 20, 0.5);
		smoother.execute();
	}

	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;
	std::cout << " " << std::endl;



	// Add the Beta Refinement
	for (auto b:m_Blocking2D.allBlocks())
	{
		int Nx = b.getNbDiscretizationI();
		int Ny = b.getNbDiscretizationJ();

		Variable<int>* var_couche = m_Blocking2D.getVariable<int, GMDS_NODE>("GMDS_Couche");
		Node n0 = b.getNode(0);
		Node n1 = b.getNode(1);
		Node n2 = b.getNode(2);
		Node n3 = b.getNode(3);

		/*
		std::cout << "---------------" << std::endl;
		std::cout << "Noeud 0 dans la couche : " << var_couche->value(n0.id()) << std::endl;
		std::cout << "Noeud 1 dans la couche : " << var_couche->value(n1.id()) << std::endl;
		std::cout << "Noeud 2 dans la couche : " << var_couche->value(n2.id()) << std::endl;
		std::cout << "Noeud 3 dans la couche : " << var_couche->value(n3.id()) << std::endl;
		 */

		int compteur(0);
		if (var_couche->value(n0.id()) == 0)
		{
			compteur++;
		}
		if (var_couche->value(n1.id()) == 0)
		{
			compteur++;
		}
		if (var_couche->value(n2.id()) == 0)
		{
			compteur++;
		}
		if (var_couche->value(n3.id()) == 0)
		{
			compteur++;
		}

		if (compteur == 2) {
			for (int i = 0; i < Nx; i++) {
				std::vector<math::Point> Points;
				for (int j = 0; j < Ny; j++) {
					Points.push_back(b(i, j).point());
				}

				RefinementBeta ref(Points, 8*pow(10, -6));
				ref.execute();
				Points = ref.GetNewPositions();

				for (int j = 1; j < Ny - 1; j++) {
					b(i, j).setPoint({Points[j]});
				}
			}
		}
	}




	// Compute the quality criterions of the blocking
	math::Utils::AnalyseQuadMeshQuality(&m_Blocking2D);

	// Ecriture finale des maillages
	std::cout << "-> Ecriture finale des maillages" << std::endl;
	t_start = clock();
	EcritureMaillage();
	t_end = clock();
	std::cout << "........................................ temps : " << 1.0*(t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;

	std::cout << " " << std::endl;
	std::cout << "INFORMATIONS COMPLEMENTAIRES :" << std::endl;
	std::cout << "Nbr de blocs : " << m_Blocking2D.getNbFaces() << std::endl;
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
	std::cout << "Ecriture 1 ..." << std::endl;

	//gmds::IGMeshIOService ioService(m_meshHex);
	gmds::IGMeshIOService ioService(&m_Blocking2D);
	gmds::VTKWriter vtkWriter_Blocking(&ioService);
	vtkWriter_Blocking.setCellOptions(gmds::N|gmds::F);
	vtkWriter_Blocking.setDataOptions(gmds::N|gmds::F);
	std::string dir(".");//TEST_SAMPLES_DIR);
	//vtkWriter_Blocking.write(m_params.output_file);
	vtkWriter_Blocking.write("AeroPipeline2D_Blocking.vtk");

	std::cout << "Ecriture 2 ..." << std::endl;

	math::Utils::BuildMesh2DFromBlocking2D(&m_Blocking2D, m_meshHex);
	math::Utils::AnalyseQuadMeshQuality(m_meshHex);
	ioService = m_meshHex;
	gmds::VTKWriter vtkWriter_HexMesh(&ioService);
	vtkWriter_HexMesh.setCellOptions(gmds::N|gmds::F);
	vtkWriter_HexMesh.setDataOptions(gmds::N|gmds::F);
	vtkWriter_HexMesh.write("AeroPipeline2D_HexMesh.vtk");

	std::cout << "Ecriture 3 ..." << std::endl;

	// Ecriture du maillage en triangles initial pour visualisation et débug
	ioService = m_meshTet;
	gmds::VTKWriter vtkWriter_TetMesh(&ioService);
	vtkWriter_TetMesh.setCellOptions(gmds::N|gmds::F);
	vtkWriter_TetMesh.setDataOptions(gmds::N|gmds::F);
	vtkWriter_TetMesh.write("AeroPipeline2D_TetMesh.vtk");

	// Ecriture du maillage au format su2
	SU2Writer writer(m_meshHex, "AeroPipeline2D_QuadMesh.su2", m_params.x_lim_SU2_inoutlet);
	SU2Writer::STATUS result = writer.execute();

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

	// Classification géométrique des arêtes de la paroi

	for (auto e_id:m_meshHex->edges())
	{
		Edge e = m_meshHex->get<Edge>(e_id);
		std::vector<Node> nodes = e.get<Node>();

		int geom_id_node_0 = m_linker_HG->getGeomId<Node>(nodes[0].id()) ;
		int geom_id_node_1 = m_linker_HG->getGeomId<Node>(nodes[1].id()) ;
		int dim_node_0 = m_linker_HG->getGeomDim<Node>(nodes[0].id()) ;
		int dim_node_1 = m_linker_HG->getGeomDim<Node>(nodes[1].id());

		if (dim_node_0 == 2)
		{
			m_linker_HG->linkNodeToCurve(e_id, geom_id_node_0);
		}
		else if (dim_node_1 == 2)
		{
			m_linker_HG->linkNodeToCurve(e_id, geom_id_node_1);
		}
		else if (dim_node_0 == 1 && dim_node_1 == 1)
		{
			cad::GeomPoint* PG_0 = m_manager->getPoint(geom_id_node_0);
			cad::GeomPoint* PG_1 = m_manager->getPoint(geom_id_node_1);
			//cad::GeomCurve* Curve = PG_0.co ;
			std::vector<cad::GeomCurve*> PG_0_Curves = PG_0->curves() ;
			std::vector<cad::GeomCurve*> PG_1_Curves = PG_1->curves() ;
			if (PG_0_Curves[0]->id() == PG_1_Curves[0]->id()
			    || PG_0_Curves[0]->id() == PG_1_Curves[1]->id() )
			{
				m_linker_HG->linkNodeToCurve(e_id, PG_0_Curves[0]->id());
			}
			else if (PG_0_Curves[1]->id() == PG_1_Curves[0]->id()
			         || PG_0_Curves[1]->id() == PG_1_Curves[1]->id() )
			{
				m_linker_HG->linkNodeToCurve(e_id, PG_0_Curves[1]->id());
			}

		}

	}

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void
AeroPipeline_2D::ConvertisseurMeshToBlocking(){

	//std::cout << "-> Conversion Maillage en Blocking2D et connectivités" << std::endl;

	m_linker_BG->setGeometry(m_manager);
	m_linker_BG->setMesh(&m_Blocking2D);

	//m_Blocking2D = Blocking2D(*m_meshHex);

	Variable<TCellID>* var_new_id = m_meshHex->newVariable<TCellID,GMDS_NODE>("New_ID");
	Variable<int>* var_couche = m_Blocking2D.newVariable<int, GMDS_NODE>("GMDS_Couche");

	for (auto n_id:m_meshHex->nodes()){
		Node n = m_meshHex->get<Node>(n_id);
		Node n_blocking = m_Blocking2D.newBlockCorner(n.point());
		var_new_id->set(n_id, n_blocking.id());
		var_couche->set(n_blocking.id(), m_couche_id->value(n_id));
		//std::cout << "Noeud " << n_blocking.id() << " dans la couche " << var_couche->value(n_blocking.id()) << std::endl;

		UpdateLinker(m_linker_HG, n, m_linker_BG, n_blocking); // Init linker_BG
	}

	for (auto f_id:m_meshHex->faces()){
		Face f = m_meshHex->get<Face>(f_id);
		std::vector<Node> quad_nodes = f.get<Node>() ;
		Blocking2D::Block B0 = m_Blocking2D.newBlock( var_new_id->value(quad_nodes[0].id()), var_new_id->value(quad_nodes[1].id()),
		                      var_new_id->value(quad_nodes[2].id()), var_new_id->value(quad_nodes[3].id()));
		//B0.setNbDiscretizationI(10);
		//B0.setNbDiscretizationJ(10);
	}

	IntervalAssignment_2D IntAss(&m_Blocking2D, m_params);
	IntAss.execute();

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
AeroPipeline_2D::UpdateLinkerLastLayer(){

	int max_layer_id(0);
	Variable<int>* couche_id = m_meshHex->getVariable<int, GMDS_NODE>("GMDS_Couche_Id");

	// Compute the id of the last layer
	for (auto n_id:m_meshHex->nodes())
	{
		if (couche_id->value(n_id) > max_layer_id)
		{
			max_layer_id = couche_id->value(n_id);
		}
	}

	for (auto n_id:m_meshHex->nodes()){
		if( m_couche_id->value(n_id) == max_layer_id ) {
			Node n = m_meshHex->get<Node>(n_id);
			math::Point p = n.point();
			TCellID closest_n_id = m_Bnd->ClosestNodeOnBnd(m_Bnd->getColorAmont(), p);
			UpdateLinker(m_linker_TG, m_meshTet->get<Node>(closest_n_id), m_linker_HG, n);
			//Node n_test = m_meshTet->get<Node>(closest_n_id);

			int geom_id = m_linker_HG->getGeomId<Node>(n.id()) ;
			int geom_dim = m_linker_HG->getGeomDim<Node>(n.id()) ;
			if (geom_dim == 2)
			{
				cad::GeomCurve* curve = m_manager->getCurve(geom_id);
				curve->project(p);
				n.setPoint(p);
			}
		}
	}

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void
AeroPipeline_2D::BlockingClassification(){

	Variable<int>* var_couche = m_Blocking2D.getVariable<int, GMDS_NODE>("GMDS_Couche");

	for (auto B0:m_Blocking2D.allBlocks())
	{
		int Nx = B0.getNbDiscretizationI()-1;
		int Ny = B0.getNbDiscretizationJ()-1;

		if ( var_couche->value( B0(0,0).id() ) == var_couche->value( B0(Nx,0).id() )) {
			int geom_id_corner_1 = m_linker_BG->getGeomId<Node>(B0(0,0).id()) ;
			int geom_id_corner_2 = m_linker_BG->getGeomId<Node>(B0(Nx,0).id()) ;
			int dim_corner_1 = m_linker_BG->getGeomDim<Node>(B0(0,0).id()) ;
			int dim_corner_2 = m_linker_BG->getGeomDim<Node>(B0(Nx,0).id());
			if( dim_corner_1 == dim_corner_2 && dim_corner_1 == 2 )
			{
				for (int i=1;i<Nx;i++)
				{
					Node n = B0(i,0);
					m_linker_BG->linkNodeToCurve(n.id(), geom_id_corner_1);
					cad::GeomCurve* curve = m_manager->getCurve(geom_id_corner_1);
					math::Point p0 = n.point();
					math::Point p = n.point();
					curve->project(p);
					n.setPoint(p);


					if (var_couche->value( B0(0,0).id() ) == 0) {
						math::Vector3d delta_p = p - p0;
						Node n_opp = B0(i, Ny);
						math::Point p_opp = n_opp.point();
						n_opp.setPoint(p_opp + delta_p);
					}

				}
			}
		}

		if ( var_couche->value( B0(0,0).id() ) == var_couche->value( B0(0,Ny).id() )) {
			int geom_id_corner_1 = m_linker_BG->getGeomId<Node>(B0(0,0).id()) ;
			int geom_id_corner_2 = m_linker_BG->getGeomId<Node>(B0(0,Ny).id()) ;
			int dim_corner_1 = m_linker_BG->getGeomDim<Node>(B0(0,0).id()) ;
			int dim_corner_2 = m_linker_BG->getGeomDim<Node>(B0(0,Ny).id());
			if( dim_corner_1 == dim_corner_2 && dim_corner_1 == 2 )
			{
				for (int j=1;j<Ny;j++){
					Node n = B0(0,j);
					m_linker_BG->linkNodeToCurve(n.id(), geom_id_corner_1);
					cad::GeomCurve* curve = m_manager->getCurve(geom_id_corner_1);
					math::Point p0 = n.point();
					math::Point p = n.point();
					curve->project(p);
					n.setPoint(p);


					if (var_couche->value( B0(0,0).id() ) == 0) {
						math::Vector3d delta_p = p - p0;
						Node n_opp = B0(Nx, j);
						math::Point p_opp = n_opp.point();
						n_opp.setPoint(p_opp + delta_p);
					}

				}
			}
		}

		if ( var_couche->value( B0(0,Ny).id() ) == var_couche->value( B0(Nx,Ny).id() )) {
			int geom_id_corner_1 = m_linker_BG->getGeomId<Node>(B0(0,Ny).id()) ;
			int geom_id_corner_2 = m_linker_BG->getGeomId<Node>(B0(Nx,Ny).id()) ;
			int dim_corner_1 = m_linker_BG->getGeomDim<Node>(B0(0,Ny).id()) ;
			int dim_corner_2 = m_linker_BG->getGeomDim<Node>(B0(Nx,Ny).id());
			if( dim_corner_1 == dim_corner_2 && dim_corner_1 == 2 )
			{
				for (int i=1;i<Nx;i++){
					Node n = B0(i,Ny);
					m_linker_BG->linkNodeToCurve(n.id(), geom_id_corner_1);
					cad::GeomCurve* curve = m_manager->getCurve(geom_id_corner_1);
					math::Point p0 = n.point();
					math::Point p = n.point();
					curve->project(p);
					n.setPoint(p);


					if (var_couche->value( B0(0,0).id() ) == 0) {
						math::Vector3d delta_p = p - p0;
						Node n_opp = B0(i, 0);
						math::Point p_opp = n_opp.point();
						n_opp.setPoint(p_opp + delta_p);
					}

				}
			}
		}

		if ( var_couche->value( B0(Nx,0).id() ) == var_couche->value( B0(Nx,Ny).id() ) ) {
			int geom_id_corner_1 = m_linker_BG->getGeomId<Node>(B0(Nx,0).id()) ;
			int geom_id_corner_2 = m_linker_BG->getGeomId<Node>(B0(Nx,Ny).id()) ;
			int dim_corner_1 = m_linker_BG->getGeomDim<Node>(B0(Nx,0).id()) ;
			int dim_corner_2 = m_linker_BG->getGeomDim<Node>(B0(Nx,Ny).id());
			if( dim_corner_1 == dim_corner_2 && dim_corner_1 == 2 )
			{
				for (int j=1;j<Ny;j++){
					Node n = B0(Nx,j);
					m_linker_BG->linkNodeToCurve(n.id(), geom_id_corner_1);
					cad::GeomCurve* curve = m_manager->getCurve(geom_id_corner_1);
					math::Point p0 = n.point();
					math::Point p = n.point();
					curve->project(p);
					n.setPoint(p);


					if (var_couche->value( B0(0,0).id() ) == 0) {
						math::Vector3d delta_p = p - p0;
						Node n_opp = B0(0, j);
						math::Point p_opp = n_opp.point();
						n_opp.setPoint(p_opp + delta_p);
					}

				}
			}
		}

	}

	// New mesh of the interior points of the blocks
	for(auto B0:m_Blocking2D.allBlocks()) {
		auto Nx = B0.getNbDiscretizationI();
		auto Ny = B0.getNbDiscretizationJ();

		Array2D<TCellID> *a = new Array2D<TCellID>(Nx, Ny);
		Array2D<math::Point> pnts(Nx, Ny);

		for (auto i = 0; i < Nx; i++) {
			pnts(i, 0) = B0(i,0).point();
			pnts(i, Ny - 1) = B0(i,Ny-1).point();
		}
		for (auto j = 0; j < Ny; j++) {
			pnts(0, j) = B0(0,j).point();
			pnts(Nx - 1, j) = B0(Nx-1,j).point();
		}

		math::TransfiniteInterpolation::computeQuad(pnts);
		for (auto i = 1; i < Nx - 1; i++) {
			for (auto j = 1; j < Ny - 1; j++) {
				B0(i,j).setPoint(pnts(i,j));
			}
		}
	}

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
void
AeroPipeline_2D::ComputeVectorFieldForExtrusion(){

	Variable<math::Vector3d>* var_VectorsForExtrusion = m_meshTet->getOrCreateVariable<math::Vector3d, GMDS_NODE>("VectorField_Extrusion");

	if (m_params.vectors_field <= 0 || m_params.vectors_field > 3)
	{
		// Compute the gradient field of the level set from the wall to the external boundary
		LeastSquaresGradientComputation grad2D(m_meshTet, m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
		                                       var_VectorsForExtrusion);
		grad2D.execute();
	}
	else if (m_params.vectors_field == 1)
	{
		// Compute the gradient field of the level set from the wall to the external boundary
		LeastSquaresGradientComputation grad2D(m_meshTet, m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance"),
		                                       var_VectorsForExtrusion);
		grad2D.execute();
	}
	else if (m_params.vectors_field == 2)
	{
		Variable<math::Vector3d>* var_VectorField_1 = m_meshTet->getOrCreateVariable<math::Vector3d, GMDS_NODE>("VectorField_1");
		Variable<math::Vector3d>* var_VectorField_2 = m_meshTet->getOrCreateVariable<math::Vector3d, GMDS_NODE>("VectorField_2");

		LeastSquaresGradientComputation grad2D_1(m_meshTet, m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
		                                       var_VectorField_1);
		grad2D_1.execute();
		LeastSquaresGradientComputation grad2D_2(m_meshTet, m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance"),
		                                         var_VectorField_2);
		grad2D_2.execute();

		for (auto n_id:m_meshTet->nodes())
		{
			math::Vector3d vec_1 = var_VectorField_1->value(n_id);
			math::Vector3d vec_2 = var_VectorField_2->value(n_id);

			vec_1.normalize();
			vec_2.normalize();

			Node n = m_meshTet->get<Node>(n_id);
			math::Point p = n.point();

			if (p.X() <= m_params.x_VectorField_Z1)
			{
				var_VectorsForExtrusion->set(n_id, vec_1);
			}
			else if (p.X() >= m_params.x_VectorField_Z2)
			{
				var_VectorsForExtrusion->set(n_id, vec_2);
			}
			else
			{
				// Compute the transition field
				double alpha = (p.X() - m_params.x_VectorField_Z1)/(m_params.x_VectorField_Z2-m_params.x_VectorField_Z1) ;
				math::Vector3d v_transit = alpha*vec_2 + (1.0-alpha)*vec_1 ;
				v_transit.normalize();
				var_VectorsForExtrusion->set(n_id, v_transit);
			}

		}

		m_meshTet->deleteVariable(GMDS_NODE, var_VectorField_1);
		m_meshTet->deleteVariable(GMDS_NODE, var_VectorField_2);

	}
	else if (m_params.vectors_field == 3)
	{
		Variable<math::Vector3d>* var_VectorField_1 = m_meshTet->getOrCreateVariable<math::Vector3d, GMDS_NODE>("VectorField_1");

		LeastSquaresGradientComputation grad2D_1(m_meshTet, m_meshTet->getVariable<double,GMDS_NODE>("GMDS_Distance_Int"),
		                                         var_VectorField_1);
		grad2D_1.execute();

		for (auto n_id:m_meshTet->nodes())
		{
			math::Vector3d vec_1 = var_VectorField_1->value(n_id);
			double angle_rad = m_params.angle_attack*M_PI/180.0 ;
			math::Vector3d vec_flow({cos(angle_rad), sin(angle_rad), 0.0}) ;

			vec_1.normalize();
			vec_flow.normalize();

			Node n = m_meshTet->get<Node>(n_id);
			math::Point p = n.point();

			if (p.X() < m_params.x_VectorField_Z1)
			{
				var_VectorsForExtrusion->set(n_id, vec_1);
			}
			else if (p.X() > m_params.x_VectorField_Z2)
			{
				var_VectorsForExtrusion->set(n_id, vec_flow);
			}
			else
			{
				// Compute the transition field
				double alpha = (p.X() - m_params.x_VectorField_Z1)/(m_params.x_VectorField_Z2-m_params.x_VectorField_Z1) ;
				math::Vector3d v_transit = alpha*vec_flow + (1.0-alpha)*vec_1 ;
				v_transit.normalize();
				var_VectorsForExtrusion->set(n_id, v_transit);
			}

		}

		m_meshTet->deleteVariable(GMDS_NODE, var_VectorField_1);

	}


	// Normalisation du champ
	for (auto n_id:m_meshTet->nodes())
	{
		math::Vector3d vec = var_VectorsForExtrusion->value(n_id);
		vec.normalize();
		var_VectorsForExtrusion->set(n_id, vec);
	}

}
/*------------------------------------------------------------------------*/
