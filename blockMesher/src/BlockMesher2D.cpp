//
// Created by calderans on 05/07/22.
//

#include "gmds/blockMesher/BlockMesher2D.h"
#include "gmds/math/DiscretizationScheme1D.h"
#include "gmds/math/TransfiniteInterpolation.h"
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/

gmds::BlockMesher2D::BlockMesher2D(Mesh *ABlocks, cad::GeomMeshLinker *ALinker, cad::FACManager* AManager):
						m_blocks(ABlocks),
						m_linker(ALinker),
						m_manager(AManager){
	m_mesh = new Mesh(MeshModel(DIM2|N|E|F|N2E|N2F|E2N|E2F|F2N|F2E));
	m_linkerMesh = new cad::GeomMeshLinker;
	m_linkerMesh->setMesh(m_mesh);
	m_linkerMesh->setGeometry(m_manager);

	m_face_grids = m_blocks->newVariable<std::vector<std::vector<TCellID>>, GMDS_FACE >("grid_ref");
	m_edge_grids = m_blocks->newVariable<std::vector<TCellID>, GMDS_EDGE >("grid_ref");

}
/*----------------------------------------------------------------------------*/
gmds::BlockMesher2D::~BlockMesher2D() {
	delete m_mesh;
}
/*----------------------------------------------------------------------------*/
void BlockMesher2D::updateClassification()
{
	std::vector<Edge> boundary_edges;
	do{
		boundary_edges.clear();
		for (auto e : m_blocks->edges()) {
			Edge edge = m_blocks->get<Edge>(e);
			// Une face voisine donc on est sur le bord
			//std::cout<<"nb faces on edges "<<edge.get<Face>().size()<<std::endl;

			if (edge.get<Face>().size() == 1 && m_linker->getGeomDim<Edge>(e) == 0) {
				//std::cout<<"Test"<<std::endl;

				boundary_edges.push_back(edge);
			}
		}
		classifyBoundary(boundary_edges);
	}while(!boundary_edges.empty());
}
/*----------------------------------------------------------------------------*/
void BlockMesher2D::classifyBoundary(std::vector<Edge> &AEdges){
	for(auto &edge : AEdges){

		TCellID n0 = edge.getIDs<Node>()[0];
		TCellID n1 = edge.getIDs<Node>()[1];

		int geom_id0 = m_linker->getGeomId<Node>(n0),
			geom_id1= m_linker->getGeomId<Node>(n1);
		int geom_dim0 = m_linker->getGeomDim<Node>(n0),
			geom_dim1 = m_linker->getGeomDim<Node>(n1);

		//std::cout<<"geom dim0 "<<geom_dim0<<std::endl;
		//std::cout<<"geom dim1 "<<geom_dim1<<std::endl;

		int curve_id = -1;

		if(geom_dim0 == 1){
			if(geom_dim1 == 1){
				// Ici on est dans le cas ou l'arête du bord a ses deux sommets qui sont sur des
				// sommets géométriques donc on suppose que l'arête est sur la courbe qui relie
				// les deux sommets

				cad::GeomPoint* p0 = m_manager->getPoint(geom_id0);
				cad::GeomPoint* p1 = m_manager->getPoint(geom_id1);

				curve_id = m_manager->getCommonCurve(p0,p1);

				//std::cout<<"Cas 1"<<std::endl;

			}
			else if(geom_dim1 == 2){
				// L'arête de bloc est sur le bord donc une courbe, comme le sommet 1 est sur
				// une courbe on récupère l'id de cette courbe
				//std::cout<<"Cas 2"<<std::endl;

				curve_id = geom_id1;
			}else if(geom_dim1 == 0){
				// Ici on est dans le cas où le sommet 1 à été créé donc il n'a aucune
				// classification, il faut donc retrouver sur quelle courbe il est

				cad::GeomPoint* p0 = m_manager->getPoint(geom_id0);
				std::vector<cad::GeomCurve*> curves = p0->curves();

				float distance = 9999999;
				math::Point p_n1 = m_blocks->get<Node>(n1).point();
				int min_dist_curve = -1;


				//Supposition forte : on projette sur les courbes géom et on prend celle qui
				// est la plus proche
				for(auto c : curves){
					math::Point p_copy(p_n1);
					c->project(p_copy);
					float current_dist = p_copy.distance(p_n1);
					if(current_dist<distance){
						distance = current_dist;
						min_dist_curve = c->id();
					}
				}
				curve_id = min_dist_curve;
				m_linker->linkNodeToCurve(n1,curve_id);
				/*std::cout<<"Cas 3"<<std::endl;

				std::cout<<"Node ? "<<n0<<std::endl;
				std::cout<<"Node "<<n1<<":"<<m_linker->getGeomDim<Node>(n1)<<std::endl;
				std::cout<<"Node "<<n1<<":"<<m_linker->getGeomId<Node>(n1)<<std::endl;
				std::cout<<"curve id "<<curve_id<<std::endl;*/


			}
		}else if(geom_dim0 == 2){
			curve_id = geom_id0;
			if(geom_dim1 == 0){
				//std::cout<<"Cas 4"<<std::endl;

				m_linker->linkNodeToCurve(n1,curve_id);
				//std::cout<<"Node "<<n1<<m_linker->getGeomDim<Node>(n1)<<std::endl;

			}
		}else if(geom_dim0 == 0){

			if(geom_dim1 == 1){

				// Ici on est dans le cas où le sommet 0 à été créé donc il n'a aucune
				// classification, il faut donc retrouver sur quelle courbe il est

				cad::GeomPoint* p1 = m_manager->getPoint(geom_id1);
				std::vector<cad::GeomCurve*> curves = p1->curves();

				float distance = 9999999;
				math::Point p_n0 = m_blocks->get<Node>(n0).point();
				int min_dist_curve = -1;


				//Supposition forte : on projette sur les courbes géom et on prend celle qui
				// est la plus proche
				for(auto c : curves){
					math::Point p_copy(p_n0);
					c->project(p_copy);
					float current_dist = p_copy.distance(p_n0);
					if(current_dist<distance){
						distance = current_dist;
						min_dist_curve = c->id();
					}
				}
				//std::cout<<"Cas 5"<<std::endl;

				curve_id = min_dist_curve;
				m_linker->linkNodeToCurve(n0,curve_id);
			}else if(geom_dim1 == 2){

				// L'arête de bloc est sur le bord donc une courbe, comme le sommet 0 est sur
				// une courbe on récupère l'id de cette courbe
				curve_id = geom_id1;
				//std::cout<<"Cas 6"<<std::endl;
				m_linker->linkNodeToCurve(n0,curve_id);
				//std::cout<<"Node "<<n1<<m_linker->getGeomDim<Node>(n1)<<std::endl;

			}else if(geom_dim1 == 0){
				//std::cout<<"Cas 7"<<std::endl;

			}
		}

		if(curve_id != -1)
			m_linker->linkEdgeToCurve(edge.id(),curve_id);
	}

}
/*----------------------------------------------------------------------------*/
void BlockMesher2D::projectNodes(){
	for(auto n : m_blocks->nodes()){
		if (m_linker->getGeomDim<Node>(n) == 2){
			cad::GeomCurve* c = m_manager->getCurve(m_linker->getGeomId<Node>(n));
			Node node = m_blocks->get<Node>(n);

			math::Point projection = node.point();
			c->project(projection);
			node.setPoint(projection);
		}
	}
}
/*----------------------------------------------------------------------------*/
bool BlockMesher2D::executeMeshing(){
	meshVertices();

	meshEdges();

	meshFaces();
	return false;
}
/*----------------------------------------------------------------------------*/
void BlockMesher2D::meshVertices()
{
	for (auto v_id : m_blocks->nodes()) {
		Node v = m_blocks->get<Node>(v_id);
		Node n = m_mesh->newNode(v.point());
		m_blockV_to_meshN[v.id()] = n.id();

		//std::cout<<"Node classification "<<m_linker->getGeomDim<Node>(v_id)<<std::endl;

	}
}
/*----------------------------------------------------------------------------*/
bool BlockMesher2D::meshEdges(){
	for (auto e_id : m_blocks->edges()) {
		Edge e = m_blocks->get<Edge>(e_id);
		auto n_ids = e.getIDs<Node>();
		if (n_ids.size() != 2) return false;

		// the edge must be discritize in the gmds id increasing direction
		TCellID id_node_0 = NullID;
		TCellID id_node_1 = NullID;
		if (n_ids[0] < n_ids[1]) {
			id_node_0 = n_ids[0];
			id_node_1 = n_ids[1];
		}
		else if (n_ids[0] > n_ids[1]) {
			id_node_0 = n_ids[1];
			id_node_1 = n_ids[0];
		}
		else {
			// the edge ends are the same
			return false;
		}
		Node mesh_n0 = m_mesh->get<Node>(m_blockV_to_meshN[id_node_0]);
		Node mesh_n1 = m_mesh->get<Node>(m_blockV_to_meshN[id_node_1]);
		math::Point p0 = mesh_n0.point();
		math::Point p1 = mesh_n1.point();
		auto N = 10;
		std::vector<TCellID> discretization_ids;
		discretization_ids.resize(10 + 1);
		discretization_ids[0] = mesh_n0.id();
		double inv = 1.0 / (double) (N);
		for (auto i = 1; i < 10; i++) {
			math::Point pi = inv * ((N - i) * p0 + i * p1);
			Node ni = m_mesh->newNode(pi);
			discretization_ids[i] = ni.id();
		}
		discretization_ids[10] = mesh_n1.id();
		// we store the edge discretization
		m_blockE_to_meshN.insert(std::make_pair(e.id(),discretization_ids));


		//std::cout<<"Edge classification "<<m_linker->getGeomDim<Edge>(e.id())<<std::endl;

		if (m_linker->getGeomDim(e) == cad::GeomMeshLinker::LINK_CURVE) {
			int c_id = m_linker->getGeomId(e);
			cad::GeomCurve* c = m_manager->getCurve(c_id);
			for (auto i = 1; i < 10; i++) {
				m_linkerMesh->linkNodeToCurve(discretization_ids[i], c_id);
				math::Point project = m_mesh->get<Node>(discretization_ids[i]).point();
				c->project(project);
				m_mesh->get<Node>(discretization_ids[i]).setPoint(project);
			}
		}

		// meshFaces();
	}
return false;
}
/*----------------------------------------------------------------------------*/
bool BlockMesher2D::meshFaces()
{
	for (int f_id = 0; f_id< m_blocks->getNbFaces();f_id++) {
		Face f = m_blocks->get<Face>(f_id);
		std::vector<Edge> f_edges = f.get<Edge>();
		std::vector<Node> f_nodes = f.get<Node>();
		// Edge locally defined from 0 to 1
		Node n0 = f_nodes[0];
		Node n1 = f_nodes[1];
		bool inverted_01 = false;
		Edge e01;
		if (!getEdgeFrom(n0.id(), n1.id(), f_edges, e01, inverted_01)) return false;

		// Edge locally defined from 1 to 2
		Node n2 = f_nodes[2];
		bool inverted_12 = false;
		Edge e12;
		if (!getEdgeFrom(n1.id(), n2.id(), f_edges, e12, inverted_12)) return false;

		// Edge locally defined from 3 to 2 (to be in the same direction as e01
		Node n3 = f_nodes[3];
		bool inverted_32 = false;
		Edge e32;
		if (!getEdgeFrom(n3.id(), n2.id(), f_edges, e32, inverted_32)) return false;

		// Edge locally defined from 0 to 3 (to be in the same direction as e12
		bool inverted_03 = false;
		Edge e03;
		if (!getEdgeFrom(n0.id(), n3.id(), f_edges, e03, inverted_03)) return false;

		// now we get the 4 edges in the expected direction. We can store
		// their respective discretization in the right way into a 2D-like array
		std::vector<std::vector<TCellID>> discretization_ids;
		std::vector<std::vector<math::Point>> discretization_pnts;
		discretization_ids.resize(10 + 1);
		for (auto &v : discretization_ids) {
			v.resize(10 + 1);
		}
		discretization_pnts.resize(10 + 1);
		for (auto &v : discretization_pnts) {
			v.resize(10 + 1);
		}

		auto N = 10 + 1;
		// Corner points retrieval
		discretization_ids[0][0] = n0.id();
		discretization_ids[N - 1][0] = n1.id();
		discretization_ids[N - 1][N - 1] = n2.id();
		discretization_ids[0][N - 1] = n3.id();
		discretization_pnts[0][0] = n0.point();
		discretization_pnts[N - 1][0] = n1.point();
		discretization_pnts[N - 1][N - 1] = n2.point();
		discretization_pnts[0][N - 1] = n3.point();

		// Boundary edges point retrieval for edge 01
		std::vector<TCellID> edge_disc = m_blockE_to_meshN[e01.id()];

		if (n0.id() != edge_disc[0]) std::reverse(edge_disc.begin(), edge_disc.end());
		for (auto i = 1; i < N - 1; i++) {
			discretization_ids[i][0] = edge_disc[i];
			discretization_pnts[i][0] = m_mesh->get<Node>(edge_disc[i]).point();
		}
		// Boundary edges point retrieval for edge 01
		edge_disc = m_blockE_to_meshN[e12.id()];
		if (n1.id() != edge_disc[0]) std::reverse(edge_disc.begin(), edge_disc.end());

		for (auto i = 1; i < N - 1; i++) {
			discretization_ids[N - 1][i] = edge_disc[i];
			discretization_pnts[N - 1][i] = m_mesh->get<Node>(edge_disc[i]).point();
		}
		// Boundary edges point retrieval for edge 32
		edge_disc = m_blockE_to_meshN[e32.id()];
		if (n3.id() != edge_disc[0]) std::reverse(edge_disc.begin(), edge_disc.end());
		for (auto i = 1; i < N - 1; i++) {
			discretization_ids[i][N - 1] = edge_disc[i];
			discretization_pnts[i][N - 1] = m_mesh->get<Node>(edge_disc[i]).point();
		}
		// Boundary edges point retrieval for edge 01
		edge_disc = m_blockE_to_meshN[e03.id()];
		if (n0.id() != edge_disc[0]) std::reverse(edge_disc.begin(), edge_disc.end());
		for (auto i = 1; i < N - 1; i++) {
			discretization_ids[0][i] = edge_disc[i];
			discretization_pnts[0][i] = m_mesh->get<Node>(edge_disc[i]).point();
		}
		// compute the points location and create node
		math::TransfiniteInterpolation::compute(discretization_pnts);
		for(auto i=1; i<N-1;i++) {
			for(auto j=1; j<N-1;j++) {
				Node nij = m_mesh->newNode(discretization_pnts[i][j]);
				discretization_ids[i][j]= nij.id();
			}
		}
		m_blockF_to_meshN[f_id] = discretization_ids;
	}
	//std::cout<<"NB Faces "<<m_blocks->getNbFaces()<<std::endl;

	Variable<int>* block = m_mesh->newVariable<int,GMDS_FACE>("blocks");

	for(TCellID f_id = 0; f_id < m_blocks->getNbFaces();f_id++){
		//std::cout<<"Face id "<<f_id<<std::endl;

		std::vector<TCellID> face_ids;
		std::vector<std::vector<TCellID> > n_ids = m_blockF_to_meshN[f_id];
		auto I = n_ids.size();
		auto J = n_ids[0].size();
		for(auto i=0;i<I-1;i++){
			for(auto j=0;j<J-1;j++){
				Face q = m_mesh->newQuad(n_ids[i][j],n_ids[i+1][j],n_ids[i+1][j+1],n_ids[i][j+1]);
				face_ids.push_back(q.id());
				block->set(q.id(), f_id);
			}
		}
	}
	return false;
}
/*----------------------------------------------------------------------------*/
bool BlockMesher2D::getEdgeFrom(const TCellID AN0, const TCellID AN1,
               						const std::vector<Edge> &AEdges, Edge &AResult,
               						bool &AInverted){
	for(auto e: AEdges){
		std::vector<TCellID > e_node_ids = e.getIDs<Node>();
		if(e_node_ids[0]==AN0 && e_node_ids[1]==AN1){
			AResult=e;
			AInverted=false;
			return true;
		} else if(e_node_ids[0]==AN1 && e_node_ids[1]==AN0){
			AResult=e;
			AInverted=true;
			return true;
		}
	}
	return false;
}
/*----------------------------------------------------------------------------*/
Mesh* BlockMesher2D::getMesh(){
	std::cout<<"Nb nodes quad "<<m_mesh->getNbNodes()<<std::endl;
	std::cout<<"Nb edges quad "<<m_mesh->getNbEdges()<<std::endl;
	std::cout<<"Nb faces quad "<<m_mesh->getNbFaces()<<std::endl;

	return m_mesh;
}