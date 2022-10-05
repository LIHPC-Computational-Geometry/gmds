/*----------------------------------------------------------------------------*/
#include <gmds/ig/Blocking2D.h>
#include <gmds/math/TransfiniteInterpolation.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
Blocking2D::Blocking2D()
  : Mesh(MeshModel(DIM3|F|E|N|F2N|E2N|F2E|E2F|N2E|N2F))
{
	m_embedding_dim = newVariable<int, GMDS_NODE>("embedding_dim");
	m_embedding_id = newVariable<TCellID , GMDS_NODE>("embedding_id");
	m_discretization_I= newVariable<int, GMDS_FACE>("discretizationI");
	m_discretization_J= newVariable<int, GMDS_FACE>("discretizationJ");
	m_face_grids = newVariable<Array2D<TCellID>*, GMDS_FACE >("grid_ref");
	m_edge_grids = newVariable<std::vector<TCellID>*, GMDS_EDGE >("grid_ref");

}
/*----------------------------------------------------------------------------*/
Blocking2D::Blocking2D(const Mesh& AMesh): Mesh(AMesh.getModel()){

	// for now, we assume that the mesh hasn't already been used for a blocking

	m_embedding_dim = newVariable<int, GMDS_NODE>("embedding_dim");
	m_embedding_id = newVariable<TCellID , GMDS_NODE>("embedding_id");
	m_discretization_I= newVariable<int, GMDS_FACE>("discretizationI");
	m_discretization_J= newVariable<int, GMDS_FACE>("discretizationJ");
	m_face_grids = newVariable<Array2D<TCellID>*, GMDS_FACE >("grid_ref");
	m_edge_grids = newVariable<std::vector<TCellID>*, GMDS_EDGE >("grid_ref");

	//Building the blocking from the mesh
	// Maybe make an independent function for the building, with different types of discretization method etc..

	//map the relation from mesh vertices to blocking vertices
	std::map<TCellID,TCellID> n2n;
	for (auto n : AMesh.nodes()) {
		Node n_b = newBlockCorner(AMesh.get<Node>(n).point());
		n2n[n] = n_b.id();
	}
	for (auto f : AMesh.faces()) {
		Face face = AMesh.get<Face>(f);
		std::vector<TCellID> f_nodes;
		face.getIDs<Node>(f_nodes);
		Block b = newBlock(n2n[f_nodes[0]],n2n[f_nodes[1]],n2n[f_nodes[2]],n2n[f_nodes[3]]);
	}



}
/*----------------------------------------------------------------------------*/
Blocking2D::~Blocking2D() {
	//by deleting the mesh, we also delete the memory space pointed by m_embedding
	// but we need to delete explicitly allocated data
	for(auto f_id:faces()){
		delete m_face_grids->value(f_id);
		m_face_grids->set(f_id,NULL);
	}
	for(auto e_id:edges()){
		delete m_edge_grids->value(e_id);
		m_edge_grids->set(e_id,NULL);
	}

}
/*----------------------------------------------------------------------------*/
Node Blocking2D::newBlockCorner(const math::Point &APnt) {
	Node n = newNode(APnt);
	m_embedding_dim->set(n.id(),0);
	m_embedding_id->set(n.id(),n.id());
	return n;
}
/*----------------------------------------------------------------------------*/
Node Blocking2D::newBlockCorner(const TCoord AX,const TCoord AY, const TCoord AZ) {
	Node n = newNode(math::Point(AX,AY,AZ));
	m_embedding_dim->set(n.id(),0);
	m_embedding_id->set(n.id(),n.id());
	return n;
}
/*----------------------------------------------------------------------------*/
Blocking2D::Block Blocking2D::newBlock(const TCellID AN1, const TCellID AN2,
                     const TCellID AN3, const TCellID AN4){

	Face b = newQuad(AN1,AN2,AN3,AN4);
	Edge e1 = get<Edge>(getEdge(AN1,AN2));
	b.add(e1); e1.add(b);
	Edge e2 = get<Edge>(getEdge(AN2,AN3));
	b.add(e2); e2.add(b);
	Edge e3 = get<Edge>(getEdge(AN3,AN4));
	b.add(e3); e3.add(b);
	Edge e4 = get<Edge>(getEdge(AN4,AN1));
	b.add(e4); e4.add(b);

	Node n1 = get<Node>(AN1);
	n1.add(b);
	Node n2 = get<Node>(AN2);
	n2.add(b);
	Node n3 = get<Node>(AN3);
	n1.add(b);
	Node n4 = get<Node>(AN4);
	n4.add(b);

	return Block(b, this);
}
/*----------------------------------------------------------------------------*/
Blocking2D::Block Blocking2D::newBlock(const Node &AN1, const Node &AN2,
                     const Node &AN3, const Node &AN4){
	return newBlock(AN1.id(),AN2.id(),AN3.id(),AN4.id());
}
/*----------------------------------------------------------------------------*/
Blocking2D::Block  Blocking2D::
   block(const TCellID AId) {
	return Block(get<Face>(AId),this);
}
/*----------------------------------------------------------------------------*/
std::vector<Blocking2D::Block> Blocking2D::allBlocks() {
	std::vector<Blocking2D::Block> blocks;
	for(auto f_id:faces()){
		blocks.push_back(block(f_id));
	}
	return blocks;
}
/*----------------------------------------------------------------------------*/
TCellID Blocking2D::getEdge(const TCellID AN1, const TCellID AN2) {
	std::vector<TCellID> edges_1 = get<Node>(AN1).getIDs<Edge>();
	std::vector<TCellID> edges_2 = get<Node>(AN2).getIDs<Edge>();

	for(auto e1:edges_1){
		for(auto e2:edges_2){
			if(e1==e2){
				//found common edge
				return e1;
			}
		}
	}
	//there is not a common edge, so we add it

	Edge e = newEdge(AN1,AN2);
	//E2N update
	get<Node>(AN1).add(e);
	get<Node>(AN2).add(e);
	return e.id();
}
/*----------------------------------------------------------------------------*/
bool Blocking2D::checkDiscretizationValidity() const {
	for(auto e_id:edges()){
		Edge ei = get<Edge>(e_id);
		std::vector<Face> adj_faces = ei.get<Face>();
		if(adj_faces.size()>1){
			//we verify that we have the same discretization scheme on both sides
			if(adj_faces.size()!=2) {
				throw GMDSException("An edge should be adjacent to 1 or 2 faces");
			}
			int nb_subdiv0 = getNbDiscretization(adj_faces[0],ei);
			int nb_subdiv1 = getNbDiscretization(adj_faces[1],ei);
			if(nb_subdiv0!=nb_subdiv1){
				return false;
			}
		}
	}

	return true;
}
/*----------------------------------------------------------------------------*/
void Blocking2D::initializeEdgesPoints()
{
	// first we check that edges are discretiezd similarly by blocks sharing
	// them
	checkDiscretizationValidity();
	// then we mesh edges
	for (auto e_id : edges()) {
		Edge ei = get<Edge>(e_id);
		int nb_subdiv = getNbDiscretization(ei.get<Face>()[0], ei);
		// We use a discrete linear interpolation along the edge
		std::vector<Node> ei_end_nodes = ei.get<Node>();
		Node n0 = ei_end_nodes[0];
		Node n1 = ei_end_nodes[1];
		math::DiscretizationScheme1DUniform d(n0.point(), n1.point(), nb_subdiv);
		std::vector<TCellID> *edge_disc = new std::vector<TCellID>();
		edge_disc->push_back(n0.id());
		for (auto i = 1; i < nb_subdiv - 1; i++) {
			Node ni = newNode(d(i));
			m_embedding_dim->set(ni.id(), 1);
			m_embedding_id->set(ni.id(), e_id);
			edge_disc->push_back(ni.id());
		}
		edge_disc->push_back(n1.id());
		m_edge_grids->set(e_id, edge_disc);
	}
}
/*----------------------------------------------------------------------------*/
void Blocking2D::initializeBlocksPoints(){
	for(auto f_id:faces()){
		Block bi(f_id,this);
		auto nb_I = bi.getNbDiscretizationI();
		auto nb_J = bi.getNbDiscretizationJ();

		Array2D<TCellID>* a = new Array2D<TCellID>(nb_I,nb_J);
		Array2D<math::Point> pnts(nb_I,nb_J);

		Node n0 = bi.getNode(0);
		Node n1 = bi.getNode(1);
		Node n2 = bi.getNode(2);
		Node n3 = bi.getNode(3);

		(*a)(0,0)           =n0.id();
		(*a)(nb_I-1,0)      =n1.id();
		(*a)(nb_I-1,nb_J-1) =n2.id();
		(*a)(0,nb_J-1)      =n3.id();

		pnts(0,0)           =n0.point();
		pnts(nb_I-1,0)      =n1.point();
		pnts(nb_I-1,nb_J-1) =n2.point();
		pnts(0,nb_J-1)      =n3.point();

		Edge e01 = bi.getEdge(0,1);
		std::vector<TCellID>* e01_nodes = m_edge_grids->value(e01.id());
		if((*e01_nodes)[0]==n1.id()){
			//need to reverse
			std::reverse(e01_nodes->begin(),e01_nodes->end());
		}

		Edge e12 = bi.getEdge(1,2);
		std::vector<TCellID>* e12_nodes = m_edge_grids->value(e12.id());
		if((*e12_nodes)[0]==n2.id()){
			//need to reverse
			std::reverse(e12_nodes->begin(),e12_nodes->end());
		}

		Edge e32 = bi.getEdge(3,2);
		std::vector<TCellID>* e32_nodes = m_edge_grids->value(e32.id());
		if((*e32_nodes)[0]==n2.id()){
			//need to reverse
			std::reverse(e32_nodes->begin(),e32_nodes->end());
		}

		Edge e03 = bi.getEdge(0,3);
		std::vector<TCellID>* e03_nodes = m_edge_grids->value(e03.id());
		if((*e03_nodes)[0]==n3.id()){
			//need to reverse
			std::reverse(e03_nodes->begin(),e03_nodes->end());
		}

		for(auto i=1; i<e01_nodes->size()-1;i++){
			(*a)(i,0) =(*e01_nodes)[i];
			pnts(i,0)= get<Node>((*e01_nodes)[i]).point();
		}
		for(auto i=1; i<e32_nodes->size()-1;i++){
			(*a)(i,nb_J-1) =(*e32_nodes)[i];
			pnts(i,nb_J-1)= get<Node>((*e32_nodes)[i]).point();
		}
		for(auto i=1; i<e03_nodes->size()-1;i++){
			(*a)(0,i) =(*e03_nodes)[i];
			pnts(0,i)= get<Node>((*e03_nodes)[i]).point();
		}
		for(auto i=1; i<e12_nodes->size()-1;i++){
			(*a)(nb_I-1,i) =(*e12_nodes)[i];
			pnts(nb_I-1,i)= get<Node>((*e12_nodes)[i]).point();
		}

		math::TransfiniteInterpolation::computeQuad(pnts);
		for(auto i=1; i<nb_I-1;i++){
			for(auto j=1; j<nb_J-1;j++) {
				Node nij = newNode(pnts(i,j));
				(*a)(i,j)=nij.id();
				m_embedding_dim->set(nij.id(),2);
				m_embedding_id->set(nij.id(),f_id);
			}
		}
		m_face_grids->set(f_id,a);
	}
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID>* Blocking2D::getEdgeGrid(TCellID AID){
	return m_edge_grids->value(AID);
}

/*----------------------------------------------------------------------------*/
void Blocking2D::initializeGridPoints() {
	// first we check that edges are discretiezd similarly by blocks sharing
	// them
	checkDiscretizationValidity();
	//then we mesh edges
	for(auto e_id:edges()){
		Edge ei = get<Edge>(e_id);
		int nb_subdiv = getNbDiscretization(ei.get<Face>()[0],ei);
		//We use a discrete linear interpolation along the edge
		std::vector<Node> ei_end_nodes = ei.get<Node>();
		Node n0 = ei_end_nodes[0];
		Node n1 = ei_end_nodes[1];
		math::DiscretizationScheme1DUniform d(n0.point(),n1.point(),nb_subdiv);
		std::vector<TCellID>* edge_disc = new std::vector<TCellID>();
		edge_disc->push_back(n0.id());
		for(auto i=1; i<nb_subdiv-1;i++){
			Node ni = newNode(d(i));
			m_embedding_dim->set(ni.id(),1);
			m_embedding_id->set(ni.id(),e_id);
			edge_disc->push_back(ni.id());
		}
		edge_disc->push_back(n1.id());
		m_edge_grids->set(e_id,edge_disc);
	}

	//and now faces. We simply use a transfinite interpolation right now
	for(auto f_id:faces()){
		Block bi(f_id,this);
		auto nb_I = bi.getNbDiscretizationI();
		auto nb_J = bi.getNbDiscretizationJ();

		Array2D<TCellID>* a = new Array2D<TCellID>(nb_I,nb_J);
		Array2D<math::Point> pnts(nb_I,nb_J);

		Node n0 = bi.getNode(0);
		Node n1 = bi.getNode(1);
		Node n2 = bi.getNode(2);
		Node n3 = bi.getNode(3);

		(*a)(0,0)           =n0.id();
		(*a)(nb_I-1,0)      =n1.id();
		(*a)(nb_I-1,nb_J-1) =n2.id();
		(*a)(0,nb_J-1)      =n3.id();

		pnts(0,0)           =n0.point();
		pnts(nb_I-1,0)      =n1.point();
		pnts(nb_I-1,nb_J-1) =n2.point();
		pnts(0,nb_J-1)      =n3.point();

		Edge e01 = bi.getEdge(0,1);
		std::vector<TCellID>* e01_nodes = m_edge_grids->value(e01.id());
		if((*e01_nodes)[0]==n1.id()){
			//need to reverse
			std::reverse(e01_nodes->begin(),e01_nodes->end());
		}

		Edge e12 = bi.getEdge(1,2);
		std::vector<TCellID>* e12_nodes = m_edge_grids->value(e12.id());
		if((*e12_nodes)[0]==n2.id()){
			//need to reverse
			std::reverse(e12_nodes->begin(),e12_nodes->end());
		}

		Edge e32 = bi.getEdge(3,2);
		std::vector<TCellID>* e32_nodes = m_edge_grids->value(e32.id());
		if((*e32_nodes)[0]==n2.id()){
			//need to reverse
			std::reverse(e32_nodes->begin(),e32_nodes->end());
		}

		Edge e03 = bi.getEdge(0,3);
		std::vector<TCellID>* e03_nodes = m_edge_grids->value(e03.id());
		if((*e03_nodes)[0]==n3.id()){
			//need to reverse
			std::reverse(e03_nodes->begin(),e03_nodes->end());
		}

		for(auto i=1; i<e01_nodes->size()-1;i++){
			(*a)(i,0) =(*e01_nodes)[i];
			pnts(i,0)= get<Node>((*e01_nodes)[i]).point();
		}
		for(auto i=1; i<e32_nodes->size()-1;i++){
			(*a)(i,nb_J-1) =(*e32_nodes)[i];
			pnts(i,nb_J-1)= get<Node>((*e32_nodes)[i]).point();
		}
		for(auto i=1; i<e03_nodes->size()-1;i++){
			(*a)(0,i) =(*e03_nodes)[i];
			pnts(0,i)= get<Node>((*e03_nodes)[i]).point();
		}
		for(auto i=1; i<e12_nodes->size()-1;i++){
			(*a)(nb_I-1,i) =(*e12_nodes)[i];
			pnts(nb_I-1,i)= get<Node>((*e12_nodes)[i]).point();
		}

		math::TransfiniteInterpolation::computeQuad(pnts);
		for(auto i=1; i<nb_I-1;i++){
			for(auto j=1; j<nb_J-1;j++) {
				Node nij = newNode(pnts(i,j));
				(*a)(i,j)=nij.id();
				m_embedding_dim->set(nij.id(),2);
				m_embedding_id->set(nij.id(),f_id);
			}
		}
		m_face_grids->set(f_id,a);
	}
}
/*----------------------------------------------------------------------------*/
int Blocking2D::getNbDiscretization(const Face& AFace, const Edge &AEdge) const {
	//find the edge position in the face
	std::vector<TCellID> face_nids = AFace.getIDs<Node>();
	std::vector<TCellID> edge_nids = AEdge.getIDs<Node>();
	auto edge_position = -1;
	for(auto i=0;i<4 && edge_position==-1;i++){
		TCellID f_ni = face_nids[i];
		TCellID f_nj = face_nids[(i+1)%4];
		if((f_ni==edge_nids[0] && f_nj==edge_nids[1]) ||
		    (f_ni==edge_nids[1] && f_nj==edge_nids[0])){
			edge_position=i;
		}
	}
	//found the edge position
	if(edge_position==-1){
		throw GMDSException("Error: edge not adjacent to the face");
	}
	if(edge_position%2==0){
		//edge 0 or 2, so subdivsion I
		return m_discretization_I->value(AFace.id());
	}
	//otherwise, edge 1 or 3 so subdivision J
	return m_discretization_J->value(AFace.id());
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> Blocking2D::getNodeNeighbors(TCellID AId){

	std::vector<TCellID> neighbors;

	int dim = m_embedding_dim->value(AId);
	TCellID id = m_embedding_id->value(AId);

	if(dim == 0){

		Node n = get<Node>(AId);

		std::cout<<"Nb faces "<<n.get<Face>().size()<<std::endl;

		if(n.nbFaces()%2 != 0){
			throw GMDSException("Node is singular, cannot get neighboring nodes.");
		}

		for(auto const &f : n.get<Face>()){
			int n_position = -1;
			std::vector<TCellID> f_nids;
			f.getIDs<Node>(f_nids);

			for(int i = 0; i < 4 && n_position == -1; i++){
				if(f_nids[i] == n.id()){
					n_position = i;
				}
			}

			Block b = block(f.id());

			switch (n_position) {
			case 0:
				neighbors.push_back(b(1,1).id());
				break;
			case 1:
				neighbors.push_back(b(b.getNbDiscretizationI()-2,1).id());
				break;
			case 2:
				neighbors.push_back(b(b.getNbDiscretizationI()-2,b.getNbDiscretizationJ()-2).id());
				break;
			case 3:
				neighbors.push_back(b(1,b.getNbDiscretizationJ()-2).id());
				break;
			default:
				throw GMDSException("Error: node not adjacent to the face");
			}
		}

		for(auto const &e : n.get<Edge>()){

			std::vector<TCellID>* e_nids = m_edge_grids->value(e.id());

			if((*e_nids)[0] == n.id()){
				neighbors.push_back((*e_nids)[1]);
			}
			else if((*e_nids)[e_nids->size()-1] == n.id()){
				neighbors.push_back((*e_nids)[e_nids->size()-2]);
			}
			else{
				throw GMDSException("Error: node not adjacent to the edge");
			}
		}
	}else if(dim == 1) {

		Edge e = get<Edge>(id);

		std::vector<TCellID> *edge_grid = m_edge_grids->value(id);

		// used to store the starting point to get the neighborhood in the adjacent faces
		int indice_on_edge;

		for (int i = 1; i < edge_grid->size(); i++) {
			if ((*edge_grid)[i] == AId) {
				neighbors.push_back((*edge_grid)[i - 1]);
				neighbors.push_back((*edge_grid)[i + 1]);
				indice_on_edge = i;
			}
		}

		for (auto const &f : e.get<Face>()) {

			std::vector<TCellID> face_eids;
			f.getIDs<Edge>(face_eids);

			auto edge_position = -1;
			for (auto i = 0; i < 4 && edge_position == -1; i++) {
				if (face_eids[i] == id) {
					edge_position = i;
				}
			}

			Block b = Block(f, this);

			switch (edge_position) {
			case 0:     // Case edge (i,0)
				neighbors.push_back(b(indice_on_edge - 1, 1).id());
				neighbors.push_back(b(indice_on_edge, 1).id());
				neighbors.push_back(b(indice_on_edge + 1, 1).id());
				break;
			case 1:     // Case edge (I,j)
				neighbors.push_back(b(b.getNbDiscretizationI() - 2, indice_on_edge - 1).id());
				neighbors.push_back(b(b.getNbDiscretizationI() - 2, indice_on_edge).id());
				neighbors.push_back(b(b.getNbDiscretizationI() - 2, indice_on_edge + 1).id());
				break;
			case 2:     // Case edge (i,J)
				neighbors.push_back(b(indice_on_edge - 1, b.getNbDiscretizationJ() - 2).id());
				neighbors.push_back(b(indice_on_edge, b.getNbDiscretizationJ() - 2).id());
				neighbors.push_back(b(indice_on_edge + 1, b.getNbDiscretizationJ() - 2).id());
				break;
			case 3:     // Case edge (0,j)
				neighbors.push_back(b(1, indice_on_edge - 1).id());
				neighbors.push_back(b(1, indice_on_edge).id());
				neighbors.push_back(b(1, indice_on_edge + 1).id());
				break;
			default:
				throw GMDSException("Error: edge not adjacent to the face");
			}
		}

	}else if(dim == 2){
		Block b = block(id);
		int nbDiscrI = b.getNbDiscretizationI();
		int nbDiscrJ = b.getNbDiscretizationJ();

		//the interval is [1, nbDiscr-1] because we know the node is on a face and not an edge
		bool found = false;
		for(int i = 1; i<nbDiscrI-1 && !found; i++){
			for(int j = 1; j<nbDiscrJ-1 && !found; j++){

				if(b(i,j).id() == AId){

					found = true;

					for(int x = -1; x<=1; x++){
						for(int y = -1; y<=1; y++){
							if(x == 0 && y == 0) continue;

							neighbors.push_back(b(i+x,j+y).id());
						}
					}
				}
			}
		}

	}else{
		throw GMDSException("Wrong node dimension for 2D");
	}

	return neighbors;
}
/*----------------------------------------------------------------------------*/
void Blocking2D::buildBlocks(const int AN){

	for (auto n : nodes()) {

		m_embedding_dim->set(n,0);
		m_embedding_id->set(n,n);
	}
	for (auto f : faces()) {

		std::vector<TCellID> nodes;
		Face face = get<Face>(f);
		face.getIDs<Node>(nodes);

		Edge e1 = get<Edge>(getEdge(nodes[0],nodes[1]));
		face.add(e1); e1.add(face);
		Edge e2 = get<Edge>(getEdge(nodes[1],nodes[2]));
		face.add(e2); e2.add(face);
		Edge e3 = get<Edge>(getEdge(nodes[2],nodes[3]));
		face.add(e3); e3.add(face);
		Edge e4 = get<Edge>(getEdge(nodes[3],nodes[0]));
		face.add(e4); e4.add(face);

		Node n1 = get<Node>(nodes[0]);
		n1.add(face);
		Node n2 = get<Node>(nodes[1]);
		n2.add(face);
		Node n3 = get<Node>(nodes[2]);
		n3.add(face);
		Node n4 = get<Node>(nodes[3]);
		n4.add(face);

		Block b(f, this);
		b.setNbDiscretizationI(AN);
		b.setNbDiscretizationJ(AN);
	}
}
/*----------------------------------------------------------------------------*/
TCellID Blocking2D::getBlockingId(TCellID AId){
	return m_embedding_id->value(AId);
}
/*----------------------------------------------------------------------------*/
int Blocking2D::getBlockingDim(TCellID AId){
	return m_embedding_dim->value(AId);
}
/*----------------------------------------------------------------------------*/

/*============================================================================*/
//Methods for the friend class Block
/*============================================================================*/
Blocking2D::Block::Block(const Face &AFace,  Blocking2D* ASupport)
  : m_face(AFace), m_support(ASupport), m_grid_view(ASupport->m_face_grids->value(AFace.id()))
{}
/*----------------------------------------------------------------------------*/
Blocking2D::Block::Block(const TCellID AFaceID, Blocking2D *ASupport) {
	m_support= ASupport;
	m_face = m_support->get<Face>(AFaceID);
	m_grid_view=m_support->m_face_grids->value(AFaceID);
}
/*----------------------------------------------------------------------------*/
TCellID Blocking2D::Block::origin() {return m_face.getIDs<Node>()[0];}
/*----------------------------------------------------------------------------*/
Node Blocking2D::Block::getNode(const int &AIndex) {
	if(AIndex<0 || AIndex>3){
		throw GMDSException("A 2D block has exactly 4 nodes in range [0,3]");
	}
	return m_face.get<Node>()[AIndex];
}
/*----------------------------------------------------------------------------*/
Node Blocking2D::Block::operator()(const int AI, const int AJ) {
	/**
	 *   J = 0    e2
	 *       ____________
	 *      |			  |
	 *  e3  |			  |  e1
	 *      |			  |
	 *      |___________|  I=0
	 *  Origin   e0
	 */

	bool minus_I = AI < 0;
	bool minus_J = AJ < 0;
	bool major_I = AI > this->getNbDiscretizationI()-1;
	bool major_J = AJ > this->getNbDiscretizationJ()-1;

	int outI = 0;
	int outJ = 0;

	int resultI;
	int resultJ;

	if(minus_I){

		outI = -(AI);

		Edge edge = this->getEdge(0,3);
		std::vector<TCellID> fids = edge.getIDs<Face>();
		if(fids.size() == 2){
			TCellID face_id = fids[0] == this->id() ? fids[1] : fids[0];
			Block b = m_support->block(face_id);

			int b_edge_index = -1;

			for(int i = 0; i<4; i++){
				if(b.getEdge(i,(i+1)%3).id() == edge.id()){
					b_edge_index = i;
					break;
				}
			}

			if(b_edge_index == 0){
				resultI = AJ;
				resultJ = outI;

			}else if(b_edge_index == 1){
				resultI = (b.getNbDiscretizationI()-1)-outI;
				resultJ = AJ;

			}else if(b_edge_index == 2){
				resultI = AJ;
				resultJ = (b.getNbDiscretizationJ()-1)-outI;

			}else if(b_edge_index == 3){
				resultI = outI;
				resultJ = AJ;

			}else{
				throw GMDSException("Face doesn't have the edge.");
			}

			return b(resultI,resultJ);

		}

	}else if(major_I){

		outI = AI - (getNbDiscretizationI()-1);

		Edge edge = this->getEdge(1,2);
		std::vector<TCellID> fids = edge.getIDs<Face>();
		if(fids.size() == 2){
			TCellID face_id = fids[0] == this->id() ? fids[1] : fids[0];
			Block b = m_support->block(face_id);

			int b_edge_index = -1;


			for(int i = 0; i<4; i++){


				int modulo = (i+1)%4;

				if(b.getEdge(i,(i+1)%4).id() == edge.id()){
					b_edge_index = i;
					break;
				}
			}

			if(b_edge_index == 0){
				resultI = AJ;
				resultJ = outI;

			}else if(b_edge_index == 1){
				resultI = (b.getNbDiscretizationI()-1)-outI;
				resultJ = AJ;

			}else if(b_edge_index == 2){
				resultI = AJ;
				resultJ = (b.getNbDiscretizationJ()-1)-outI;

			}else if(b_edge_index == 3){
				resultI = outI;
				resultJ = AJ;

			}else{
				throw GMDSException("Face doesn't have the edge.");
			}

			return b(resultI,resultJ);

		}
	}else{
		if(minus_J){
			Edge edge = this->getEdge(0,1);

			outJ = -(AI);

			std::vector<TCellID> fids = edge.getIDs<Face>();
			if(fids.size() == 2){
				TCellID face_id = fids[0] == this->id() ? fids[1] : fids[0];
				Block b = m_support->block(face_id);

				int b_edge_index = -1;

				for(int i = 0; i<4; i++){
					if(b.getEdge(i,(i+1)%3).id() == edge.id()){
						b_edge_index = i;
						break;
					}
				}

				if(b_edge_index == 0){
					resultI = AI;
					resultJ = outJ;

				}else if(b_edge_index == 1){
					resultI = (b.getNbDiscretizationI()-1)-outJ;
					resultJ = AI;

				}else if(b_edge_index == 2){
					resultI = AI;
					resultJ = (b.getNbDiscretizationJ()-1)-outJ;

				}else if(b_edge_index == 3){
					resultI = outJ;
					resultJ = AI;

				}else{
					throw GMDSException("Face doesn't have the edge.");
				}

				return b(resultI,resultJ);

			}
		}else if (major_J){
			Edge edge = this->getEdge(2,3);

			outJ = AJ - (getNbDiscretizationJ()-1);

			std::vector<TCellID> fids = edge.getIDs<Face>();
			if(fids.size() == 2){
				TCellID face_id = fids[0] == this->id() ? fids[1] : fids[0];
				Block b = m_support->block(face_id);

				int b_edge_index = -1;

				for(int i = 0; i<4; i++){
					if(b.getEdge(i,(i+1)%3).id() == edge.id()){
						b_edge_index = i;
						break;
					}
				}

				if(b_edge_index == 0){
					resultI = AI;
					resultJ = outI;

				}else if(b_edge_index == 1){
					resultI = (b.getNbDiscretizationI()-1)-outJ;
					resultJ = AI;

				}else if(b_edge_index == 2){
					resultI = AI;
					resultJ = (b.getNbDiscretizationJ()-1)-outJ;

				}else if(b_edge_index == 3){
					resultI = outJ;
					resultJ = AI;

				}else{
					throw GMDSException("Face doesn't have the edge.");
				}

				return b(resultI,resultJ);

			}
		}else{
			return m_support->get<Node>((*m_grid_view)(AI,AJ));
		}
	}
}
/*----------------------------------------------------------------------------*/
Edge Blocking2D::Block::getEdgeI() {
	return getEdge(0,1);
}
/*----------------------------------------------------------------------------*/
Edge Blocking2D::Block::getEdgeJ() {
	return getEdge(0,3);
}
/*----------------------------------------------------------------------------*/
Edge Blocking2D::Block::getEdge(const int AI, const int AJ)  {
	if(AI==AJ)
		throw GMDSException("Error in local numbering for getting an edge: same node value");
	if(AI<0 || AJ<0)
		throw GMDSException("Error in local numbering for getting an edge: negative local numbering");
	if(AI>3 || AJ>3)
		throw GMDSException("Error in local numbering for getting an edge: a local number can not exceed 3");
	std::vector<TCellID> nids = m_face.getIDs<Node>();
	std::vector<TCellID> edges_1 = m_support->get<Node>(nids[AI]).getIDs<Edge>();
	std::vector<TCellID> edges_2 = m_support->get<Node>(nids[AJ]).getIDs<Edge>();

	for(auto e1:edges_1){
		for(auto e2:edges_2){
			if(e1==e2){
				//found common edge
				return m_support->get<Edge>(e1);
			}
		}
	}
	throw GMDSException("error no edge found!");
}
/*----------------------------------------------------------------------------*/
void Blocking2D::Block::setNbDiscretizationI(const int AN) {
	if(AN<=0)
		throw GMDSException("discretization value must be strictly positive");

	m_support->m_discretization_I->set(m_face.id(),AN);
}
/*----------------------------------------------------------------------------*/
void Blocking2D::Block::setNbDiscretizationJ(const int AN) {
	if(AN<=0)
		throw GMDSException("discretization value must be strictly positive");

	m_support->m_discretization_J->set(m_face.id(),AN);
}
/*----------------------------------------------------------------------------*/
int Blocking2D::Block::getNbDiscretizationI() const{
	return m_support->m_discretization_I->value(m_face.id());
}
/*----------------------------------------------------------------------------*/
int Blocking2D::Block::getNbDiscretizationJ() const{
	return m_support->m_discretization_J->value(m_face.id());
}
/*----------------------------------------------------------------------------*/
math::Vector3d Blocking2D::Block::getUnitVectorI() {
	std::vector<Node> nids = m_face.get<Node>();
	math::Vector3d v=nids[1].point()-nids[0].point();
	return v.getNormalize();
}
/*----------------------------------------------------------------------------*/
math::Vector3d Blocking2D::Block::getUnitVectorJ() {
	std::vector<Node> nids = m_face.get<Node>();
	math::Vector3d v=nids[3].point()-nids[0].point();
	return v.getNormalize();
}
/*----------------------------------------------------------------------------*/
std::pair<int,int> Blocking2D::Block::getIndices(const TCellID AID){
	//Maybe kind of heavy method
	//Should find a smarter way to find the point
	for(int i = 0; i<getNbDiscretizationI(); i++){
		for(int j = 0; j<getNbDiscretizationI(); j++){
			if((*m_grid_view)(i,j) == AID){
				return {i,j};
			}
		}
	}
	throw GMDSException("ERROR: Node not in the block grid.");
}
/*----------------------------------------------------------------------------*/
bool Blocking2D::Block::isEdgeOnI(const TCellID AID) {
	std::vector<TCellID> e_nodes = m_support->get<Edge>(AID).getIDs<Node>();
	std::vector<TCellID> nids = m_face.getIDs<Node>();
	//Same as getEdge, we use the node indices (more suitable for class changes?)
	return (e_nodes[0] == nids[0] && e_nodes[1] == nids[1]) || (e_nodes[0] == nids[1] && e_nodes[1] == nids[0])
	       || (e_nodes[0] == nids[2] && e_nodes[1] == nids[3]) || (e_nodes[0] == nids[3] && e_nodes[1] == nids[2]);
	//CAUTION : if the edge does not belong the bloc it will return false instead of an error
}
/*----------------------------------------------------------------------------*/
bool Blocking2D::Block::isEdgeOnJ(const TCellID AID) {
	std::vector<TCellID> e_nodes = m_support->get<Edge>(AID).getIDs<Node>();
	std::vector<TCellID> nids = m_face.getIDs<Node>();
	//Same as getEdge, we use the node indices (more suitable for class changes?)
	return (e_nodes[0] == nids[1] && e_nodes[1] == nids[2]) || (e_nodes[0] == nids[2] && e_nodes[1] == nids[1])
	       || (e_nodes[0] == nids[0] && e_nodes[1] == nids[3]) || (e_nodes[0] == nids[3] && e_nodes[1] == nids[0]);
	//CAUTION : if the edge does not belong the bloc it will return false instead of an error
}
/*----------------------------------------------------------------------------*/
void Blocking2D::Block::computeT(){
	std::vector<Edge> edges = m_face.get<Edge>();
	for(int i_e = 0; i_e<4; i_e++){
		std::vector<TCellID> faces = edges[i_e].getIDs<Face>();
		if(faces.size() == 2){
			TCellID f_voisine = faces[0] == id() ? faces[1] : faces[0];

			Block voisin = m_support->block(f_voisine);

			int i_e_voisin = 0;
			for(int indice_e = 0; indice_e<4; indice_e++){
				TCellID e_voisin = voisin.m_face.getIDs<Edge>()[indice_e];
				if(e_voisin == edges[i_e].id()) i_e_voisin = indice_e;
			}

			//i_e et i_e_voisin sont les indices de l'arete interface entre les deux blocs pour chaque bloc
			//Maintenant on va vérifier si elles sont dans le même sens, en connaissant le sens et l'indice de l'arête
			//dans chaque bloc on peut retrouver la matrice de transformation

			TCellID n0 = edges[i_e].getIDs<Node>()[0];
			TCellID n1 = edges[i_e].getIDs<Node>()[1];

			int indice_0 = -1;
			int indice_1 = -1;

			for(int i = 0; i<4; i++){
				if(getNode(i).id() == n0) indice_0 = i;
				if(getNode(i).id() == n1) indice_1 = i;
			}
			//On ordonne les noeuds selon les indices
			if(indice_0 > indice_1){
				int tmp_indice = indice_0;
				indice_0 = indice_1;
				indice_1 = tmp_indice;

				TCellID tmp_id = n0;
				n0 = n1;
				n1 = tmp_id;
			}

			int indice_v0 = -1;
			int indice_v1 = -1;

			for(int i_voisin = 0; i_voisin<4; i_voisin++){
				if(voisin.getNode(i_voisin).id() == n0) indice_v0 = i_voisin;
				if(voisin.getNode(i_voisin).id() == n1) indice_v1 = i_voisin;
			}

			//Ici on a les quatre indices
			//Pour le moment on considère que la structure de blocs est conforme, les indices sont donc
			// 0,Imax,Jmax

			std::vector<int> indices = {indice_0,indice_1,indice_v0,indice_v1};

			constexpr int dim = 2; //Dans l'idée dim prend la valeur de la dimension de la structure de bloc, la on est forcément en 2D

			int coords[4][dim];

			for(int i = 0; i<4;i++){
				if(indices[i] == 0){
					coords[i][0] = 0;
					coords[i][1] = 0;
				}else if(indices[i] == 1){
					coords[i][0] = 1;
					coords[i][1] = 0;
				}else if(indices[i] == 2){
					coords[i][0] = 1;
					coords[i][1] = 1;
				}else if(indices[i] == 3){
					coords[i][0] = 0;
					coords[i][1] = 1;
				}
			}

			std::vector<int> T;
			T.resize(dim);

			for(int i = 0; i<2; i++){
				if(coords[0][i] == coords[1][i]){
					int val1 = coords[0][i];
					//On cherche la position de l'arête dans l'autre bloc
					for(int i2 = 0; i2<2; i2++){
						if(coords[2][i2] == coords[3][i2]){
							T[i] = i2 + 1;
							int val2 = coords[2][i2];
							if(val1 == val2){
								T[i] = -T[i];
							}
						}
					}
				}else{
					int val1 = coords[0][i];
					//On cherche la direction de l'arête dans l'autre bloc
					for(int i2 = 0; i2<2; i2++){
						if(coords[2][i2] != coords[3][i2]){
							T[i] = i2 + 1;
							int val2 = coords[2][i2];
							if(val1 != val2){
								T[i] = -T[i];
							}
						}
					}
				}
			}

			m_b_T.emplace(voisin.id(), T);
			std::vector<int> e_indices;
			e_indices.push_back(i_e);
			e_indices.push_back(i_e_voisin);
			m_interface_info.emplace(voisin.id(),e_indices);
		}
	}
}
/*----------------------------------------------------------------------------*/
std::map<int,std::vector<int>> Blocking2D::Block::getT(){
	return m_b_T;
}
/*----------------------------------------------------------------------------*/
std::map<int,std::vector<int>> Blocking2D::Block::getInterfaceInfo(){
	return m_interface_info;
}
