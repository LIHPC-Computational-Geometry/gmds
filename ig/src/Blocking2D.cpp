/*----------------------------------------------------------------------------*/
#include <gmds/ig/Blocking2D.h>
#include <gmds/math/TransfiniteInterpolation.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
Blocking2D::Blocking2D() 
: Mesh(MeshModel(DIM3|F|E|N|F2N|E2N|F2E|E2F|N2E))
{
    m_embedding_dim = newVariable<int, GMDS_NODE>("embedding_dim");
    m_embedding_id = newVariable<TCellID , GMDS_NODE>("embedding_id");
    m_discretization_I= newVariable<int, GMDS_FACE>("discretizationI");
    m_discretization_J= newVariable<int, GMDS_FACE>("discretizationJ");
    m_face_grids = newVariable<Array2D<TCellID>*, GMDS_FACE >("grid_ref");
    m_edge_grids = newVariable<std::vector<TCellID>*, GMDS_EDGE >("grid_ref");

}
/*----------------------------------------------------------------------------*/
Blocking2D::Blocking2D(const Mesh& AMesh): Mesh(MeshModel(AMesh.getModel())){

	// for now we assume that the mesh hasnt already been used for a blocking

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
		b.setNbDiscretizationI(11); //Setting default discretization to 11 for all edges
		b.setNbDiscretizationJ(11);
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
    return Block(b, this);
}
/*----------------------------------------------------------------------------*/
Blocking2D::Block Blocking2D::newBlock(const Node &AN1, const Node &AN2,
                                       const Node &AN3, const Node &AN4){
    return newBlock(AN1.id(),AN2.id(),AN3.id(),AN4.id());
}
/*----------------------------------------------------------------------------*/
Blocking2D::Block  Blocking2D::block(const TCellID AId) {
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
    bool found_common=false;
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
    return m_support->get<Node>((*m_grid_view)(AI,AJ));
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
    bool found_common=false;
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
    math::Vector3d v(nids[0].point(),nids[1].point());
    return v.getNormalize();
}
/*----------------------------------------------------------------------------*/
math::Vector3d Blocking2D::Block::getUnitVectorJ() {
    std::vector<Node> nids = m_face.get<Node>();
    math::Vector3d v(nids[0].point(),nids[3].point());
    return v.getNormalize();
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

	return std::vector<TCellID>();
}
/*----------------------------------------------------------------------------*/
