//
// Created by rochec on 29/11/23.
//
/*----------------------------------------------------------------------------*/
#include <gmds/claire/Blocking3D.h>
//#include <gmds/math/TransfiniteInterpolation.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
Blocking3D::Blocking3D()
  : Mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F | F2R |
                 F2E | E2F | R2E | E2R | N2R | N2F | N2E))
{
	m_embedding_dim = newVariable<int, GMDS_NODE>("embedding_dim");
	m_embedding_id = newVariable<TCellID , GMDS_NODE>("embedding_id");
	m_discretization_I= newVariable<int, GMDS_REGION>("discretizationI");
	m_discretization_J= newVariable<int, GMDS_REGION>("discretizationJ");
	m_discretization_K= newVariable<int, GMDS_REGION>("discretizationK");
	m_region_grids = newVariable<Array3D<TCellID>*, GMDS_REGION >("grid_ref");
	m_face_grids = newVariable<Array2D<TCellID>*, GMDS_FACE >("grid_ref");
	m_edge_grids = newVariable<std::vector<TCellID>*, GMDS_EDGE >("grid_ref");

}
/*----------------------------------------------------------------------------*/
Blocking3D::Blocking3D(const Mesh& AMesh): Mesh(AMesh.getModel()){

	// for now, we assume that the mesh hasn't already been used for a blocking

	m_embedding_dim = newVariable<int, GMDS_NODE>("embedding_dim");
	m_embedding_id = newVariable<TCellID , GMDS_NODE>("embedding_id");
	m_discretization_I= newVariable<int, GMDS_REGION>("discretizationI");
	m_discretization_J= newVariable<int, GMDS_REGION>("discretizationJ");
	m_discretization_K= newVariable<int, GMDS_REGION>("discretizationK");
	m_region_grids = newVariable<Array3D<TCellID>*, GMDS_REGION >("grid_ref");
	m_face_grids = newVariable<Array2D<TCellID>*, GMDS_FACE >("grid_ref");
	m_edge_grids = newVariable<std::vector<TCellID>*, GMDS_EDGE >("grid_ref");

	// Building the blocking from the mesh
	// Maybe make an independent function for the building, with different types of discretization method etc..

	//map the relation from mesh vertices to blocking vertices
	std::map<TCellID,TCellID> n2n;
	for (auto n : AMesh.nodes()) {
		Node n_b = newBlockCorner(AMesh.get<Node>(n).point());
		n2n[n] = n_b.id();
	}
	for (auto r : AMesh.regions()) {
		Region region = AMesh.get<Region>(r);
		std::vector<TCellID> r_nodes;
		region.getIDs<Node>(r_nodes);
		Block b = newBlock(n2n[r_nodes[0]],n2n[r_nodes[1]],n2n[r_nodes[2]],n2n[r_nodes[3]],
		                   n2n[r_nodes[4]],n2n[r_nodes[5]],n2n[r_nodes[6]],n2n[r_nodes[7]]);
	}


}
/*----------------------------------------------------------------------------*/
Blocking3D::~Blocking3D()
{
	//by deleting the mesh, we also delete the memory space pointed by m_embedding
	// but we need to delete explicitly allocated data
	for(auto r_id:regions()){
		delete m_region_grids->value(r_id);
		m_region_grids->set(r_id,NULL);
	}
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
Node Blocking3D::newBlockCorner(const math::Point &APnt) {
	Node n = newNode(APnt);
	m_embedding_dim->set(n.id(),0);
	m_embedding_id->set(n.id(),n.id());
	return n;
}
/*----------------------------------------------------------------------------*/
Node Blocking3D::newBlockCorner(const TCoord AX,const TCoord AY, const TCoord AZ) {
	Node n = newNode(math::Point(AX,AY,AZ));
	m_embedding_dim->set(n.id(),0);
	m_embedding_id->set(n.id(),n.id());
	return n;
}
/*----------------------------------------------------------------------------*/
Blocking3D::Block Blocking3D::newBlock(const TCellID AN1, const TCellID AN2,
                     const TCellID AN3, const TCellID AN4,
                     const TCellID AN5, const TCellID AN6,
                     const TCellID AN7, const TCellID AN8)
{
	Region b = newHex(AN1, AN2, AN3, AN4, AN5, AN6, AN7, AN8);

	// F2R and R2F
	Face f1 = get<Face>(getFace(AN1, AN2, AN3, AN4)) ;
	b.add(f1); f1.add(b);
	Face f2 = get<Face>(getFace(AN5, AN6, AN7, AN8)) ;
	b.add(f2); f2.add(b);
	Face f3 = get<Face>(getFace(AN1, AN2, AN6, AN5)) ;
	b.add(f3); f3.add(b);
	Face f4 = get<Face>(getFace(AN4, AN3, AN7, AN8)) ;
	b.add(f4); f4.add(b);
	Face f5 = get<Face>(getFace(AN1, AN5, AN8, AN4)) ;
	b.add(f5); f5.add(b);
	Face f6 = get<Face>(getFace(AN2, AN6, AN7, AN3)) ;
	b.add(f6); f6.add(b);

	// E2R and R2E
	Edge e1 = get<Edge>(getEdge(AN1,AN2));
	b.add(e1); e1.add(b);
	Edge e2 = get<Edge>(getEdge(AN2,AN3));
	b.add(e2); e2.add(b);
	Edge e3 = get<Edge>(getEdge(AN3,AN4));
	b.add(e3); e3.add(b);
	Edge e4 = get<Edge>(getEdge(AN4,AN1));
	b.add(e4); e4.add(b);
	Edge e5 = get<Edge>(getEdge(AN5,AN6));
	b.add(e5); e5.add(b);
	Edge e6 = get<Edge>(getEdge(AN6,AN7));
	b.add(e6); e6.add(b);
	Edge e7 = get<Edge>(getEdge(AN7,AN8));
	b.add(e7); e7.add(b);
	Edge e8 = get<Edge>(getEdge(AN8,AN5));
	b.add(e8); e8.add(b);
	Edge e9 = get<Edge>(getEdge(AN1,AN5));
	b.add(e9); e9.add(b);
	Edge e10 = get<Edge>(getEdge(AN2,AN6));
	b.add(e10); e10.add(b);
	Edge e11 = get<Edge>(getEdge(AN3,AN7));
	b.add(e11); e11.add(b);
	Edge e12 = get<Edge>(getEdge(AN4,AN8));
	b.add(e12); e12.add(b);

	// Update N2R
	get<Node>(AN1).add(b);
	get<Node>(AN2).add(b);
	get<Node>(AN3).add(b);
	get<Node>(AN4).add(b);
	get<Node>(AN5).add(b);
	get<Node>(AN6).add(b);
	get<Node>(AN7).add(b);
	get<Node>(AN8).add(b);

	return Block(b, this);
}
/*----------------------------------------------------------------------------*/
Blocking3D::Block Blocking3D::newBlock(const Node &AN1, const Node &AN2,
                     const Node &AN3, const Node &AN4,
                     const Node &AN5, const Node &AN6,
                     const Node &AN7, const Node &AN8){
	return newBlock(AN1.id(),AN2.id(),AN3.id(),AN4.id(),
	                AN5.id(),AN6.id(),AN7.id(),AN8.id());
}
/*----------------------------------------------------------------------------*/
Blocking3D::Block  Blocking3D::
   block(const TCellID AId) {
	return Block(get<Region>(AId),this);
}
/*----------------------------------------------------------------------------*/
std::vector<Blocking3D::Block> Blocking3D::allBlocks() {
	std::vector<Blocking3D::Block> blocks;
	for(auto r_id:regions()){
		blocks.push_back(block(r_id));
	}
	return blocks;
}
/*----------------------------------------------------------------------------*/
TCellID Blocking3D::getEdge(const TCellID AN1, const TCellID AN2) {
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
TCellID Blocking3D::getFace(const TCellID AN1, const TCellID AN2, const TCellID AN3, const TCellID AN4) {
	std::vector<TCellID> faces_1 = get<Node>(AN1).getIDs<Face>();
	std::vector<TCellID> faces_2 = get<Node>(AN2).getIDs<Face>();
	std::vector<TCellID> faces_3 = get<Node>(AN3).getIDs<Face>();

	for(auto f1:faces_1){
		for(auto f2:faces_2){
			for(auto f3:faces_3) {
				if (f1 == f2 && f2 == f3) {
					// found common face
					return f1;
				}
			}
		}
	}
	//there is not a common face, so we add it

	Face f = newQuad(AN1,AN2,AN3,AN4);

	// F2E and E2F
	Edge e1 = get<Edge>(getEdge(AN1,AN2));
	f.add(e1); e1.add(f);
	Edge e2 = get<Edge>(getEdge(AN2,AN3));
	f.add(e2); e2.add(f);
	Edge e3 = get<Edge>(getEdge(AN3,AN4));
	f.add(e3); e3.add(f);
	Edge e4 = get<Edge>(getEdge(AN4,AN1));
	f.add(e4); e4.add(f);

	// F2N updates
	get<Node>(AN1).add(f);
	get<Node>(AN2).add(f);
	get<Node>(AN3).add(f);
	get<Node>(AN4).add(f);

	return f.id();
}
/*----------------------------------------------------------------------------*/
int Blocking3D::getNbDiscretization(const Region& ARegion, const Edge &AEdge) const {
	//find the edge position in the region
	std::vector<TCellID> region_nids = ARegion.getIDs<Node>();
	std::vector<TCellID> edge_nids = AEdge.getIDs<Node>();

	if ( (edge_nids[0] == region_nids[0] && edge_nids[0] == region_nids[1])
	    || (edge_nids[0] == region_nids[1] && edge_nids[0] == region_nids[0])
	    || (edge_nids[0] == region_nids[3] && edge_nids[0] == region_nids[2])
	    || (edge_nids[0] == region_nids[2] && edge_nids[0] == region_nids[3])
	    || (edge_nids[0] == region_nids[4] && edge_nids[0] == region_nids[5])
	    || (edge_nids[0] == region_nids[5] && edge_nids[0] == region_nids[4])
	    || (edge_nids[0] == region_nids[7] && edge_nids[0] == region_nids[6])
	    || (edge_nids[0] == region_nids[6] && edge_nids[0] == region_nids[7]))
	{
		return m_discretization_I->value(ARegion.id());
	}
	else if ( (edge_nids[0] == region_nids[0] && edge_nids[0] == region_nids[3])
	    || (edge_nids[0] == region_nids[3] && edge_nids[0] == region_nids[0])
	    || (edge_nids[0] == region_nids[1] && edge_nids[0] == region_nids[2])
	    || (edge_nids[0] == region_nids[2] && edge_nids[0] == region_nids[1])
	    || (edge_nids[0] == region_nids[4] && edge_nids[0] == region_nids[7])
	    || (edge_nids[0] == region_nids[7] && edge_nids[0] == region_nids[4])
	    || (edge_nids[0] == region_nids[5] && edge_nids[0] == region_nids[6])
	    || (edge_nids[0] == region_nids[6] && edge_nids[0] == region_nids[5]))
	{
		return m_discretization_J->value(ARegion.id());
	}
	else if ((edge_nids[0] == region_nids[0] && edge_nids[0] == region_nids[4])
	         || (edge_nids[0] == region_nids[4] && edge_nids[0] == region_nids[0])
	         || (edge_nids[0] == region_nids[1] && edge_nids[0] == region_nids[5])
	         || (edge_nids[0] == region_nids[5] && edge_nids[0] == region_nids[1])
	         || (edge_nids[0] == region_nids[3] && edge_nids[0] == region_nids[7])
	         || (edge_nids[0] == region_nids[7] && edge_nids[0] == region_nids[3])
	         || (edge_nids[0] == region_nids[2] && edge_nids[0] == region_nids[6])
	         || (edge_nids[0] == region_nids[6] && edge_nids[0] == region_nids[2]))
	{
		return m_discretization_K->value(ARegion.id());
	}
	else
	{
		throw GMDSException("Error: edge not adjacent to the region");
	}
}
/*----------------------------------------------------------------------------*/
bool Blocking3D::checkDiscretizationValidity() const {
   for(auto e_id:edges()){
      Edge ei = get<Edge>(e_id);
      std::vector<Region> adj_regions = ei.get<Region>();
      if(adj_regions.size()>1){
         //we verify that we have the same discretization scheme on each regions
         int nb_subdiv0 = getNbDiscretization(adj_regions[0],ei);
			for (auto i=1;i<adj_regions.size();i++)
			{
				int nb_subdiv1 = getNbDiscretization(adj_regions[i],ei);
				if(nb_subdiv0!=nb_subdiv1){
					return false;
				}
			}
      }
   }

return true;
}
/*----------------------------------------------------------------------------*/




/*============================================================================*/
//Methods for the friend class Block
/*============================================================================*/
Blocking3D::Block::Block(const Region &ARegion,  Blocking3D* ASupport)
  : m_region(ARegion), m_support(ASupport), m_grid_view(ASupport->m_region_grids->value(ARegion.id()))
{}
/*----------------------------------------------------------------------------*/
Blocking3D::Block::Block(const TCellID ARegionID, Blocking3D *ASupport) {
	m_support= ASupport;
	m_region = m_support->get<Region>(ARegionID);
	m_grid_view=m_support->m_region_grids->value(ARegionID);
}
/*----------------------------------------------------------------------------*/
TCellID Blocking3D::Block::origin() {return m_region.getIDs<Node>()[0];}
/*----------------------------------------------------------------------------*/
Node Blocking3D::Block::getNode(const int &AIndex) {
	if(AIndex<0 || AIndex>7){
		throw GMDSException("A 3D block has exactly 8 nodes in range [0,7]");
	}
	return m_region.get<Node>()[AIndex];
}
/*----------------------------------------------------------------------------*/
Node Blocking3D::Block::operator()(const int AI, const int AJ, const int AK) {
	/**
	 *   J = 0    e2
	 *       ____________
	 *      |			  |
	 *  e3  |			  |  e1
	 *      |			  |
	 *      |___________|  I=0
	 *  Origin   e0
	 */
	return m_support->get<Node>((*m_grid_view)(AI,AJ,AK));
}
/*----------------------------------------------------------------------------*/
Edge Blocking3D::Block::getEdgeI() {
	return getEdge(0,1);
}
/*----------------------------------------------------------------------------*/
Edge Blocking3D::Block::getEdgeJ() {
	return getEdge(0,3);
}
/*----------------------------------------------------------------------------*/
Edge Blocking3D::Block::getEdgeK() {
	return getEdge(0,4);
}
/*----------------------------------------------------------------------------*/
Edge Blocking3D::Block::getEdge(const int AI, const int AJ)  {
	if(AI==AJ)
		throw GMDSException("Error in local numbering for getting an edge: same node value");
	if(AI<0 || AJ<0)
		throw GMDSException("Error in local numbering for getting an edge: negative local numbering");
	if(AI>7 || AJ>7)
		throw GMDSException("Error in local numbering for getting an edge: a local number can not exceed 7");
	std::vector<TCellID> nids = m_region.getIDs<Node>();
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
Face Blocking3D::Block::getFace(const int AI, const int AJ, const int AK, const int AL)
{
	if(AI==AJ || AI==AK || AJ==AK)
		throw GMDSException("Error in local numbering for getting an edge: same node value");
	if(AI<0 || AJ<0 || AK<0 || AL<0)
		throw GMDSException("Error in local numbering for getting an edge: negative local numbering");
	if(AI>7 || AJ>7 || AK>7 || AL>7)
		throw GMDSException("Error in local numbering for getting an edge: a local number can not exceed 7");
	std::vector<TCellID> nids = m_region.getIDs<Node>();
	std::vector<TCellID> faces_1 = m_support->get<Node>(nids[AI]).getIDs<Face>();
	std::vector<TCellID> faces_2 = m_support->get<Node>(nids[AJ]).getIDs<Face>();
	std::vector<TCellID> faces_3 = m_support->get<Node>(nids[AK]).getIDs<Face>();

	for(auto f1:faces_1){
		for(auto f2:faces_2){
			for(auto f3:faces_3) {
				if (f1 == f2 && f2 == f3) {
					// found common face
					return m_support->get<Face>(f1);
				}
			}
		}
	}
	throw GMDSException("error no edge found!");
}
/*----------------------------------------------------------------------------*/
void Blocking3D::Block::setNbDiscretizationI(const int AN) {
	if(AN<=0)
		throw GMDSException("discretization value must be strictly positive");

	m_support->m_discretization_I->set(m_region.id(),AN);
}
/*----------------------------------------------------------------------------*/
void Blocking3D::Block::setNbDiscretizationJ(const int AN) {
	if(AN<=0)
		throw GMDSException("discretization value must be strictly positive");

	m_support->m_discretization_J->set(m_region.id(),AN);
}
/*----------------------------------------------------------------------------*/
void Blocking3D::Block::setNbDiscretizationK(const int AN) {
	if(AN<=0)
		throw GMDSException("discretization value must be strictly positive");

	m_support->m_discretization_K->set(m_region.id(),AN);
}
/*----------------------------------------------------------------------------*/
int Blocking3D::Block::getNbDiscretizationI() const{
	return m_support->m_discretization_I->value(m_region.id());
}
/*----------------------------------------------------------------------------*/
int Blocking3D::Block::getNbDiscretizationJ() const{
	return m_support->m_discretization_J->value(m_region.id());
}
/*----------------------------------------------------------------------------*/
int Blocking3D::Block::getNbDiscretizationK() const{
	return m_support->m_discretization_K->value(m_region.id());
}
/*----------------------------------------------------------------------------*/
math::Vector3d Blocking3D::Block::getUnitVectorI() {
	std::vector<Node> nids = m_region.get<Node>();
	math::Vector3d v=nids[1].point()-nids[0].point();
	return v.getNormalize();
}
/*----------------------------------------------------------------------------*/
math::Vector3d Blocking3D::Block::getUnitVectorJ() {
	std::vector<Node> nids = m_region.get<Node>();
	math::Vector3d v=nids[3].point()-nids[0].point();
	return v.getNormalize();
}
/*----------------------------------------------------------------------------*/
math::Vector3d Blocking3D::Block::getUnitVectorK() {
	std::vector<Node> nids = m_region.get<Node>();
	math::Vector3d v=nids[4].point()-nids[0].point();
	return v.getNormalize();
}
/*----------------------------------------------------------------------------*/
std::tuple<int,int,int> Blocking3D::Block::getIndices(const TCellID AID){
	//Maybe kind of heavy method
	//Should find a smarter way to find the point
	for(int i = 0; i<getNbDiscretizationI(); i++){
		for(int j = 0; j<getNbDiscretizationJ(); j++){
			for(int k = 0; k<getNbDiscretizationK(); j++) {
				if ((*m_grid_view)(i, j, k) == AID) {
					return {i, j, k};
				}
			}
		}
	}
	throw GMDSException("ERROR: Node not in the block grid.");
}
/*----------------------------------------------------------------------------*/
bool Blocking3D::Block::isEdgeOnI(const TCellID AID) {
	std::vector<TCellID> e_nodes = m_support->get<Edge>(AID).getIDs<Node>();
	std::vector<TCellID> nids = m_region.getIDs<Node>();
	//Same as getEdge, we use the node indices (more suitable for class changes?)
	return (e_nodes[0] == nids[0] && e_nodes[1] == nids[1]) || (e_nodes[0] == nids[1] && e_nodes[1] == nids[0])
	       || (e_nodes[0] == nids[2] && e_nodes[1] == nids[3]) || (e_nodes[0] == nids[3] && e_nodes[1] == nids[2])
	       || (e_nodes[0] == nids[4] && e_nodes[1] == nids[5]) || (e_nodes[0] == nids[5] && e_nodes[1] == nids[4])
	       || (e_nodes[0] == nids[6] && e_nodes[1] == nids[7]) || (e_nodes[0] == nids[7] && e_nodes[1] == nids[6]);
	//CAUTION : if the edge does not belong the bloc it will return false instead of an error
}
/*----------------------------------------------------------------------------*/
bool Blocking3D::Block::isEdgeOnJ(const TCellID AID) {
	std::vector<TCellID> e_nodes = m_support->get<Edge>(AID).getIDs<Node>();
	std::vector<TCellID> nids = m_region.getIDs<Node>();
	//Same as getEdge, we use the node indices (more suitable for class changes?)
	return (e_nodes[0] == nids[1] && e_nodes[1] == nids[2]) || (e_nodes[0] == nids[2] && e_nodes[1] == nids[1])
	       || (e_nodes[0] == nids[0] && e_nodes[1] == nids[3]) || (e_nodes[0] == nids[3] && e_nodes[1] == nids[0])
	       || (e_nodes[0] == nids[5] && e_nodes[1] == nids[6]) || (e_nodes[0] == nids[6] && e_nodes[1] == nids[5])
	       || (e_nodes[0] == nids[4] && e_nodes[1] == nids[7]) || (e_nodes[0] == nids[7] && e_nodes[1] == nids[4]);
	//CAUTION : if the edge does not belong the bloc it will return false instead of an error
}
/*----------------------------------------------------------------------------*/
bool Blocking3D::Block::isEdgeOnK(const TCellID AID) {
	std::vector<TCellID> e_nodes = m_support->get<Edge>(AID).getIDs<Node>();
	std::vector<TCellID> nids = m_region.getIDs<Node>();
	//Same as getEdge, we use the node indices (more suitable for class changes?)
	return (e_nodes[0] == nids[0] && e_nodes[1] == nids[4]) || (e_nodes[0] == nids[4] && e_nodes[1] == nids[0])
	       || (e_nodes[0] == nids[1] && e_nodes[1] == nids[5]) || (e_nodes[0] == nids[5] && e_nodes[1] == nids[1])
	       || (e_nodes[0] == nids[3] && e_nodes[1] == nids[7]) || (e_nodes[0] == nids[7] && e_nodes[1] == nids[3])
	       || (e_nodes[0] == nids[2] && e_nodes[1] == nids[6]) || (e_nodes[0] == nids[6] && e_nodes[1] == nids[2]);
	//CAUTION : if the edge does not belong the bloc it will return false instead of an error
}
/*----------------------------------------------------------------------------*/