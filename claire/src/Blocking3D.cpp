//
// Created by rochec on 29/11/23.
//
/*----------------------------------------------------------------------------*/
#include <gmds/claire/Blocking3D.h>
#include <gmds/claire/Utils.h>
#include <gmds/math/TransfiniteInterpolation.h>
#include <gmds/claire/TransfiniteInterpolation_3D.h>
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
Blocking3D::Block
Blocking3D::block(const TCellID AId) {
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
Blocking3D::BlockFace
Blocking3D::blockFace(const TCellID AId) {
	return BlockFace(get<Face>(AId),this);
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

	if ( (edge_nids[0] == region_nids[0] && edge_nids[1] == region_nids[1])
	    || (edge_nids[0] == region_nids[1] && edge_nids[1] == region_nids[0])
	    || (edge_nids[0] == region_nids[3] && edge_nids[1] == region_nids[2])
	    || (edge_nids[0] == region_nids[2] && edge_nids[1] == region_nids[3])
	    || (edge_nids[0] == region_nids[4] && edge_nids[1] == region_nids[5])
	    || (edge_nids[0] == region_nids[5] && edge_nids[1] == region_nids[4])
	    || (edge_nids[0] == region_nids[7] && edge_nids[1] == region_nids[6])
	    || (edge_nids[0] == region_nids[6] && edge_nids[1] == region_nids[7]))
	{
		return m_discretization_I->value(ARegion.id());
	}
	else if ( (edge_nids[0] == region_nids[0] && edge_nids[1] == region_nids[3])
	    || (edge_nids[0] == region_nids[3] && edge_nids[1] == region_nids[0])
	    || (edge_nids[0] == region_nids[1] && edge_nids[1] == region_nids[2])
	    || (edge_nids[0] == region_nids[2] && edge_nids[1] == region_nids[1])
	    || (edge_nids[0] == region_nids[4] && edge_nids[1] == region_nids[7])
	    || (edge_nids[0] == region_nids[7] && edge_nids[1] == region_nids[4])
	    || (edge_nids[0] == region_nids[5] && edge_nids[1] == region_nids[6])
	    || (edge_nids[0] == region_nids[6] && edge_nids[1] == region_nids[5]))
	{
		return m_discretization_J->value(ARegion.id());
	}
	else if ((edge_nids[0] == region_nids[0] && edge_nids[1] == region_nids[4])
	         || (edge_nids[0] == region_nids[4] && edge_nids[1] == region_nids[0])
	         || (edge_nids[0] == region_nids[1] && edge_nids[1] == region_nids[5])
	         || (edge_nids[0] == region_nids[5] && edge_nids[1] == region_nids[1])
	         || (edge_nids[0] == region_nids[3] && edge_nids[1] == region_nids[7])
	         || (edge_nids[0] == region_nids[7] && edge_nids[1] == region_nids[3])
	         || (edge_nids[0] == region_nids[2] && edge_nids[1] == region_nids[6])
	         || (edge_nids[0] == region_nids[6] && edge_nids[1] == region_nids[2]))
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
void Blocking3D::initializeEdgesPoints()
{
	// first we check that edges are discretiezd similarly by blocks sharing
	// them
	checkDiscretizationValidity();
	// then we mesh edges
	for (auto e_id : edges()) {
			Edge ei = get<Edge>(e_id);
			int nb_subdiv = getNbDiscretization(ei.get<Region>()[0], ei);
			// We use a discrete linear interpolation along the edge
			std::vector<Node> ei_end_nodes = ei.get<Node>();
			Node n0 = ei_end_nodes[0];
			Node n1 = ei_end_nodes[1];
			math::DiscretizationScheme1DUniform d(n0.point(), n1.point(), nb_subdiv);
			auto *edge_disc = new std::vector<TCellID>();
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
void Blocking3D::initializeFacesPoints(){
	for(auto f_id:faces()){
		   Face f = get<Face>(f_id);
		   Block bi(f.get<Region>()[0], this);

		   Node n0 = f.get<Node>()[0];
		   Node n1 = f.get<Node>()[1];
		   Node n2 = f.get<Node>()[2];
		   Node n3 = f.get<Node>()[3];

		   TCellID e01_id = math::Utils::CommonEdge(this, n0.id(), n1.id());
		   TCellID e03_id = math::Utils::CommonEdge(this, n0.id(), n3.id());
		   TCellID e12_id = math::Utils::CommonEdge(this, n1.id(), n2.id());
		   TCellID e23_id = math::Utils::CommonEdge(this, n2.id(), n3.id());

		   int nb_I;
		   int nb_J;

		   if (bi.isEdgeOnI(e01_id))
		   {
			   nb_I = bi.getNbDiscretizationI();
		   }
		   else if (bi.isEdgeOnJ(e01_id))
		   {
			   nb_I = bi.getNbDiscretizationJ();
		   }
		   else
		   {
			   nb_I = bi.getNbDiscretizationK();
		   }

		   if (bi.isEdgeOnI(e03_id))
		   {
			   nb_J = bi.getNbDiscretizationI();
		   }
		   else if (bi.isEdgeOnJ(e03_id))
		   {
			   nb_J = bi.getNbDiscretizationJ();
		   }
		   else
		   {
			   nb_J = bi.getNbDiscretizationK();
		   }

		   auto* a = new Array2D<TCellID>(nb_I,nb_J);
		   Array2D<math::Point> pnts(nb_I,nb_J);

		   (*a)(0,0)           =n0.id();
		   (*a)(nb_I-1,0)      =n1.id();
		   (*a)(nb_I-1,nb_J-1) =n2.id();
		   (*a)(0,nb_J-1)      =n3.id();

		   pnts(0,0)           =n0.point();
		   pnts(nb_I-1,0)      =n1.point();
		   pnts(nb_I-1,nb_J-1) =n2.point();
		   pnts(0,nb_J-1)      =n3.point();

		   //Edge e01 = bi.getEdge(0,1);
		   Edge e01 = get<Edge>(math::Utils::CommonEdge(this, n0.id(), n1.id()));
		   std::vector<TCellID>* e01_nodes = m_edge_grids->value(e01.id());
		   if((*e01_nodes)[0]==n1.id()){
			   //need to reverse
			   std::reverse(e01_nodes->begin(),e01_nodes->end());
		   }

		   //Edge e12 = bi.getEdge(1,2);
		   Edge e12 = get<Edge>(e12_id);
		   std::vector<TCellID>* e12_nodes = m_edge_grids->value(e12.id());
		   if((*e12_nodes)[0]==n2.id()){
			   //need to reverse
			   std::reverse(e12_nodes->begin(),e12_nodes->end());
		   }

		   //Edge e32 = bi.getEdge(3,2);
		   Edge e32 = get<Edge>(e23_id);
		   std::vector<TCellID>* e32_nodes = m_edge_grids->value(e32.id());
		   if((*e32_nodes)[0]==n2.id()){
			   //need to reverse
			   std::reverse(e32_nodes->begin(),e32_nodes->end());
		   }

		   //Edge e03 = bi.getEdge(0,3);
		   Edge e03 = get<Edge>(e03_id);
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
void Blocking3D::initializeBlocksPoints()
{
	for(auto r_id:regions())
	{
		   Block bi(r_id, this);
		   auto nb_I = bi.getNbDiscretizationI();
		   auto nb_J = bi.getNbDiscretizationJ();
		   auto nb_K = bi.getNbDiscretizationK();

		   auto *a = new Array3D<TCellID>(nb_I, nb_J, nb_K);
		   Array3D<math::Point> pnts(nb_I, nb_J, nb_K);

		   Node n0 = bi.getNode(0);
		   Node n1 = bi.getNode(1);
		   Node n2 = bi.getNode(2);
		   Node n3 = bi.getNode(3);
		   Node n4 = bi.getNode(4);
		   Node n5 = bi.getNode(5);
		   Node n6 = bi.getNode(6);
		   Node n7 = bi.getNode(7);

		   (*a)(0, 0,0) = n0.id();
		   (*a)(nb_I-1, 0, 0) = n1.id();
		   (*a)(nb_I-1, nb_J-1, 0) = n2.id();
		   (*a)(0, nb_J-1, 0) = n3.id();
			(*a)(0, 0,nb_K-1) = n4.id();
		   (*a)(nb_I-1, 0, nb_K-1) = n5.id();
		   (*a)(nb_I-1, nb_J-1, nb_K-1) = n6.id();
		   (*a)(0, nb_J-1, nb_K-1) = n7.id();

		   pnts(0, 0, 0) = n0.point();
		   pnts(nb_I-1, 0, 0) = n1.point();
		   pnts(nb_I-1, nb_J-1, 0) = n2.point();
		   pnts(0, nb_J-1, 0) = n3.point();
		   pnts(0, 0, nb_K-1) = n4.point();
		   pnts(nb_I-1, 0, nb_K-1) = n5.point();
		   pnts(nb_I-1, nb_J-1, nb_K-1) = n6.point();
		   pnts(0, nb_J-1, nb_K-1) = n7.point();

		   Face f0123 = bi.getFace(0,1,2,3);
		   Array2D<TCellID> face_grid_K0 = reorientFaceGrid(f0123.id(), n0.id(),n1.id(),n2.id(),n3.id());
		   Face f4567 = bi.getFace(4,5,6,7);
		   Array2D<TCellID> face_grid_Kmax = reorientFaceGrid(f4567.id(), n4.id(),n5.id(),n6.id(),n7.id());
		   for (auto i=0;i<nb_I;i++)
		   {
			   for (auto j=0;j<nb_J;j++)
			   {
				(*a)(i,j,0) = face_grid_K0(i,j);
				pnts(i,j,0) = get<Node>(face_grid_K0(i,j)).point();
				(*a)(i,j,nb_K-1) = face_grid_Kmax(i,j);
				pnts(i,j,nb_K-1) = get<Node>(face_grid_Kmax(i,j)).point();
			   }
		   }

		   Face f0154 = bi.getFace(0,1,5,4);
		   Array2D<TCellID> face_grid_J0 = reorientFaceGrid(f0154.id(), n0.id(),n1.id(),n5.id(),n4.id());
		   Face f3267 = bi.getFace(3,2,6,7);
		   Array2D<TCellID> face_grid_Jmax = reorientFaceGrid(f3267.id(), n3.id(),n2.id(),n6.id(),n7.id());
		   for (auto i=0;i<nb_I;i++)
		   {
			   for (auto k=0;k<nb_K;k++)
			   {
				(*a)(i,0,k) = face_grid_J0(i,k);
				pnts(i,0,k) = get<Node>(face_grid_J0(i,k)).point();
				(*a)(i,nb_J-1,k) = face_grid_Jmax(i,k);
				pnts(i,nb_J-1,k) = get<Node>(face_grid_Jmax(i,k)).point();
			   }
		   }

		   Face f0374 = bi.getFace(0,3,7,4);
		   Array2D<TCellID> face_grid_I0 = reorientFaceGrid(f0374.id(), n0.id(),n3.id(),n7.id(),n4.id());
		   Face f1265 = bi.getFace(1,2,6,5);
		   Array2D<TCellID> face_grid_Imax = reorientFaceGrid(f1265.id(), n1.id(),n2.id(),n6.id(),n5.id());
		   for (auto j=0;j<nb_J;j++)
		   {
			   for (auto k=0;k<nb_K;k++)
			   {
				(*a)(0,j,k) = face_grid_I0(j,k);
				pnts(0,j,k) = get<Node>(face_grid_I0(j,k)).point();
				(*a)(nb_I-1,j,k) = face_grid_Imax(j,k);
				pnts(nb_I-1,j,k) = get<Node>(face_grid_Imax(j,k)).point();
			   }
		   }

		   TransfiniteInterpolation_3D::computeHex(pnts);

		   for(auto i=1; i<nb_I-1;i++){
			   for(auto j=1; j<nb_J-1;j++) {
					for (auto k=1; k<nb_K-1;k++) {
						Node nijk = newNode(pnts(i, j, k));
						(*a)(i, j, k) = nijk.id();
						m_embedding_dim->set(nijk.id(), 3);
						m_embedding_id->set(nijk.id(), r_id);
					}
			   }
		   }
		   m_region_grids->set(r_id,a);
	}
}
/*----------------------------------------------------------------------------*/
void Blocking3D::initializeGridPoints()
{
	initializeEdgesPoints();
	initializeFacesPoints();
	initializeBlocksPoints();
}
/*----------------------------------------------------------------------------*/
Array2D<TCellID>
Blocking3D::reorientFaceGrid(const TCellID Af_id, const TCellID An0_id, const TCellID An1_id, const TCellID An2_id, const TCellID An3_id)
{
	Array2D<TCellID>* f_grid = m_face_grids->value(Af_id);

	Face f = get<Face>(Af_id);
	std::vector<Node> f_nodes = f.get<Node>() ;

	int nb_I;
	int nb_J;

	if ( (An0_id == f_nodes[0].id() && An1_id == f_nodes[1].id())
	    || (An0_id == f_nodes[1].id() && An1_id == f_nodes[0].id())
	    || (An0_id == f_nodes[2].id() && An1_id == f_nodes[3].id())
	    || (An0_id == f_nodes[3].id() && An1_id == f_nodes[2].id()))
	{
		   nb_I = f_grid->nbLines();
		   nb_J = f_grid->nbColumns();
	}
	else
	{
		   nb_I = f_grid->nbColumns();
		   nb_J = f_grid->nbLines();
	}

	Array2D<TCellID> reOrient(nb_I, nb_J);

	if (An0_id == f_nodes[0].id())
	{
		   if (An1_id == f_nodes[1].id())
		   {
			   for (auto i=0;i<nb_I;i++)
			   {
					for (auto j=0;j<nb_J;j++)
					{
						reOrient(i,j) = (*f_grid)(i,j) ;
				   }
			   }
		   }
		   else
		   {
			   for (auto i=0;i<nb_I;i++)
			   {
					for (auto j=0;j<nb_J;j++)
					{
						reOrient(i,j) = (*f_grid)(j,i) ;
				   }
			   }
		   }
	}
	else if (An0_id == f_nodes[1].id())
	{
		   if (An1_id == f_nodes[2].id())
		   {
			   for (auto i=0;i<nb_I;i++)
			   {
					for (auto j=0;j<nb_J;j++)
					{
						reOrient(i,j) = (*f_grid)((nb_J-1)-j,i) ;
				   }
			   }
		   }
		   else
		   {
			   for (auto i=0;i<nb_I;i++)
			   {
					for (auto j=0;j<nb_J;j++)
					{
						reOrient(i,j) = (*f_grid)((nb_I-1)-i,j) ;
				   }
			   }
		   }
	}
	else if (An0_id == f_nodes[2].id())
	{
		if (An1_id == f_nodes[3].id())
		{
			for (auto i=0;i<nb_I;i++)
			{
				for (auto j=0;j<nb_J;j++)
				{
					reOrient(i,j) = (*f_grid)((nb_I-1)-i,(nb_J-1)-j) ;
				}
			}
		}
		else
		{
			for (auto i=0;i<nb_I;i++)
			{
				for (auto j=0;j<nb_J;j++)
				{
					reOrient(i,j) = (*f_grid)((nb_J-1)-j,(nb_I-1)-i) ;
				}
			}
		}

	}
	else if (An0_id == f_nodes[3].id())
	{
		if (An1_id == f_nodes[0].id())
		{
			for (auto i=0;i<nb_I;i++)
			{
				for (auto j=0;j<nb_J;j++)
				{
					reOrient(i,j) = (*f_grid)(j,(nb_I-1)-i) ;
				}
			}
		}
		else
		{
			for (auto i=0;i<nb_I;i++)
			{
				for (auto j=0;j<nb_J;j++)
				{
					reOrient(i,j) = (*f_grid)(i,(nb_J-1)-j) ;
				}
			}
		}

	}

	return reOrient;

}
/*----------------------------------------------------------------------------*/
std::vector<math::Point>
Blocking3D::getEdgeNodesPoints(const TCellID Ae_id)
{
	Edge e = get<Edge>(Ae_id);
	std::vector<TCellID>* n_ids = m_edge_grids->value(Ae_id);
	std::vector<math::Point> pts(m_edge_grids->value(Ae_id)->size());
	for (int i=0;i<m_edge_grids->value(Ae_id)->size();i++)
	{
		pts[i] = get<Node>((*n_ids)[i]).point();
	}
	return pts;
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID>
Blocking3D::computeSheet(TCellID e_id)
{
	std::vector<TCellID> sheet;
	TInt mark_isTreated = (*this).newMark<Edge>();

	// Propagation à tous les noeuds connexes qui sont sur le bord
	std::vector<TCellID> next_edges;
	next_edges.push_back(e_id);

	while (!next_edges.empty())
	{
		// On récupère un noeud de la liste de noeuds next à traiter
		TCellID current_e_id = next_edges.back();
		next_edges.pop_back();
		sheet.push_back(current_e_id);

		Edge current_e = (*this).get<Edge>(current_e_id);
		(*this).mark(current_e, mark_isTreated);	// Mark the current edge as treated

		// Get the opposite edges of the edge current_e
		std::vector<Face> current_e_faces = current_e.get<Face>();
		for (auto const &f:current_e_faces)
		{
			Edge e_opp = math::Utils::oppositeEdgeInFace(this, current_e_id, f.id());
			if (!(*this).isMarked(e_opp, mark_isTreated))
			{
				next_edges.push_back(e_opp.id());
			}
		}

	}

	(*this).unmarkAll<Edge>(mark_isTreated);
	(*this).freeMark<Edge>(mark_isTreated);

	return sheet;
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID>
Blocking3D::computeSheetBlocks(TCellID e_id)
{
	std::vector<TCellID> sheet_faces;
	std::vector<TCellID> sheet = computeSheet(e_id);
	TInt mark_isFaceOnSheet = (*this).newMark<Face>();

	for (auto e_id:sheet)
	{
		Edge e = (*this).get<Edge>(e_id);
		std::vector<Face> e_faces = e.get<Face>();
		for (auto const& f:e_faces)
		{
			(*this).mark(f, mark_isFaceOnSheet);
		}
	}

	for (auto f_id:(*this).faces())
	{
		Face b = (*this).get<Face>(f_id) ;
		if ((*this).isMarked(b, mark_isFaceOnSheet))
		{
			sheet_faces.push_back(f_id);
		}
	}

	(*this).unmarkAll<Face>(mark_isFaceOnSheet);
	(*this).freeMark<Face>(mark_isFaceOnSheet);

	return sheet;
}
/*----------------------------------------------------------------------------*/
std::vector<std::vector<TCellID>>
Blocking3D::computeAllSheets()
{
	std::vector<std::vector<TCellID>> sheets;
	TInt mark_isTreated = (*this).newMark<Edge>();

	for (auto e_id: (*this).edges())
	{
		// If an edge is marked as not treated
		if (!(*this).isMarked<Edge>(e_id, mark_isTreated))
		{
			// We compute the sheet
			std::vector<TCellID> sheet_edges_ids = computeSheet(e_id);
			// We mark all of the edges of the sheet as treated
			for (auto sheet_edge_id:sheet_edges_ids)
			{
				Edge sheet_edge = (*this).get<Edge>(sheet_edge_id);
				(*this).mark(sheet_edge, mark_isTreated);
			}
			// We add the new sheet to the vector of sheets
			sheets.push_back(sheet_edges_ids);
		}
	}

	(*this).unmarkAll<Edge>(mark_isTreated);
	(*this).freeMark<Edge>(mark_isTreated);

	return sheets;

}
/*----------------------------------------------------------------------------*/


/*============================================================================*/
//Methods for the friend class Block
/*============================================================================*/
Blocking3D::Block::Block(const Region &ARegion,  Blocking3D* ASupport)
  : m_region(ARegion), m_support(ASupport), m_grid_view(ASupport->m_region_grids->value(ARegion.id()))
{}
/*----------------------------------------------------------------------------*/
Blocking3D::Block::Block(const TCellID ARegionID, Blocking3D *ASupport)
{
	m_support= ASupport;
	m_region = m_support->get<Region>(ARegionID);
	m_grid_view=m_support->m_region_grids->value(ARegionID);
}
/*----------------------------------------------------------------------------*/
TCellID Blocking3D::Block::origin() {return m_region.getIDs<Node>()[0];}
/*----------------------------------------------------------------------------*/
Node Blocking3D::Block::getNode(const int &AIndex)
{
	if(AIndex<0 || AIndex>7){
		throw GMDSException("A 3D block has exactly 8 nodes in range [0,7]");
	}
	return m_region.get<Node>()[AIndex];
}
/*----------------------------------------------------------------------------*/
Node Blocking3D::Block::operator()(const int AI, const int AJ, const int AK)
{
	/**
	 *   J = 0    e2
	 *       ____________
	 *      |			  |
	 *  e3  |			  |  e1
	 *      |			  |
	 *      |___________|  I=0
	 *  Origin   e0
	 */
	if(AI<0 || AJ<0 || AK<0){
		throw GMDSException("Index has to be at least 0");
	}
	else if (AI>getNbDiscretizationI()-1 || AJ>getNbDiscretizationJ()-1 || AK>getNbDiscretizationK()-1)
	{
		throw GMDSException("Index out of range");
	}
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
bool Blocking3D::Block::isFaceI0(TCellID AfID)
{
	Face f = (*this).getFace(0,3,4,7);
	return (AfID == f.id());
}
/*----------------------------------------------------------------------------*/
bool Blocking3D::Block::isFaceImax(TCellID AfID)
{
	Face f = (*this).getFace(1,2,5,6);
	return (AfID == f.id());
}
/*----------------------------------------------------------------------------*/
bool Blocking3D::Block::isFaceJ0(TCellID AfID)
{
	Face f = (*this).getFace(0,1,4,5);
	return (AfID == f.id());
}
/*----------------------------------------------------------------------------*/
bool Blocking3D::Block::isFaceJmax(TCellID AfID)
{
	Face f = (*this).getFace(2,3,6,7);
	return (AfID == f.id());
}
/*----------------------------------------------------------------------------*/
bool Blocking3D::Block::isFaceK0(TCellID AfID)
{
	Face f = (*this).getFace(0,1,2,3);
	return (AfID == f.id());
}
/*----------------------------------------------------------------------------*/
bool Blocking3D::Block::isFaceKmax(TCellID AfID)
{
	Face f = (*this).getFace(4,5,6,7);
	return (AfID == f.id());
}
/*----------------------------------------------------------------------------*/
void Blocking3D::Block::computeFaceNodesPoints(const int AI, const int AJ, const int AK, const int AL)
{
	if (AI >= AJ || AJ >= AK || AK >= AL)
	{
		throw GMDSException("Error: please set the node index in an ordered way");
	}

	int nb_I = (*this).getNbDiscretizationI();
	int nb_J = (*this).getNbDiscretizationJ();
	int nb_K = (*this).getNbDiscretizationK();

	if (AI==0 && AJ == 1 && AK == 2 && AL == 3)	// Face K=0
	{
		Array2D<math::Point> pnts(nb_I, nb_J);
		for (auto i=0;i<nb_I;i++)
		{
			pnts(i,0) = (*this)(i,0,0).point();
			pnts(i,nb_J-1) = (*this)(i,nb_J-1,0).point();
		}
		for (auto j=0;j<nb_J;j++)
		{
			pnts(0,j) = (*this)(0,j,0).point();
			pnts(nb_I-1,j) = (*this)(nb_I-1,j,0).point();
		}
		math::TransfiniteInterpolation::computeQuad(pnts);
		for (auto i=1;i<nb_I-1;i++)
		{
			for (auto j=1;j<nb_J-1;j++)
			{
				(*this)(i,j,0).setPoint(pnts(i,j));
			}
		}

	}
	else if (AI==4 && AJ == 5 && AK == 6 && AL == 7)	// Face K=max
	{
		Array2D<math::Point> pnts(nb_I, nb_J);
		for (auto i=0;i<nb_I;i++)
		{
			pnts(i,0) = (*this)(i,0,nb_K-1).point();
			pnts(i,nb_J-1) = (*this)(i,nb_J-1,nb_K-1).point();
		}
		for (auto j=0;j<nb_J;j++)
		{
			pnts(0,j) = (*this)(0,j,nb_K-1).point();
			pnts(nb_I-1,j) = (*this)(nb_I-1,j,nb_K-1).point();
		}
		math::TransfiniteInterpolation::computeQuad(pnts);
		for (auto i=1;i<nb_I-1;i++)
		{
			for (auto j=1;j<nb_J-1;j++)
			{
				(*this)(i,j,nb_K-1).setPoint(pnts(i,j));
			}
		}

	}
	else if (AI==0 && AJ == 1 && AK == 4 && AL == 5)	// Face J=0
	{
		Array2D<math::Point> pnts(nb_I, nb_K);
		for (auto i=0;i<nb_I;i++)
		{
			pnts(i,0) = (*this)(i,0,0).point();
			pnts(i,nb_K-1) = (*this)(i,0,nb_K-1).point();
		}
		for (auto k=0;k<nb_K;k++)
		{
			pnts(0,k) = (*this)(0,0,k).point();
			pnts(nb_I-1,k) = (*this)(nb_I-1,0,k).point();
		}
		math::TransfiniteInterpolation::computeQuad(pnts);
		for (auto i=1;i<nb_I-1;i++)
		{
			for (auto k=1;k<nb_K-1;k++)
			{
				(*this)(i,0,k).setPoint(pnts(i,k));
			}
		}

	}
	else if (AI==2 && AJ == 3 && AK == 6 && AL == 7)	// Face J=max
	{
		Array2D<math::Point> pnts(nb_I, nb_K);
		for (auto i=0;i<nb_I;i++)
		{
			pnts(i,0) = (*this)(i,nb_J-1,0).point();
			pnts(i,nb_K-1) = (*this)(i,nb_J-1,nb_K-1).point();
		}
		for (auto k=0;k<nb_K;k++)
		{
			pnts(0,k) = (*this)(0,nb_J-1,k).point();
			pnts(nb_I-1,k) = (*this)(nb_I-1,nb_J-1,k).point();
		}
		math::TransfiniteInterpolation::computeQuad(pnts);
		for (auto i=1;i<nb_I-1;i++)
		{
			for (auto k=1;k<nb_K-1;k++)
			{
				(*this)(i,nb_J-1,k).setPoint(pnts(i,k));
			}
		}

	}
	else if (AI==0 && AJ == 3 && AK == 4 && AL == 7)	// Face I=0
	{
		Array2D<math::Point> pnts(nb_J, nb_K);
		for (auto j=0;j<nb_J;j++)
		{
			pnts(j,0) = (*this)(0,j,0).point();
			pnts(j,nb_K-1) = (*this)(0,j,nb_K-1).point();
		}
		for (auto k=0;k<nb_K;k++)
		{
			pnts(0,k) = (*this)(0,0,k).point();
			pnts(nb_J-1,k) = (*this)(0,nb_J-1,k).point();
		}
		math::TransfiniteInterpolation::computeQuad(pnts);
		for (auto j=1;j<nb_J-1;j++)
		{
			for (auto k=1;k<nb_K-1;k++)
			{
				(*this)(0,j,k).setPoint(pnts(j,k));
			}
		}

	}
	else if (AI==1 && AJ == 2 && AK == 5 && AL == 6)	// Face I=max
	{
		Array2D<math::Point> pnts(nb_J, nb_K);
		for (auto j=0;j<nb_J;j++)
		{
			pnts(j,0) = (*this)(nb_I-1,j,0).point();
			pnts(j,nb_K-1) = (*this)(nb_I-1,j,nb_K-1).point();
		}
		for (auto k=0;k<nb_K;k++)
		{
			pnts(0,k) = (*this)(nb_I-1,0,k).point();
			pnts(nb_J-1,k) = (*this)(nb_I-1,nb_J-1,k).point();
		}
		math::TransfiniteInterpolation::computeQuad(pnts);
		for (auto j=1;j<nb_J-1;j++)
		{
			for (auto k=1;k<nb_K-1;k++)
			{
				(*this)(nb_I-1,j,k).setPoint(pnts(j,k));
			}
		}

	}
}
/*----------------------------------------------------------------------------*/
void Blocking3D::Block::computeInnerBlockNodesPoints()
{
	auto nb_I = (*this).getNbDiscretizationI();
	auto nb_J = (*this).getNbDiscretizationJ();
	auto nb_K = (*this).getNbDiscretizationK();
	Array3D<math::Point> pnts(nb_I, nb_J, nb_K);
	for (auto i=0;i<nb_I;i++)
	{
		for (auto j=0;j<nb_J;j++)
		{
			pnts(i,j,0) = (*this)(i,j,0).point();
			pnts(i,j,nb_K-1) = (*this)(i,j,nb_K-1).point();
		}
	}
	for (auto i=0;i<nb_I;i++)
	{
		for (auto k=0;k<nb_K;k++)
		{
			pnts(i,0,k) = (*this)(i,0,k).point();
			pnts(i,nb_J-1,k) = (*this)(i,nb_J-1,k).point();
		}
	}
	for (auto j=0;j<nb_J;j++)
	{
		for (auto k=0;k<nb_K;k++)
		{
			pnts(0,j,k) = (*this)(0,j,k).point();
			pnts(nb_I-1,j,k) = (*this)(nb_I-1,j,k).point();
		}
	}
	TransfiniteInterpolation_3D::computeHex(pnts);
	for(auto i=1; i<nb_I-1;i++)
	{
		for(auto j=1; j<nb_J-1;j++)
		{
			for (auto k=1; k<nb_K-1;k++)
			{
				(*this)(i,j,k).setPoint(pnts(i,j,k));
			}
		}
	}
}
/*----------------------------------------------------------------------------*/



/*============================================================================*/
//Methods for the friend class BlockFace
/*============================================================================*/
Blocking3D::BlockFace::BlockFace(const Face &AFace,  Blocking3D* ASupport):
  m_face(AFace),
  m_support(ASupport),
  m_grid_view(ASupport->m_face_grids->value(AFace.id()))
{}
/*----------------------------------------------------------------------------*/
Blocking3D::BlockFace::BlockFace(const TCellID AFaceID, Blocking3D *ASupport)
{
	m_support= ASupport;
	m_face = m_support->get<Face>(AFaceID);
	m_grid_view=m_support->m_face_grids->value(AFaceID);
}
/*----------------------------------------------------------------------------*/
TCellID Blocking3D::BlockFace::origin() {return m_face.getIDs<Node>()[0];}
/*----------------------------------------------------------------------------*/
Node Blocking3D::BlockFace::getNode(const int &AIndex)
{
	if(AIndex<0 || AIndex>3){
		throw GMDSException("A 3D block face has exactly 4 nodes in range [0,3]");
	}
	return m_face.get<Node>()[AIndex];
}
/*----------------------------------------------------------------------------*/
Node Blocking3D::BlockFace::operator()(const int AI, const int AJ)
{
	/**
	 *   J = 0    e2
	 *       ____________
	 *      |			  |
	 *  e3  |			  |  e1
	 *      |			  |
	 *      |___________|  I=0
	 *  Origin   e0
	 */
	if(AI<0 || AJ<0){
		throw GMDSException("Index has to be at least 0");
	}
	else if (AI>getNbDiscretizationI()-1 || AJ>getNbDiscretizationJ()-1)
	{
		throw GMDSException("Index out of range");
	}
	return m_support->get<Node>((*m_grid_view)(AI,AJ));
}
/*----------------------------------------------------------------------------*/
Edge Blocking3D::BlockFace::getEdge(const int AI, const int AJ)  {
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
Edge Blocking3D::BlockFace::getEdgeI()
{
	return getEdge(0,1);
}
/*----------------------------------------------------------------------------*/
Edge Blocking3D::BlockFace::getEdgeJ()
{
	return getEdge(0,3);
}
/*----------------------------------------------------------------------------*/
bool Blocking3D::BlockFace::isEdgeOnI(const TCellID AID)
{
	std::vector<TCellID> e_nodes = m_support->get<Edge>(AID).getIDs<Node>();
	std::vector<TCellID> nids = m_face.getIDs<Node>();
	//Same as getEdge, we use the node indices (more suitable for class changes?)
	return (e_nodes[0] == nids[0] && e_nodes[1] == nids[1]) || (e_nodes[0] == nids[1] && e_nodes[1] == nids[0])
	       || (e_nodes[0] == nids[2] && e_nodes[1] == nids[3]) || (e_nodes[0] == nids[3] && e_nodes[1] == nids[2]);
	//CAUTION : if the edge does not belong the bloc it will return false instead of an error
}
/*----------------------------------------------------------------------------*/
bool Blocking3D::BlockFace::isEdgeOnJ(const TCellID AID)
{
	std::vector<TCellID> e_nodes = m_support->get<Edge>(AID).getIDs<Node>();
	std::vector<TCellID> nids = m_face.getIDs<Node>();
	//Same as getEdge, we use the node indices (more suitable for class changes?)
	return (e_nodes[0] == nids[1] && e_nodes[1] == nids[2]) || (e_nodes[0] == nids[2] && e_nodes[1] == nids[1])
	       || (e_nodes[0] == nids[0] && e_nodes[1] == nids[3]) || (e_nodes[0] == nids[3] && e_nodes[1] == nids[0]);
	//CAUTION : if the edge does not belong the bloc it will return false instead of an error
}
/*----------------------------------------------------------------------------*/
int Blocking3D::BlockFace::getNbDiscretizationI() const{
	return m_grid_view->nbLines();
}
/*----------------------------------------------------------------------------*/
int Blocking3D::BlockFace::getNbDiscretizationJ() const{
	return m_grid_view->nbColumns();
}
/*----------------------------------------------------------------------------*/



/*============================================================================*/
//Methods for the friend class BlockEdge
/*============================================================================*/
Blocking3D::BlockEdge::BlockEdge(const Edge &AEdge,  Blocking3D* ASupport):
  m_edge(AEdge),
  m_support(ASupport),
  m_grid_view(ASupport->m_edge_grids->value(AEdge.id()))
{}
/*----------------------------------------------------------------------------*/
Blocking3D::BlockEdge::BlockEdge(const TCellID AEdgeID, Blocking3D *ASupport)
{
	m_support= ASupport;
	m_edge = m_support->get<Edge>(AEdgeID);
	m_grid_view=m_support->m_edge_grids->value(AEdgeID);
}
/*----------------------------------------------------------------------------*/
TCellID Blocking3D::BlockEdge::origin() {return m_edge.getIDs<Node>()[0];}
/*----------------------------------------------------------------------------*/
Node Blocking3D::BlockEdge::getNode(const int &AIndex)
{
	if(AIndex<0 || AIndex>1){
		throw GMDSException("A 3D block edge has exactly 2 nodes in range [0,1]");
	}
	return m_edge.get<Node>()[AIndex];
}
/*----------------------------------------------------------------------------*/
Node Blocking3D::BlockEdge::operator()(const int AI)
{
	/**
	 *   J = 0    e2
	 *       ____________
	 *      |			  |
	 *  e3  |			  |  e1
	 *      |			  |
	 *      |___________|  I=0
	 *  Origin   e0
	 */
	if(AI<0){
		throw GMDSException("Index has to be at least 0");
	}
	else if (AI>getNbDiscretizationI()-1)
	{
		throw GMDSException("Index out of range");
	}
	return m_support->get<Node>((*m_grid_view)[AI]);
}
/*----------------------------------------------------------------------------*/
int Blocking3D::BlockEdge::getNbDiscretizationI() const{
	return m_grid_view->size();
}
/*----------------------------------------------------------------------------*/