/*------------------------------------------------------------------------*/
#include <gmds/claire/FastLocalize.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/
FastLocalize::FastLocalize(Mesh *AMesh): m_mesh(AMesh) {
	int nb_pnts 	= m_mesh->getNbNodes();
	int	k			= 20;      // max number of nearest neighbors
	int	dim		= 2;       // dimension
	int	maxPts	= nb_pnts; // maximum number of data points

	int					nPts;						// actual number of data points
	ANNpointArray		dataPts;					// data points
	ANNpoint				queryPt;					// query point

	m_queryPt = annAllocPt(2);


	dataPts = annAllocPts(maxPts, dim);			// allocate data points
	m_nnIdx = new ANNidx[k];						// allocate near neigh indices
	m_dists = new ANNdist[k];						// allocate near neighbor dists

	//========================================================
	// (1) Fill in the  ANN structure for storing points
	//
	// Important: Points in APnts and dataPnts are stored with
	// same index.
	//========================================================
	nPts = 0;
	for(auto n_id: m_mesh->nodes()){
		math::Point p = m_mesh->get<Node>(n_id).point();
		dataPts[nPts][0] = p.X();
		dataPts[nPts][1] = p.Y();
		m_ann2gmds_id[nPts]=n_id;
		nPts++;
	};

	//========================================================
	// (2) Build the search structure
	//========================================================
	m_kdTree = new ANNkd_tree(dataPts,	// the data points
	                          nPts,		// number of points
	                          dim);		// dimension of space
}
/*------------------------------------------------------------------------*/
FastLocalize::~FastLocalize(){
	delete [] m_nnIdx;							// clean things up
	delete [] m_dists;
	delete m_kdTree;
	annDeallocPt(m_queryPt);
	annClose();									// done with ANN

}
/*------------------------------------------------------------------------*/
Cell::Data
FastLocalize::find(const math::Point &APoint)
{
	int	k			= 5;      // max number of nearest neighbors
	m_queryPt[0]=APoint.X();
	m_queryPt[1]=APoint.Y();

	m_kdTree->annkSearch(		// search
	   m_queryPt,// query point
	   k,
	   m_nnIdx,
	   m_dists,
	   0.01);
	return Cell::Data(0,m_ann2gmds_id[m_nnIdx[0]]);
}
