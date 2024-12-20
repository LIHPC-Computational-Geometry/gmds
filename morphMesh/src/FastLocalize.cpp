/*------------------------------------------------------------------------*/
#include "gmds/morphMesh/FastLocalize.h"
//#include <gmds/claire/Utils.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/
FastLocalize::FastLocalize(const std::vector<Node> &ANodes): m_nodes(ANodes){
	int nb_pnts 	= m_nodes.size();
	int	k			= 20;      // max number of nearest neighbors
	int	dim		= 3;       // dimension
	int	maxPts	= nb_pnts; // maximum number of data points

	int					nPts;						// actual number of data points

	m_queryPt = annAllocPt(3);
	m_dataPts = annAllocPts(maxPts, dim);			// allocate data points
	m_nnIdx = new ANNidx[k];						// allocate near neigh indices
	m_dists = new ANNdist[k];						// allocate near neighbor dists

	//========================================================
	// (1) Fill in the  ANN structure for storing points
	//
	// Important: Points in APnts and dataPnts are stored with
	// same index.
	//========================================================
	nPts = 0;
	for(auto const &n: m_nodes){
		math::Point p = n.point();
		m_dataPts[nPts][0] = p.X();
		m_dataPts[nPts][1] = p.Y();
		m_dataPts[nPts][2] = p.Z();
		m_ann2gmds_id[nPts]=n.id();
		nPts++;
	}

	//========================================================
	// (2) Build the search structure
	//========================================================
	m_kdTree = new ANNkd_tree(m_dataPts,	// the data points
	                          nPts,		// number of points
	                          dim);		// dimension of space
}
/*------------------------------------------------------------------------*/
FastLocalize::~FastLocalize(){
	delete [] m_nnIdx;							// clean things up
	delete [] m_dists;
	delete m_kdTree;
	annDeallocPt(m_queryPt);
	annDeallocPts(m_dataPts);
	annClose();									// done with ANN

}
/*------------------------------------------------------------------------*/
TCellID
FastLocalize::find(const math::Point &APoint)
{
	int	k			= 5;      // max number of nearest neighbors
	m_queryPt[0]=APoint.X();
	m_queryPt[1]=APoint.Y();
	m_queryPt[2]=APoint.Z();

	m_kdTree->annkSearch(		// search
	   m_queryPt,// query point
	   k,
	   m_nnIdx,
	   m_dists,
	   0.01);
	return m_ann2gmds_id[m_nnIdx[0]];
}
/*------------------------------------------------------------------------*/