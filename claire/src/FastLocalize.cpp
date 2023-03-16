/*------------------------------------------------------------------------*/
#include <gmds/claire/FastLocalize.h>
#include <gmds/claire/Utils.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/
FastLocalize::FastLocalize(Mesh *AMesh): m_mesh(AMesh) {
	int nb_pnts 	= m_mesh->getNbNodes();
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
	for(auto n_id: m_mesh->nodes()){
		math::Point p = m_mesh->get<Node>(n_id).point();
		m_dataPts[nPts][0] = p.X();
		m_dataPts[nPts][1] = p.Y();
		m_dataPts[nPts][2] = p.Z();
		m_ann2gmds_id[nPts]=n_id;
		nPts++;
	};

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
Cell::Data
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
	return Cell::Data(0,m_ann2gmds_id[m_nnIdx[0]]);
}
/*------------------------------------------------------------------------*/
TCellID
FastLocalize::findTetra(const math::Point &APoint)
{
	TCellID tetra_id(NullID);

	Cell::Data data = (*this).find(APoint) ;
	if (data.dim==0)
	{
		Node closest_node = m_mesh->get<Node>(data.id);
		std::vector<Region> closest_regions = closest_node.get<Region>();
		if (closest_regions.empty())
		{
			std::cout << "ATTENTION FastLocalize: in findTetra, the closest node to the point is connected to 0 regions." << std::endl;
		}
		bool tetraFound(false);

		// Test on regions of the first ring
		for (auto r:closest_node.get<Region>())
		{
			if (!tetraFound)
			{
				tetraFound = math::Utils::isInTetra(r.get<Node>()[0].point(), r.get<Node>()[1].point(), r.get<Node>()[2].point(),
				                                    r.get<Node>()[3].point(), APoint);
				if (tetraFound)
				{
					tetra_id = r.id();
				}
			}
		}

		// Test on regions of the second ring, if the tetra is not found on the first ring
		if (!tetraFound) {
			std::vector<TCellID> regions_to_test;
			int mark_RegionsTreated = m_mesh->newMark<Region>();

			// Compute the regions to test. Here, we choose the regions on the 2 ring of the node.
			for (auto r : closest_node.get<Region>())
			{
				regions_to_test.push_back(r.id());
				m_mesh->mark(r, mark_RegionsTreated);
				for (auto r_node : r.get<Node>())
				{
					for (auto r_2 : r_node.get<Region>())
					{
						if (!m_mesh->isMarked(r_2, mark_RegionsTreated))
						{
							regions_to_test.push_back(r_2.id());
							m_mesh->mark(r_2, mark_RegionsTreated);
						}
					}
				}
			}

			m_mesh->unmarkAll<Region>(mark_RegionsTreated);
			m_mesh->freeMark<Region>(mark_RegionsTreated);

			// std::cout << "Regions dans le 1 ring: " << closest_node.get<Region>().size() << std::endl;
			// std::cout << "Regions dans le 2 ring: " << regions_to_test.size() << std::endl;

			for (auto r_id : regions_to_test)
			{
				Region r = m_mesh->get<Region>(r_id);
				if (!tetraFound)
				{
					tetraFound =
					   math::Utils::isInTetra(r.get<Node>()[0].point(), r.get<Node>()[1].point(), r.get<Node>()[2].point(), r.get<Node>()[3].point(), APoint);
					if (tetraFound)
					{
						tetra_id = r.id();
					}
				}
			}
		}

	}

	return tetra_id;
}
/*------------------------------------------------------------------------*/