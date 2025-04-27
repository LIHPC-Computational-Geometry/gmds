/*----------------------------------------------------------------------------*/
#include <gmds/smoothy/LaplacianSmoother3UC.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::smoothy;
/*----------------------------------------------------------------------------*/
LaplacianSmoother3UC::LaplacianSmoother3UC(gmds::Mesh *AMesh) : AbstractSmoother(AMesh) {}
/*----------------------------------------------------------------------------*/
bool LaplacianSmoother3UC::isValid() const
{
	MeshModel model = m_mesh->getModel();

	if (!model.has(N2R) && !model.has(R2N)) return false;

	return true;
}
/*----------------------------------------------------------------------------*/
int LaplacianSmoother3UC::smooth()
{
	// First we build the N2N connectivity for the nodes to be smoothed
	std::map<TCellID, std::set<TCellID>> n2n;
	for (auto nid : m_nodes) {
		Node ni = m_mesh->get<Node>(nid);
		std::vector<Region> ni_regions = ni.get<Region>();
		for (const auto& r : ni_regions) {
			TCellID id1, id2, id3;
			getAdjacentNodes(r, nid, id1, id2, id3);
			n2n[nid].insert(id1);
			n2n[nid].insert(id2);
			n2n[nid].insert(id3);
		}
	}

	for (auto i = 0; i < m_nb_iterations; i++) {

		// We smooth each node in a naive way
		for (auto nid : m_nodes) {
			Node ni = m_mesh->get<Node>(nid);
			std::set<TCellID> adj_n = n2n[nid];
			math::Point pi(0, 0, 0);
			for (auto adj_id : adj_n) {
				math::Point pj = m_mesh->get<Node>(adj_id).point();
				pi = pi + pj;
			}
			pi = pi * (1.0 / adj_n.size());
			ni.setPoint(pi);
		}
	}

	return 1;
}

/*----------------------------------------------------------------------------*/
void LaplacianSmoother3UC::getAdjacentNodes(const Region& ARegion, const TCellID & ANode1, TCellID& ANode2, TCellID& ANode3, TCellID& ANode4){
	std::vector<TCellID> nodes = ARegion.getIDs<Node>();
		int index1 = -1, index2, index3, index4;
		for (int i = 0; i < 8; i++) {
			if (nodes[i] == ANode1)
				index1 = i;
		}
		switch (index1) {
		case -1:
#ifdef _DEBUG
			std::cout << "Index error" << std::endl;
#endif  //_DEBUG
			throw GMDSException("getAdjacentNodes: node 1 is not adjacent to the face");
			break;
		case 0:
			index2 = 3;
			index3 = 1;
		   index4 = 4;
			break;
		case 1:
			index2 = 0;
			index3 = 2;
		   index4 = 5;
			break;
		case 2:
			index2 = 1;
			index3 = 3;
		   index4 = 6;
			break;
		case 3:
			index2 = 2;
			index3 = 0;
		   index4 = 7;
			break;
	   case 4:
		   index2 = 7;
		   index3 = 5;
		   index4 = 0;
		   break;
	   case 5:
		   index2 = 4;
		   index3 = 6;
		   index4 = 1;
		   break;
	   case 6:
		   index2 = 5;
		   index3 = 7;
		   index4 = 2;
		   break;
	   case 7:
		   index2 = 6;
		   index3 = 4;
		   index4 = 3;
		   break;
		}

		ANode2 = nodes[index2];
		ANode3 = nodes[index3];
	   ANode4 = nodes[index4];
}