/*----------------------------------------------------------------------------*/
#include <gmds/smoothy/LaplacianSmoother2UC.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::smoothy;
using namespace gmds::cad;
/*----------------------------------------------------------------------------*/
LaplacianSmoother2UC::LaplacianSmoother2UC(gmds::Mesh *AMesh) : AbstractSmoother(AMesh) {}
/*----------------------------------------------------------------------------*/
bool LaplacianSmoother2UC::isValid() const
{
	MeshModel model = m_mesh->getModel();
	if (model.has(R)) return false;

	if (!model.has(N2F) && !model.has(F2N)) return false;

	return true;
}
/*----------------------------------------------------------------------------*/
int LaplacianSmoother2UC::smooth()
{
	// First we build the N2N connectivity for the nodes to be smoothed
	std::map<TCellID, std::set<TCellID>> n2n;
	for (auto nid : m_nodes) {
		Node ni = m_mesh->get<Node>(nid);
		std::vector<Face> ni_faces = ni.get<Face>();
		for (auto f : ni_faces) {
			TCellID id1, id2;
			f.getAdjacentNodes(nid, id1, id2);
			n2n[nid].insert(id1);
			n2n[nid].insert(id2);
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
