/*----------------------------------------------------------------------------*/
#include <gmds/igalgo/GridBuilder.h>
/*----------------------------------------------------------------------------*/
#include <sstream>
#include <set>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
GridBuilder::GridBuilder(Mesh* AMesh, const TInt ADim)
        :m_mesh(AMesh), m_dim(ADim)
{}
/*----------------------------------------------------------------------------*/
GridBuilder::~GridBuilder()
{}
/*----------------------------------------------------------------------------*/
bool GridBuilder::isValid() const
{
    if(m_dim==3)
        return (m_mesh->getModel()==(DIM3|R|N|R2N));
    else if(m_dim==2)
        return (m_mesh->getModel()==(DIM3|F|N|F2N) ||
                m_mesh->getModel()==(DIM2|F|N|F2N));

    //dimension error
    return false;
}
/*----------------------------------------------------------------------------*/
void GridBuilder::execute(const gmds::TInt AXNb, const gmds::TCoord AXStep, const gmds::TInt AYNb,
                          const gmds::TCoord AYStep, const gmds::TInt AZNb, const gmds::TCoord AZStep) {
    m_mesh->clear();
    if (m_dim == 2){
        build2D(AXNb, AXStep, AYNb, AYStep);
    }
    else if(m_dim==3){
        build3D(AXNb,AXStep,AYNb,AYStep,AZNb,AZStep);
    }
}
/*----------------------------------------------------------------------------*/
void GridBuilder::build2D(const gmds::TInt AXNb,
                          const gmds::TCoord AXStep,
                          const gmds::TInt AYNb,
                          const gmds::TCoord AYStep)
{
	std::vector<TCellID> node_ids = std::vector<TCellID>();
	const gmds::TInt N = AXNb * AYNb;
	node_ids.reserve(N);

	for (auto x = 0; x < AXNb; x++) {
		for (auto y = 0; y < AYNb; y++) {
			Node n = m_mesh->newNode(x*AXStep, y*AYStep, 0);
			node_ids.push_back(n.id());
		}
	}
	for (auto k = 0; k < N - AYNb; k++) {
		if ((k + 1) % AYNb == 0) {
			continue;
		}

		m_mesh->newQuad(node_ids.at(k),
			node_ids.at(k + AYNb),			// [x + 1][y],
			node_ids.at(k + AYNb + 1),		// [x + 1][y + 1],
			node_ids.at(k + 1));			// [x][y + 1],
	}
}
/*----------------------------------------------------------------------------*/
void GridBuilder::build3D(const gmds::TInt AXNb,
	const gmds::TCoord AXStep,
	const gmds::TInt AYNb,
	const gmds::TCoord AYStep,
	const gmds::TInt AZNb,
	const gmds::TCoord AZStep)
{
	std::vector<TCellID> node_ids = std::vector<TCellID>();
	const gmds::TInt N = AXNb * AYNb * AZNb;
	node_ids.reserve(N);

	for (auto x = 0; x < AXNb; x++) {
		for (auto y = 0; y < AYNb; y++) {
			for (auto z = 0; z < AZNb; z++) {
				Node n = m_mesh->newNode(x*AXStep, y*AYStep, z*AZStep);
				node_ids.push_back(n.id());
			}
		}
	}
	gmds::TInt AZYNb = AYNb * AZNb;
	for (auto k = 0; k < N - AZYNb - AZNb; k++) {
		if (((k + 1) % (AZNb) == 0) || (k%AZYNb >= AZYNb - AZNb)) {
			continue;
		}
		m_mesh->newHex(node_ids.at(k),
			node_ids.at(k + AZYNb),				// [x + 1][y][z],
			node_ids.at(k + AZYNb + AZNb),		// [x + 1][y + 1][z],
			node_ids.at(k + AZNb),				// [x][y + 1][z],
			node_ids.at(k + 1),					// [x][y][z + 1],
			node_ids.at(k + AZYNb + 1),			// [x + 1][y][z + 1],
			node_ids.at(k + AZYNb + AZNb + 1),	// [x + 1][y + 1][z + 1],
			node_ids.at(k + AZNb + 1));			// [x][y + 1][z + 1]);
	}
}
