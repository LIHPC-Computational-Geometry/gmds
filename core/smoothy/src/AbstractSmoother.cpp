/*----------------------------------------------------------------------------*/
#include <gmds/smoothy/AbstractSmoother.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::smoothy;
using namespace gmds::cad;
/*----------------------------------------------------------------------------*/
AbstractSmoother::AbstractSmoother(Mesh* AMesh)
:m_mesh(AMesh),m_nb_iterations(0)
{}
/*----------------------------------------------------------------------------*/
void
AbstractSmoother::setNbIterations(const int ANbIterations)
{
	if(ANbIterations<0)
		throw GMDSException("Smoothing iteration number must be positive!");
	m_nb_iterations = ANbIterations;
}
/*----------------------------------------------------------------------------*/
void
AbstractSmoother::setNodes(std::vector<TCellID> &ANodeIds)
{
	m_nodes = ANodeIds;
}
/*----------------------------------------------------------------------------*/
void
AbstractSmoother::setNodesAll()
{
	m_nodes.reserve(m_mesh->getNbNodes());
	for(auto n_id:m_mesh->nodes())
		m_nodes.push_back(n_id);
}
/*----------------------------------------------------------------------------*/
