//
// Created by rochec on 23/03/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AbstractAeroBoundaries.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
AbstractAeroBoundaries::AbstractAeroBoundaries(Mesh *AMesh) :
  m_mesh(AMesh),
  m_isImmerged(true)
{
	m_markNodesParoi = m_mesh->newMark<gmds::Node>();
	m_markNodesAmont = m_mesh->newMark<gmds::Node>();
	m_markBoundaryNodes = m_mesh->newMark<Node>();
}
/*------------------------------------------------------------------------*/

AbstractAeroBoundaries::~AbstractAeroBoundaries()
{
	m_mesh->unmarkAll<Node>(m_markBoundaryNodes);
	m_mesh->freeMark<Node>(m_markBoundaryNodes);
	m_mesh->unmarkAll<Node>(m_markNodesParoi);
	m_mesh->freeMark<Node>(m_markNodesParoi);
	m_mesh->unmarkAll<Node>(m_markNodesAmont);
	m_mesh->freeMark<Node>(m_markNodesAmont);
}
/*------------------------------------------------------------------------*/