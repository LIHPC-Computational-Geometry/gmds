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
  m_mesh(AMesh)
{
	m_markNodesParoi = m_mesh->newMark<gmds::Node>();
	m_markNodesAmont = m_mesh->newMark<gmds::Node>();
}
/*------------------------------------------------------------------------*/