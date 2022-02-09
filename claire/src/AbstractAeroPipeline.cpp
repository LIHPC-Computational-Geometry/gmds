//
// Created by rochec on 09/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AbstractAeroPipeline.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AbstractAeroPipeline::AbstractAeroPipeline() {
	m_markFrontNodesParoi = m_mesh->newMark<gmds::Node>();
	m_markFrontNodesExt = m_mesh->newMark<gmds::Node>();
}