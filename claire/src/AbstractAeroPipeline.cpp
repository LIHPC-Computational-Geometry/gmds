//
// Created by rochec on 09/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AbstractAeroPipeline.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AbstractAeroPipeline::AbstractAeroPipeline(ParamsAero Aparams, gmds::Mesh && Am) :
  m_params(Aparams),
  m_m(std::move(Am))
{
	m_mesh = &m_m;
	m_isOver = false;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
bool AbstractAeroPipeline::getIsOver(){
	return m_isOver;
}
/*------------------------------------------------------------------------*/