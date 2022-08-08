//
// Created by rochec on 08/08/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/RefinementBetaBlocking.h>
#include <gmds/ig/Mesh.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

RefinementBetaBlocking::RefinementBetaBlocking(Blocking2D *ABlocking2D, ParamsAero Aparams_aero) {
	m_blocking = ABlocking2D;
	m_params_aero = Aparams_aero;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
RefinementBetaBlocking::STATUS
RefinementBetaBlocking::execute()
{


	return RefinementBetaBlocking::SUCCESS;
}
/*------------------------------------------------------------------------*/