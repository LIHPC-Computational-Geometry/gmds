//
// Created by rochec on 14/04/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AeroException.h>
#include <gmds/claire/AeroExtrusion_2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

AeroExtrusion_2D::AeroExtrusion_2D(Mesh *AMesh) {
	m_mesh = AMesh;
}


/*------------------------------------------------------------------------*/
AeroExtrusion_2D::STATUS
AeroExtrusion_2D::execute()
{
	if(m_mesh==NULL)
		throw AeroException("ERROR: Invalid mesh pointer");

	return AeroExtrusion_2D::SUCCESS;
}
/*------------------------------------------------------------------------*/
