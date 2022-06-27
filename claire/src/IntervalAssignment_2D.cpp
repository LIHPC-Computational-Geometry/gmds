//
// Created by rochec on 27/06/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/IntervalAssignment_2D.h>
#include <gmds/ig/Mesh.h>
#include <gmds/claire/AeroExtrusion_2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

IntervalAssignment_2D::IntervalAssignment_2D(Blocking2D* ABlocking2D, ParamsAero Aparams_aero) {
	m_blocking = ABlocking2D;
	m_params_aero = Aparams_aero;
}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
IntervalAssignment_2D::STATUS
IntervalAssignment_2D::execute()
{

	return IntervalAssignment_2D::SUCCESS;
}
/*------------------------------------------------------------------------*/


/*-------------------------------------------------------------------*/
std::map<int, std::vector<TCellID>>
IntervalAssignment_2D::ComputeChords(){

	std::map<int, std::vector<TCellID>> map_chords;

	return map_chords;

}
/*-------------------------------------------------------------------*/