//
// Created by rochec on 26/01/2022.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/PointFollowingVectorField2D.h>
#include <limits>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/

PointFollowingVectorField2D::PointFollowingVectorField2D(Mesh *AMesh) {
	m_mesh = AMesh;
}


/*------------------------------------------------------------------------*/
PointFollowingVectorField2D::STATUS PointFollowingVectorField2D::execute()
{

	return PointFollowingVectorField2D::SUCCESS;
}
/*------------------------------------------------------------------------*/