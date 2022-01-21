//
// Created by rochec on 21/01/2022.
//
/*------------------------------------------------------------------------*/
#include <gmds/claire/GradientComputation3D.h>
#include <limits>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/


GradientComputation3D::GradientComputation3D(Mesh *AMesh, Variable<double>* Adistance) {
	m_mesh = AMesh;
	m_distance = Adistance;
	m_gradient3D = m_mesh->newVariable<math::Vector3d ,GMDS_FACE>("gradient_3D");
}


/*------------------------------------------------------------------------*/
GradientComputation3D::STATUS GradientComputation3D::execute()
{

	for(auto face_id:m_mesh->faces()) {
		//math::Vector3d Gradient = computeGradientOnSimpleFace(face_id);
		//m_gradient3D->set(face_id, Gradient) ;
	}

	return GradientComputation3D::SUCCESS;
}
/*------------------------------------------------------------------------*/
