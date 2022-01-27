//
// Created by rochec on 27/01/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/LeastSquaresGradientComputation2D.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/


LeastSquaresGradientComputation2D::LeastSquaresGradientComputation2D(Mesh *AMesh, Variable<double>* Adistance) {
	m_mesh = AMesh;
	m_distance = Adistance;
	m_gradient2D = m_mesh->newVariable<math::Vector3d ,GMDS_FACE>("gradient_2D");
}


/*------------------------------------------------------------------------*/
LeastSquaresGradientComputation2D::STATUS LeastSquaresGradientComputation2D::execute()
{

	for(auto face_id:m_mesh->faces()) {
		//math::Vector3d Gradient = computeGradientOnSimpleFace(face_id);
		//m_gradient2D->set(face_id, Gradient) ;
	}

	return LeastSquaresGradientComputation2D::SUCCESS;
}
/*------------------------------------------------------------------------*/