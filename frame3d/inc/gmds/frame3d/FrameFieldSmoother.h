/*---------------------------------------------------------------------------*/
#ifndef FRAME_FIELD_SMOOTHER_2D_H_
#define FRAME_FIELD_SMOOTHER_2D_H_
/*---------------------------------------------------------------------------*/
#include <iostream>
/*---------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/math/Quaternion.h>
#include <gmds/math/Matrix.h>
/*---------------------------------------------------------------------------*/
namespace  gmds {
/*---------------------------------------------------------------------------*/
	class FrameFieldSmoother {

	public:

		// the algorithm can be applied on a mesh where quaternions are defined
		// for some nodes and outward normal is given for nodes on a boundary
		// surface
		FrameFieldSmoother(gmds::Mesh *AMesh,
						   gmds::Variable<gmds::math::Quaternion> *&AField,
						   gmds::Variable<gmds::math::Vector3d> *&ANormal);

		/* All the node with AMark must be involved in the smoothing process
        */
		void selectNodes(const int AMark);

		// nodes classified on curves, surfaces and points must be marked
		// the option is to provide this mark or to recompute them
		// Only the first option is provided right now.
		void initBoundaryMarks(const int AMarkPnt, const int AMarkCurve,
							   const int AMarkSurf);

		void execute();

		void initCandidates();

		int initOptimizationData(double *&eulerAngles);

		void rebuildQuaternions(double *&eulerAngles);

	private:
		void performSmoothing();

	private:

		gmds::Mesh *m_mesh;
		gmds::Variable<gmds::math::Quaternion> *m_frame_field;
		gmds::Variable<gmds::math::Vector3d> *m_surf_normal;
		//mark to identify the nodes to work on
		int m_mark_candidates;
		// boolean indicating if the candidate mark has been defined
		bool m_candidate_mark_initialized;
		// boolean indicating if the candidate mark has been defined
		bool m_boundary_marks_initialized;
		//we keep in mind the number of call to HLBFGS
		int nb_HLBFGS_calls;

	};
}
/*---------------------------------------------------------------------------*/
#endif //FRAME_FIELD_SMOOTHER_2D_H_
/*---------------------------------------------------------------------------*/
