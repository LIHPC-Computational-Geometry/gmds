/*----------------------------------------------------------------------------*/
#ifndef FRAME_FIELD_LAPLACE_SMOOTHER_2D_H_
#define FRAME_FIELD_LAPLACE_SMOOTHER_2D_H_
/*---------------------------------------------------------------------------*/
#include <map>
/*---------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/math/Quaternion.h>
#include <gmds/math/Matrix.h>
/*---------------------------------------------------------------------------*/
namespace  gmds {
/*---------------------------------------------------------------------------*/
    class FrameFieldLaplacianSmoothing {

    public:

        // the algorithm i can be applied on a mesh where Quaternions are defined
        // for some nodes and outward normal is given for nodes on a boundary
        // surface
        FrameFieldLaplacianSmoothing(Mesh *AMesh,
                                     Variable<math::Quaternion> *&AField,
                                     Variable<math::Vector3d> *&ANormal);

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

    private:
        void smooth();

        void smooth(Node &ANode);

    private:

        Mesh *m_mesh;
        Variable<math::Quaternion> *m_frame_field;
        Variable<math::Vector3d> *m_surf_normal;
        //mark to identify the nodes to work on
        int m_mark_candidates;
        // boolean indicating if the candidate mark has been defined
        bool m_candidate_mark_initialized;
        // mark to identify the nodes that are classified on points
        int m_mark_point;
        // mark to identify the nodes that are classified on curves
        int m_mark_curve;
        // mark to identify the nodes that are classified on surfaces
        int m_mark_surface;
        // boolean indicating if the candidate mark has been defined
        bool m_boundary_marks_initialized;
        //we keep in mind the number of call to HLBFGS
        int nb_HLBFGS_calls;

        //the ordered list of candidate nodes
        std::vector<Node> m_candidates;
        //the index of a candidate node in m_candidates
        std::map<TCellID, int> m_candidates_index;

        std::vector<std::vector<Node> > m_candidates_neighborhood;
    };
}
/*---------------------------------------------------------------------------*/
#endif //FRAME_FIELD_LAPLACE_SMOOTHER_2D_H_
/*---------------------------------------------------------------------------*/
