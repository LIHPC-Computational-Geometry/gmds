/*----------------------------------------------------------------------------*/
#ifndef SH_FIELD_SMOOTHER_H_
#define SH_FIELD_SMOOTHER_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <vector>
#include <map>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
#include <gmds/ig/Mesh.h>
#include <gmds/math/Chart.h>
#include <gmds/math/AxisAngleRotation.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
/** \class  FieldSmoother
 *  \brief  Starting from a frame field, this algorithm smooth the frame following
 *          the algorithm described in Lui's paper
 */
    class EXPORT_GMDS FieldSmoother{

            public:
            /*------------------------------------------------------------------------*/
            /** \brief Constructor.
             *
             * \param[in] AMesh              the mesh where we work on
             * \param[in] ARotField          rotation field we have to smooth.
             * \param[in] ANormal            normal for each boundary node
             * \param[in] AMarkNodeOnPoint   mark for nodes classified on points
             * \param[in] AMarkNodeOnCurve   mark for nodes classified on curves
             * \param[in] AMarkNodeOnSurface mark for nodes classified on surfaces
             */
            FieldSmoother(gmds::IGMesh* AMesh,
            gmds::Variable<gmds::math::AxisAngleRotation>* ARotField,
            std::map<gmds::TCellID, gmds::math::Vector3d>& ANormal,
            const int AMarkNodeOnPoint,
            const int AMarkNodeOnCurve,
            const int AMarkNodeOnSurface);


            /*------------------------------------------------------------------------*/
            /** \brief Function to be called for smoothing the frame field
             */
            void execute();



            /*------------------------------------------------------------------------*/
            /** \brief Select nodes marked \p AMark as the ones to work on.
             *
             * \param[in] AMark a mark number on nodes
             */
            void selectNodes(const int AMark);


            protected:

            void writeSolution();
            /*------------------------------------------------------------------------*/
            /** \brief Initialize the static list of cells the algorithm is going to
             *         work on.
             */
            void initCandidates();

            /*------------------------------------------------------------------------*/
            /** \brief Initialize the Euler angles \p AEuler from the AxisAngleRotation
             *         variable. \p AEuler must be already initialized before calling
             *         this method
             *
             * \param[in] AEuler a 1D-array of Euler angles (in a C-style).
             */
            void buildEulerAngles(double*& AEuler);
            /*------------------------------------------------------------------------*/
            /** \brief Rebuild the the AxisAngleRotation variable from the Euler angles
             *         stored in \p AEuler.
             *
             * \param[in] AEuler a 1D-array of Euler angles (in a C-style).
             */
            void rebuildAxisAngleRotations(double*& AEuler);

            /*------------------------------------------------------------------------*/
            /** \brief Static function used for calling HLFBGS library.
             *
             * \param[in]  AN     nb variables to handle
             *                    (3 euler angles per vertex * nb. vertices)
             * \param[in]  AX     AN-sized array containing Euler angles(3 per vertex)
             * \param[in]  APrevX prev x value (mandatory fro HLBFGS library)
             *                    unused in our case
             * \param[out] AFunc  pointer on the computed result function
             * \param[out] AGrad  grad(func)
             */

            static void evalF(int AN,
            double* AX, double *APrevX,
            double* AFunc, double* AGrad);

            /*------------------------------------------------------------------------*/
            /** \brief Static function used in evalF for calling HLFBGS library.
             *
             * \param[in]  ANode  the node we work on
             * \param[in]  ACos   reference on pre-calculated cosinus values
             * \param[in]  ASin   reference on pre-calculated sinus   values
             * \param[out] AMix   first  matrix for first Euler angle
             * \param[out] AMiy   second matrix for second Euler angle
             * \param[out] AMiz   third  matrix for third Euler angle
             * \param[out] ADMix  first  derivative matrix for first Euler angle
             * \param[out] ADMiy  second derivative matrix for second Euler angle
             * \param[out] ADMiz  third  derivative matrix for third Euler angle
             */
            static void initMatrix(gmds::Node& ANode,
            std::vector<gmds::TCoord>& ACos,
            std::vector<gmds::TCoord>& ASin,
            gmds::math::Matrix<3, 3, double>& AMix,
            gmds::math::Matrix<3, 3, double>& AMiy,
            gmds::math::Matrix<3, 3, double>& AMiz,
            gmds::math::Matrix<3, 3, double>& ADMix,
            gmds::math::Matrix<3, 3, double>& ADMiy,
            gmds::math::Matrix<3, 3, double>& ADMiz);

            /*------------------------------------------------------------------------*/
            /** \brief Static function used for calling HLFBGS library.
             *
             * \param[in]  ANbIter   the number of iterations
             * \param[in]  AIter     the current itration number
             * \param[in]  AX        the current unknown vector
             * \param[in]  AFunc     the current function value
             * \param[in]  AGrad     the current grad value
             * \param[in]  AGradNorm the current grad norm error
             */

            static void newiteration(int ANbIter,
            int AIter,
            double* AX,
            double* AFunc,
            double* AGrad,
            double* AGradNorm);

            protected:

            /** the mesh we work on */
            static gmds::IGMesh* m_mesh;

            /** rotation field previously computed for each node. */
            gmds::Variable<gmds::math::AxisAngleRotation>*    m_rotation_field;

            /** Defines a normal constraint alignment for boundary nodes
             * we work on */
            std::map<gmds::TCellID, gmds::math::Vector3d> m_normal;

            /* mark for boundary nodes classified on points */
            static int m_mark_node_on_point;
            /* mark for boundary nodes classified on curves */
            static int m_mark_node_on_curve;
            /* mark for boundary nodes classified on surfaces */
            static int m_mark_node_on_surface;
            //mark to identify the nodes to work on
            int m_mark_candidates;
            // boolean indicating if the candidate mark has been defined
            bool m_candidate_mark_initialized;

            //========================================================
            // STATIC VARIABLES FOR THE CALL TO HLBFGS
            //========================================================
            //the ordered list of candidate nodes
            static std::vector<gmds::Node> m_candidates;
            //the list of edges that connect candidates
            // these edges will be used to compute smoothing values between
            // end nodes
            static std::vector<gmds::Edge> m_edges;
            //the index of a candidate node in m_candidates
            static std::map<gmds::TCellID, int> m_candidates_index;
            //We keep in mind a chart for each boundary node we work on
            static std::map<gmds::TCellID, gmds::math::Chart> m_boundary_triad;

    };
}
/*----------------------------------------------------------------------------*/
#endif /* SH_FIELD_SMOOTHER_H_ */
/*----------------------------------------------------------------------------*/
