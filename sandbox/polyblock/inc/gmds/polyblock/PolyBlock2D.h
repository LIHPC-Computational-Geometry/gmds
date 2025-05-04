/*----------------------------------------------------------------------------*/
#ifndef GMDS_POLYBLOCK2D_H
#define GMDS_POLYBLOCK2D_H
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/math/Cross2D.h>
#include <Eigen/Sparse>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
/** @class PolyBlock2D
 *         This class provides a simple implementation for 2D polycube-based
 *         blocking
 */
    class PolyBlock2D {
    public:
        /**@brief Basic constructor
         *
         * @param AMesh 2D mesh we work on
         * @param ABndNodeMark node mark on @p AMesh that indicates which nodes
         *                     are on the boundary
         * @param ABndEdgeMark edge mark on @p AMesh that indicates which edges
         *                     are on the boundary
         */
        PolyBlock2D(Mesh* AMesh, const int ABndNodeMark, const int ABndEdgeMark);

        virtual ~PolyBlock2D();
        void execute();
        static const math::Vector3d REF[6];

    private:
        /**@brief update the content of the boundary normal variables.
         */
        void computeBndNormals();
        void computeTargetXYZ();
        void computeCrossField();
        void solveCrossField(const double AAlpha);

        /** @brief udpate the cross field values for faces adjacent to
         * the boundary to be aligned with the m_bnd_normal_edge value
         *
         */
        void updateCrossFieldBoundary();

    private:
        Mesh* m_mesh;
        int m_bnd_node_mark;
        int m_bnd_edge_mark;
        int m_bnd_face_mark;
        Variable<math::Vector3d>* m_bnd_normal_node;
        Variable<math::Vector3d>* m_bnd_normal_edge;
        Variable<math::Vector3d>* m_target_XYZ_node;
        Variable<math::Vector3d>* m_target_XYZ_edge;
        Variable<math::Cross2D>* m_cross_face;

        Variable<math::Vector3d>* m_v1_face;
        Variable<math::Vector3d>* m_v2_face;

        Variable<double>* m_x_value;
        Variable<double>* m_y_value;
        Variable<int>* m_patch_number;

        Variable<double>* m_U;
        Variable<double>* m_V;
        void writeDebugFile();

        void computeUVParam();

        void solveUV();

        void computeBoundaryPatchs();

        /**
         *
         * @param AEdgeID
         * @return 0 for an inner edge, 1 for a bnd U edge, 2 for a bnd V edge
         *
         * A U edge is supposed to be constant in U along the edge
         */
        int computeUVBndEdge(const TCellID AEdgeID);

        std::vector<Edge> getAdjBndEdges(const TCellID AEdgeID);

        void deformOmega();

        void buildB(Eigen::VectorXd &AB);

        int getMaxNodeValence() const;

        void buildA(Eigen::SparseMatrix<double> &AM);

        void projectInUV(const TCellID &AFaceID, math::Point &AZ1, math::Point& AZ2, math::Point& AZ3);
    };
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_POLYBLOCK2D_H
/*----------------------------------------------------------------------------*/
