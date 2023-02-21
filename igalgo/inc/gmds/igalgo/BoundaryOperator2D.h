/*----------------------------------------------------------------------------*/
#ifndef GMDS_BOUNDARY_OPERATOR_2D_H_
#define GMDS_BOUNDARY_OPERATOR_2D_H_
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include "GMDSIgAlgo_export.h"
/*----------------------------------------------------------------------------*/
#include <vector>
#include <map>
/*----------------------------------------------------------------------------*/
namespace gmds{

    /*----------------------------------------------------------------------------*/
    /** @class  BoundaryOperator2D
     *  @brief  Class gathering operations used to get and mark the cells that
     *          belongs to the boundary of a 2D mesh. This class uses a 2D mesh
     *          having F2N and N2F connections. If edges and E2N is available,
     *          some extra operations are possible.
     */
    class GMDSIgAlgo_API BoundaryOperator2D
    {
    public:

        /*------------------------------------------------------------------------*/
        /** @brief Constructor.
         *
         *  @param AMesh the mesh where sheet operations are performed.
         */
        BoundaryOperator2D(Mesh* AMesh);

        /*------------------------------------------------------------------------*/
        /** @brief  Destructor.	*/
        virtual ~BoundaryOperator2D();


        bool isValid() const;

        /*------------------------------------------------------------------------*/
        /** @brief  Mark all the boundary cells.
         *			Marks are created inside the method, but must be free by the
         *          caller.
         *
         * 	@param AMarkEOnCurv edges on boundary curves are marked with it
         * 	@param AMarkNOnCurv nodes on boundary curves are marked with it
         * 	@param AMarkNOnPnt  nodes on geometric point are marked with it
         * 	@param AMarkAN isolated nodes (connected to nothing) are marked with it
         */

        void markCellOnGeometry(int AMarkEOnCurve,
                                int AMarkNOnCurve,
                                int AMarkNOnPnt,
                                int AMarkAN);

        /*------------------------------------------------------------------------*/
        /** @brief  Color the boundary edges with different colors to get one
         *          color per geometric curve
         *
         * @param[in] AMarkEOnCurv all the boundary edges are marked with it
         * @param[in] AMarkNOnPnt  all the nodes classified on points are marked
         *            with it
         * @param[in] AColor color variable to be used. Otherwise, variable named
         *            "BND_CURVE_COLOR" will be retrieved or created.
         */
        void colorEdges(int AMarkEOnCurv, int AMarkNOnPnt,
                        Variable<int>* AColor=nullptr);

        /*------------------------------------------------------------------------*/
        /** @brief  Color the boundary  nodes to get one color per geometric vertex
         *
         * @param[in] AMarkNOnPnt  all the nodes classified on points are marked
         *            with it
         * @param[in] AColor color variable to be used. Otherwise, variable named
         *            "BND_VERTEX_COLOR" will be retrieved or created.
         */
        void colorNodes(int AMarkNOnPnt,
                        Variable<int>* AColor=nullptr);

        /*------------------------------------------------------------------------*/
        /** @brief Return the boundary nodes
         *
         * @param  ANodeIDs the ids of boundary nodes
         */
        void getBoundaryNodes(std::vector<TCellID>& ANodeIDs);


        /*------------------------------------------------------------------------*/
        /** @brief Mark the boundary edges that fit geometric curves for surface
         *         mesh
         * @param[in]  AMarkCE edges on geom curves
         * @param[in]  AMarkCN nodes on geom curves
         */
        void markCellsOnCurves( int AMarkCE,  int AMarkCN);


        /*------------------------------------------------------------------------*/
        /** @brief  Mark the boundary nodes that fit geometric vertices
         * @param[in]  AMarkCE edges on geom curves
         * @param[in]  AMarkCN nodes on geom curves
         * @param[out]  AMarkPN nodes on geom points
         */
        void markNodesOnPoint( int AMarkCE, int AMarkCN,
                               int AMarkPN);

        /*------------------------------------------------------------------------*/
        /** @brief  Mark the boundary nodes that do not have adjacent edges
         * @param  AMarkAlone  mark of alone nodes
         */
        void markAloneNodes( int AMarkAlone);

    protected:
        /* a mesh */
        Mesh* m_mesh;
    };
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_BOUNDARY_OPERATOR_2D_H_ */
/*----------------------------------------------------------------------------*/
