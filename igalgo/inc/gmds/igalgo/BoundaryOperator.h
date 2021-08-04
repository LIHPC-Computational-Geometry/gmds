/*----------------------------------------------------------------------------*/
#ifndef GMDS_BOUNDARY_OPERATOR_H_
#define GMDS_BOUNDARY_OPERATOR_H_
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
#include <vector>
#include <map>
/*----------------------------------------------------------------------------*/
namespace gmds{

    /*----------------------------------------------------------------------------*/
    /** @class  BoundaryOperator
     *  @brief  Class gathering operations used to get and mark the cells that
     belongs to the boundary of a 2D, 3D mesh
     */
    class EXPORT_GMDS BoundaryOperator
    {
    public:

        /*------------------------------------------------------------------------*/
        /** @brief Constructor.
         *
         *  @param AMesh the mesh where sheet operations are performed.
         */
        BoundaryOperator(Mesh* AMesh,const double AAngle=0.72);

        /*------------------------------------------------------------------------*/
        /** @brief  Destructor.	*/
        virtual ~BoundaryOperator();

        /*------------------------------------------------------------------------*/
        /** @brief  Setter to the angle value used to determine if an edge is sharp
         *          enough to separate two surfaces*/
        void setSurfaceAngleDot(const double AD);

        /*------------------------------------------------------------------------*/
        /** @brief  Getter to the angle value used to determine if an edge is sharp
         *          enough to separate two surfaces*/
        double getSurfaceAngleDot();


        bool isValid() const;
        /*------------------------------------------------------------------------*/
        /** @brief  Mark all the boundary cells.
         *			Marks are created inside the method, but must be free by the
         *          caller. This function also initialize and fill in the variable
         *          BND_SURFACE_COLOR, which is assigned to faces.
         *
         * 	@param AMarkFOnSurf faces on boundary surfaces are marked with it
         * 	@param AMarkEOnSurf edges on boundary surfaces are marked with it
         * 	@param AMarkNOnSurf nodes on boundary surfaces are marked with it
         * 	@param AMarkEOnCurv edges on boundary curves are marked with it
         * 	@param AMarkNOnCurv nodes on boundary curves are marked with it
         * 	@param AMarkNOnPnt  nodes on geometric point are marked with it
         * 	@param AMarkAN isolated nodes (connected to nothing) are marked with it
         */

        void markCellOnGeometry(const int AMarkFOnSurf,
                                const int AMarkEOnSurf,
                                const int AMarkNOnSurf,
                                const int AMarkEOnCurve,
                                const int AMarkNOnCurve,
                                const int AMarkNOnPnt,
                                const int AMarkAN);

        /*------------------------------------------------------------------------*/
        /** @brief  Mark cells on the boundary surface
         *
         * 	@param AMarkBF faces on the boundary surfs are marked with it
         * 	@param AMarkBE edges on the boundary surfs are marked with it
         * 	@param AMarkBE edges on the boundary surfs are marked with it
         */

        void markCellsOnSurfaces(const int AMarkBF, const int AMarkBE, const int AMarkBN);

        /*------------------------------------------------------------------------*/
        /** @brief  Mark the boundary edges that fit geometric curves
         * @param  AMarkBF IN faces on geom surf
         * @param  AMarkBE IN edges on geom surf
         * @param  AMarkCE OUT edges on geom curves
         * @param  AMarkCN OUT nodes on geom curves
         */
        void markCellsOnCurves(const int AMarkBF , const int AMarkBE,
                               const int AMarkCE, const int AMarkCN);

        /*------------------------------------------------------------------------*/
        /** @brief Mark the boundary edges that fit geometric curves for surface
         *         mesh
         * @param  AMarkCE OUT edges on geom curves
         * @param  AMarkCN OUT nodes on geom curves
         */
        void markCellsOnCurves(const int AMarkCE, const int AMarkCN);


        /*------------------------------------------------------------------------*/
        /** @brief  Mark the boundary nodes that fit geometric vertices
         * @param  AMarkCE IN  edges on geom curves
         * @param  AMarkCN IN  nodes on geom curves
         * @param  AMarkPN OUT nodes on geom points
         */
        void markNodesOnPoint(const int AMarkCE, const int AMarkCN,
                              const int AMarkPN);

        /*------------------------------------------------------------------------*/
        /** @brief  Mark the boundary nodes that do not have adjacent edges
         * @param  AMarkAlone IN mark of alone nodes
         */
        void markAloneNodes(const int AMarkAlone);

        /*------------------------------------------------------------------------*/
        /** @brief  Color the boundary faces with different colors to get one color
         *          per geometric surface
         *
         * @param[in] AMarkFOnSurf all the boundary faces are marked with it
         * @param[in] AMarkEOnCurv all the edges classified on curves are marked with
         *            it
         * @param[in] AColor color variable to be used. Otherwise, variable named
         *            "BND_SURFACE_COLOR" will be retrieved or created.
         */
        void colorFaces(const int AMarkFOnSurf, const int AMarkEOnCurv,
                        Variable<int>* AColor=NULL);

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
        void colorEdges(const int AMarkEOnCurv, const int AMarkNOnPnt,
                        Variable<int>* AColor=NULL);

        /*------------------------------------------------------------------------*/
        /** @brief  Color the boundary  nodes to get one color per geometric vertex
         *
         * @param[in] AMarkNOnPnt  all the nodes classified on points are marked
         *            with it
         * @param[in] AColor color variable to be used. Otherwise, variable named
         *            "BND_VERTEX_COLOR" will be retrieved or created.
         */
        void colorNodes(const int AMarkNOnPnt,
                        Variable<int>* AColor=NULL);
        /*------------------------------------------------------------------------*/
        /** @brief Return the boundary nodes
         *
         * @param  ANodeIDs the ids of boundary nodes
         */
        void getBoundaryNodes(std::vector<TCellID>& ANodeIDs);

        /*------------------------------------------------------------------------*/
        /** @brief Compute the normal to AFace going out of ARegion
         * @param  AFace a face
         * @param  ARegion a region adj. to AFace
         * @return the normal to AFace going out of ARegion
         */
        math::Vector3d getOutputNormal(Face& AFace, Region& ARegion);

        /*------------------------------------------------------------------------*/
        /** @brief Compute the normal to AFace going out of a mesh
         * @param  AFace a face
         * @return the normal to AFace going out
         */
        math::Vector3d getOutputNormalOfABoundaryFace(Face& AFace);
        /*------------------------------------------------------------------------*/
        /** @brief Compute the normal to ANode as the average of the normal to the
         *         adjacent boundary faces (weighted by the face area)
         * @param  ANode a face
         * @return the normal to AFace going out
         */
        math::Vector3d getOutputNormalOfABoundaryNode(const Node& ANode);

    protected:
        /* a mesh */
        Mesh* m_mesh;

        double m_surface_angle_dot;
    };
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_BOUNDARY_OPERATOR_H_ */
/*----------------------------------------------------------------------------*/
