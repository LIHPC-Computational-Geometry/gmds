/*----------------------------------------------------------------------------*/
#ifndef GMDS_PILLOW2D_H
#define GMDS_PILLOW2D_H
/*----------------------------------------------------------------------------*/
#include <gmds/sheet/Operator2D.h>
#include "GMDSSheet_export.h"
/*----------------------------------------------------------------------------*/
#include <map>
#include <vector>
/*----------------------------------------------------------------------------*/
namespace gmds{

    /*----------------------------------------------------------------------------*/
    /** @class  SheetPillow2D
     *  @brief  Class that provides a way to insert a full layer of quads around
     *          a set of quads, called a shrink set.
     *
     *          The mesh that is transformed must only have F and N and the F2N
     *          fields. Geometric classification is not handled yet. It must also
     *          be well-oriented (ce MeshDoctor::orient2D)
     *
     */
    class GMDSSheet_API Pillow2D : public Operator2D
    {
    public:

        /*------------------------------------------------------------------------*/
        /** @brief Constructor.
         *
         *  @param AMesh the mesh the pillowing operation will be performed on
         *  @param AEscapeBnd boolean value indicating the pillowing behaviour on the
         *                    boundary. True means that the pillow layer goes through
         *                    it. False means that the pillow layer goes along
         *                    boundary edges.
         */
        Pillow2D(Mesh* AMesh, const bool AEscapeBnd=true);

        /*------------------------------------------------------------------------*/
        /** @brief  Destructor.	*/
        virtual ~Pillow2D();

        /*------------------------------------------------------------------------*/
        /** @brief  Performs the pillowing
         * @param[in]  AFaceIDs the set of regions we want to pillow
         * @param[in]  AWithCheck if true check that the set of curves surrounding
         *          @AFaceIDs is valid
         *
         * @return true if the execution succeeded, false otherwise
         */
        bool execute(const std::vector<TCellID>& AFaceIDs, bool AWithCheck=false);
        /*------------------------------------------------------------------------*/
        /** @brief Performs the pillowing
         * @param  AEdge a set of virtual edges defined by nodes.
         * @param AWithCheck if true check that the set of edges \p AEdges is valid
         *
         * @return true if the execution succeeded, false otherwise
         */
        bool execute(std::vector<VirtualEdge>& AEdges, bool AWithCheck=false);
        /*------------------------------------------------------------------------*/
        /**@brief Setter of the boundary behaviour
         *
         * @param AEscape true means that the inserted sheet must execute through the
         *                boundary, while false means the sheet remains along the
         *                boundary curves.
         */
        void setBoundaryBehaviour(const bool AEscape);
        /*------------------------------------------------------------------------*/
        /**@brief Check if the curve encloses a shrink set
         * @param AEdges a set of edges defining a curve in the mesh
         * @return true if it encloses, false otherwise
         */
        bool checkCurve(std::vector<VirtualEdge>& AEdge);

    private:
        /*------------------------------------------------------------------------*/
        /**@brief   Orient edges of @AEdges in a consistent manner
         *
         *          Does not require the N2F connectivity, only works on the edges
         *
         * @param AEdges
         */
        void orient(std::vector<VirtualEdge>& AEdges);

        /*------------------------------------------------------------------------*/
        /** @brief Check if a node is on the mesh boundary
         * @param ANodeID the id of the node to be checked
         */
         bool isABoundaryNode(const TCellID ANodeID);
        /*------------------------------------------------------------------------*/
        /** @brief Check if an edge is on the mesh boundary
         * @param AN1 id of the first edge node
         * @param AN2 id of the second edge node
         */
        bool isABoundaryEdge(const TCellID AN1,const TCellID AN2);
        /*------------------------------------------------------------------------*/
        /** @brief Check if an edge is on the mesh boundary
         * @param AE A virtual edge
         */
        bool isABoundaryEdge(const VirtualEdge& AF);


        /*------------------------------------------------------------------------*/
        /** @brief Get the faces sharing an edge with face \p AID
         * @param AID a face id
         * @return the adjacent faces (1 or 2)
         */
        std::vector<TCellID> getFacesSharingAnEdge(const TCellID AID);


    private:
        /** options to indicate the behaviour of the operation along the boundary*/
        bool m_escape_bnd;
    };
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif //GMDS_PILLOW3D_H
