/*----------------------------------------------------------------------------*/
#ifndef GMDS_SHEET_SELECTOR_2D_H
#define GMDS_SHEET_SELECTOR_2D_H
/*----------------------------------------------------------------------------*/
#include <gmds/sheet/Operator2D.h>
#include "GMDSSheet_export.h"
/*----------------------------------------------------------------------------*/
#include <map>
#include <vector>
/*----------------------------------------------------------------------------*/
namespace gmds{

    /*----------------------------------------------------------------------------*/
    /** @class  Selector2D
     *  @brief  Class that allows to select a sheet in a quad mesh
     *          The mesh we work on  must have F and N and the F2N fields and no
     *          regions
     */
    class GMDSSheet_API Selector2D: public Operator2D
    {
    public:

        /*------------------------------------------------------------------------*/
        /** @brief Selector2D.
         *
         *  @param AMesh the mesh we work on
         *  @param ALinker linker to the geom model
         */
        Selector2D(Mesh* AMesh);

        /*------------------------------------------------------------------------*/
        /** @brief  Destructor.
         */
        virtual ~Selector2D();

        /*------------------------------------------------------------------------*/
        /** @brief  Performs the sheet selection starting from edge [AN1,AN2]
         * @param[in] AN1 id of the first edge extremity
         * @param[in] AN2 id of the second edge extremity
         */
        void execute(const TCellID AN1, const TCellID AN2);
        /*------------------------------------------------------------------------*/
        /** @brief Gives access to the quads stored during the last sheet retrieval
         * @param[out] ACells the set of cells that forms the sheet
         */
        void getSheetCells(std::vector<TCellID>& ACells) const;
        /*------------------------------------------------------------------------*/
        /** @brief Gives access to the cells stored during the last sheet retrieval
         * @return the set of cells that forms the sheet
         */
        std::vector<TCellID> getSheetCells() const;
        /*------------------------------------------------------------------------*/
        /** @brief Gives the set of edges traversed by the sheet
         * @return the set of cells that forms the sheet
         */
        std::vector<VirtualEdge> getSheetTraversedEdges() const;
    private:

        typedef struct {
            TCellID quad_id; /// global id of a quad cell
            TInt e1;        /// local id of the 1st edge extremity
            TInt e2;        /// local id of the 2nd edge extremity
        } QuadSheetInfo;


        int getLocalEdgeIndex(const QuadSheetInfo& AInfo);
    private:


        /** quad info in the sheet*/
        std::vector<QuadSheetInfo> m_sheet_quad_infos;

    };
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif //GMDS_SHEET_SELECTOR_3D_H
