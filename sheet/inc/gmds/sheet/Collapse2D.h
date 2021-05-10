/*----------------------------------------------------------------------------*/
#ifndef GMDS_SHEET_COLLAPSE_2D_H
#define GMDS_SHEET_COLLAPSE_2D_H
/*----------------------------------------------------------------------------*/
#include <gmds/sheet/Operator2D.h>
/*----------------------------------------------------------------------------*/
#include <map>
#include <vector>
/*----------------------------------------------------------------------------*/
namespace gmds{


    namespace cad {
        class GeomMeshLinker;
    }
    /*----------------------------------------------------------------------------*/
    /** @class  Collapse2D
     *  @brief  Class that allows to collapse/remove a complete sheet from a quad
     *          mesh
     *
     *          The mesh we work on  must have F and N and the F2N fields and not R
     */
    class EXPORT_GMDS Collapse2D:public Operator2D
    {
    public:

        /*------------------------------------------------------------------------*/
        /** @brief Constructor.
         *
         *  @param AMesh the mesh we work on
         *  @param ALinker linker to the geom model
         */
        Collapse2D(Mesh* AMesh,  cad::GeomMeshLinker* ALinker);

        /*------------------------------------------------------------------------*/
        /** @brief  Destructor.
         */
        virtual ~Collapse2D();
        /*------------------------------------------------------------------------*/
        /** @brief  Performs the sheet selection starting from edge [AN1,AN2]
         * @param[in] AN1 id of the first edge extremity
         * @param[in] AN2 id of the second edge extremity
         */
        void execute(const TCellID AN1, const TCellID AN2);

    private:
        /**
         *
         * @param[in] AEdges     set of edges traversed by a sheet
         * @param[in] ASideColor mapping that indicate the side number of each
         *                       end point of edges in @p AEdges
         * @return true if the sheet that traverse edges of @p AEdges can be
         *         collapsed, otherwise false if a geometric classification issue
         *         is encountered
         */
        bool checkGeometricClassification(const std::vector<VirtualEdge>& AEdges,
                                          const std::map<TCellID,int>& ASideColor);

        void collapseAndReconnect(const std::map<TCellID , std::vector<TCellID> >& AToCollapse,
                                  const std::vector<TCellID>& ASheetQuads,
                                  const int AMarkSheetQuads);
    private:


        cad::GeomMeshLinker* m_geom_linker;
        /** quads in the sheet*/
        std::vector<TCellID> m_sheet_quads;
    };
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif //GMDS_SHEET_COLLAPSE_H
