/*----------------------------------------------------------------------------*/
#ifndef GMDS_SHEET_SELECTOR_3D_H
#define GMDS_SHEET_SELECTOR_3D_H
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
#include <map>
#include <vector>
/*----------------------------------------------------------------------------*/
namespace gmds{

    /*----------------------------------------------------------------------------*/
    /** @class  Selector3D
     *  @brief  Class that allows to select a sheet in a hex mesh
     *          The mesh we work on  must have R and N and the R2N fields
     */
    class EXPORT_GMDS Selector3D
    {
    public:

        /*------------------------------------------------------------------------*/
        /** @brief Selector3D.
         *
         *  @param AMesh the mesh we work on
         */
        Selector3D(Mesh* AMesh);

        /*------------------------------------------------------------------------*/
        /** @brief  Destructor.
         */
        virtual ~Selector3D();

        /*------------------------------------------------------------------------*/
        /** @brief Check if the mesh fits algorithm requirements, which are:
         *         - to be a 3D mesh,
         *         - to be a full hex mesh
         */
        bool isValid() const;
        /*------------------------------------------------------------------------*/
        /** @brief  Performs the sheet selection starting from edge [AN1,AN2]
         * @param[in] AN1 id of the first edge extremity
         * @param[in] AN2 id of the second edge extremity
         */
        void execute(const TCellID AN1, const TCellID AN2);
        /*------------------------------------------------------------------------*/
        /** @brief Gives access to the hexes stored during the last sheet retrieval
         * @param[out] AHexes the set of hexes that forms the sheet
         */
        void getSheetCells(std::vector<TCellID>& AHexes) const;
        /*------------------------------------------------------------------------*/
        /** @brief Gives access to the hexes stored during the last sheet retrieval
         * @return the set of hexes that forms the sheet
         */
        std::vector<TCellID> getSheetCells() const;
    private:

        typedef struct {
            TCellID hex_id; /// global id of a hex cell
            TInt e1;        /// local id of the 1st edge extremity
            TInt e2;        /// local id of the 2nd edge extremity
        } HexSheetInfo;
        /*------------------------------------------------------------------------*/
        /** @brief Fill the local connectivity N2R
         */
        void buildLocalN2R();

        /*------------------------------------------------------------------------*/
        /** @brief  Check if there exist an edge connecting AN1 and AN2
         * @param[in] AN1 id of the first edge extremity
         * @param[in] AN2 id of the second edge extremity
         */
        bool isAnEdge(const TCellID AN1, const TCellID AN2);

        /*------------------------------------------------------------------------*/
        /** @brief  get the set of hexes that are both adjacent to AN1 and AN2
         * @param[in] AN1 A first node id
         * @param[in] AN2 A second node id
         */
        std::vector<TCellID> getHexesSharingEdge(const TCellID AN1, const TCellID AN2);

        int getLocalEdgeIndex(const HexSheetInfo& AInfo);
    private:
        /** a mesh */
        Mesh* m_mesh;
        /** hexes in the sheet*/
        std::vector<TCellID> m_sheet_hexes;
        /** Inverse connectivity N23 build on purpose for this algorithm*/
        std::map<TCellID, std::vector<TCellID> > m_N2R;

    };
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif //GMDS_SHEET_SELECTOR_3D_H
