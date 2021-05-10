/*----------------------------------------------------------------------------*/
#ifndef GMDS_SHEET_COLLAPSE_3D_H
#define GMDS_SHEET_COLLAPSE_3D_H
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
#include <map>
#include <vector>
/*----------------------------------------------------------------------------*/
namespace gmds{

    /*----------------------------------------------------------------------------*/
    /** @class  Collapse2D
     *  @brief  Class that allows to collapse/remove a complete sheet from a hex
     *          mesh
     *
     *          The mesh we work on  must have R and N and the R2N fields
     */
    class EXPORT_GMDS Collapse3D
    {
    public:

        /*------------------------------------------------------------------------*/
        /** @brief Constructor.
         *
         *  @param AMesh the mesh we work on
         */
        Collapse3D(Mesh* AMesh);

        /*------------------------------------------------------------------------*/
        /** @brief  Destructor.
         */
        virtual ~Collapse3D();

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
    private:

        /*------------------------------------------------------------------------*/
        /** @brief Fill the local connectivity N2R
         */
        void buildLocalN2R();

    private:
        /** a mesh */
        Mesh* m_mesh;
        /** hexes in the sheet*/
        std::vector<TCellID> m_sheet_hexes;
        /** Inverse connectivity N2R build on purpose for this algorithm*/
        std::map<TCellID, std::vector<TCellID> > m_N2R;

    };
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif //GMDS_SHEET_COLLAPSE_3D_H
