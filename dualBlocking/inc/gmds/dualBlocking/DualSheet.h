//
// Created by calderans on 16/12/19.
//

#ifndef GMDS_DUALSHEET_H
#define GMDS_DUALSHEET_H

#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Vector.h>
#include <map>
#include <gmds/ig/Mesh.h>

namespace db {
    class DualSheet {

    public:
        DualSheet(std::vector<gmds::TCellID> ASurface, int AID, gmds::Mesh* AMesh);
        DualSheet(int AID, gmds::Mesh* AMesh);
        DualSheet(const DualSheet &ASheet);
        virtual ~DualSheet();

        /*------------------------------------------------------------------------*/
        /** \brief Fill the surface vector
         */
        void setSurface(std::vector<gmds::TCellID> ASurface);
        void setSurface(std::map<gmds::TCellID,std::vector<gmds::math::Vector3d>> ASurface);

        /*------------------------------------------------------------------------*/
        /** \brief Get a vector that contains the id of the tet in the sheet
         */
        std::vector<gmds::TCellID> getSurface();

        /*------------------------------------------------------------------------*/
        /** \brief Get the id of the dual sheet
         */
        int getID();

        void boundarySheet();
        bool isBoundary();

        /*------------------------------------------------------------------------*/
        /** \brief Build the mesh used to visualize the surface in plan mode
         */
        void buildSurface(gmds::Mesh* AMesh);

        gmds::Mesh* getSurfaceMesh();

    protected:
        std::vector<gmds::TCellID> m_surface;
        gmds::Mesh* m_mesh_surface;
        int m_sheet_ID;
        bool boundary = false;

    };
}


#endif //GMDS_DUALSHEET_H
