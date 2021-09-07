//
// Created by calderans on 07/03/19.
//

#ifndef GMDS_BOUNDARYSURFACECREATOR_H
#define GMDS_BOUNDARYSURFACECREATOR_H

#include "DualSheetCreator.h"

#include <gmds/math/Chart.h>
#include <Predicates_psm.h>

namespace db {
    class BoundarySurfaceCreator: public DualSheetCreator{

    public:
        /*------------------------------------------------------------------------*/
        /** \brief Constructor.
         *
         * \param[in] AMesh background tetrahedral mesh we start from
         */
        BoundarySurfaceCreator(gmds::Mesh* AMesh,
                //gmds::TCellID AID,
                //int ASheetID,
                               gmds::Variable<std::vector<DualSheetCreator::intersectInfo>>* info_Var,
                               gmds::Variable<gmds::math::Chart>* vertex_chart_Var,
                               gmds::Variable<int>* propagation_round_Var,
                               gmds::Variable<int>* singularities_Var,
                               int wave_precedent_mark,
                               int face_treated_mark,
                               int wave_tet_mark,
                               double AMin_length,
                               gmds::cad::GeomMeshLinker &ALinker,
                               gmds::cad::FACManager &AManager);

        /*------------------------------------------------------------------------*/
        /** \brief Destructor.
         */
        virtual ~BoundarySurfaceCreator();

        /*------------------------------------------------------------------------*/
        /** \brief Execute the algorithm of boundary surface creation starting in tetra tetID
         * and following the face on the geometric model.
         */
        bool execute();

        /*------------------------------------------------------------------------*/
        /** \brief Set the id of the geom surface link to the sheet
         */
        void setSurfaceID(int Aid);

        /*------------------------------------------------------------------------*/
        /** \brief Get the id of the geom surface link to the sheet
         */
        int getSurfaceID();

        gmds::Mesh* buildSurfaceSheet();

    protected:

        int surface_id;

    };
}

#endif //GMDS_BOUNDARYSURFACECREATOR_H
