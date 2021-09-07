//
// Created by calderans on 07/03/19.
//

#ifndef GMDS_DUALSURFACECREATOR_H
#define GMDS_DUALSURFACECREATOR_H

#include <gmds/utils/CommonTypes.h>
#include <gmds/math/Segment.h>
#include "DualSheetCreator.h"

namespace db {
    class DualSurfaceCreator : public DualSheetCreator{

    public:

        /*------------------------------------------------------------------------*/
        /** \brief Constructor.
         *
         * \param[in] AMesh background tetrahedral mesh we start from
         */
        DualSurfaceCreator(gmds::Mesh* AMesh,
                //gmds::TCellID AID,
                //int ASheetID,
                           gmds::Variable<std::vector<intersectInfo>>* info_Var,
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
        virtual ~DualSurfaceCreator();

        /*------------------------------------------------------------------------*/
        /** \brief Setting the starting normal vector of the dual surface.
         *
         * \param norm vector used as normal
         */
        void setNorm(const gmds::math::Vector3d &norm);

        /*------------------------------------------------------------------------*/
        /** \brief Get the direction vector of the sheet.
         *
         * \return A 3D vector.
         */
        gmds::math::Vector3d getNorm();

        /*------------------------------------------------------------------------*/
        /** \brief Execute the algorithm of boundary surface creation starting in tetra tetID
         * and following the face on the geometric model.
         */
        virtual bool execute();

        gmds::Mesh* buildSurfaceSheet();

        void operator=(const DualSurfaceCreator& ADS) {

            this->sheet_ID = ADS.sheet_ID;
            this->tetID = ADS.tetID;
            this->surface = ADS.surface;
            this->m_mesh = ADS.m_mesh;
            this->min_length = ADS.min_length;
            this->linker = ADS.linker;
            this->manager = ADS.manager;

            this->var_info = ADS.var_info;
            this->m_propagation_round = ADS.m_propagation_round;
            this->m_singularities = ADS.m_singularities;
            this->wave_precedent = ADS.wave_precedent;
            this->face_treated = ADS.face_treated;
            this->mark_wave_tet = ADS.mark_wave_tet;
            this->start_norm = ADS.start_norm;
        }
    protected:

        gmds::math::Vector3d start_norm;
    };
}


#endif //GMDS_DUALSURFACECREATOR_H
