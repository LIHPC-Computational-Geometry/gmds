//
// Created by calderans on 06/03/19.
//

#ifndef GMDS_DUAL_H
#define GMDS_DUAL_H

#include <gmds/ig/Mesh.h>
#include "DualSheetCreator.h"
#include "DualSheet.h"

#include <gmds/dualBlocking/DualSurfaceCreator.h>
#include <gmds/dualBlocking/BoundarySurfaceCreator.h>


#include <gmds/cad/FACManager.h>
#include <gmds/cad/GeomMeshLinker.h>


namespace db {
    class DualBlockingSession {

    public:
        /*------------------------------------------------------------------------*/
        /** \brief Constructor.
         *
         * \param[in] AMesh background tetrahedral mesh we start from
         */
        DualBlockingSession(gmds::Mesh* AMesh, gmds::Mesh* AHMesh);

        virtual ~DualBlockingSession();

        void init();

        std::vector<gmds::TCellID> createSurface(gmds::TCellID AID, int AAxis, int &ASheetID, gmds::Mesh* ASurface_mesh);

        std::vector<gmds::TCellID> createBoundary(gmds::TCellID AID, int &ASheetID, gmds::Mesh* ASurface_mesh);

        bool checkValidity();

        void getSheetVisualisation(int ASheet);

        std::vector<gmds::TCellID> getSheet(int ASheet);

        /*------------------------------------------------------------------------*/
        /** \brief Color the tetrahedral mesh to create dual zone between dual surface
         */
        bool colorDual(int &ANbZones);

        /*------------------------------------------------------------------------*/
        /** \brief Classify the dual zones of id \AID :
         *          -Corner = 0
         *          -Edge = 1
         *          -Face = 2
         *          -Volume = 3
         */
        void classifyZone(int AID);

        /*------------------------------------------------------------------------*/
        /** \brief Refine the surface sheet to create a gap between the sheet and the geometric surface
         */
        void refineSurfaceSheet();

        void smoothing();

        /*------------------------------------------------------------------------*/
        /** \brief Create the hex blocks in the blocks structure mesh
         */
        int createBlock();

        void smoothBlocks();
        void resetBlockSmoothing(){;};
        void interpolation(int Ait);

        void interpolationV2(int Ait);

        /*------------------------------------------------------------------------*/
        /** \brief Remove the sheet of id AID
         */
        bool deleteSheet(int AID);

        int getSheetOfTet(gmds::TCellID AID);

        void getFrameAxis(double result[][3],gmds::TCellID AID);

        gmds::Mesh* getSurfaceMesh(int ASheetID);

        /*------------------------------------------------------------------------*/
        /** \brief Uncolor every tet and clean the block structure mesh
         */
        void resetDual();

        /*------------------------------------------------------------------------*/
        /** \brief Clean the block structure, deleting the regions of the mesh but not the vertices
         */
        void resetBlocks();

        /*------------------------------------------------------------------------*/
        /** \brief Get the singular graph of the tet mesh as a list of tet id
         */
        std::vector<gmds::TCellID> getSingularGraph();


    protected:
        /*------------------------------------------------------------------------*/
        /** \brief Color the tet AId in the closest X,Y,Z frame direction of the vector AV with the id of ASheet
         */
        int colorTet(gmds::TCellID AId, gmds::math::Vector3d AV, int ASheet);

        /*------------------------------------------------------------------------*/
        /** \brief Test if a tetra is in a sheet
         * \param AId is the tetra to test
         * \param ASheet is the id of the sheet
         * \return true if the tetra is in the sheet esle false
         */
        bool isColoredIn(gmds::TCellID AId, int ASheet);

        bool isColored(gmds::TCellID AId);

        /*------------------------------------------------------------------------*/
        /** \brief Compute the chart at the center of a tetra
         */
        gmds::math::Chart computeChart(std::vector<gmds::Node> ANodes, gmds::math::Point AP);

        /*------------------------------------------------------------------------*/
        /** \brief Find the closest vector of \p AV in the chart \AChart
         */
        const gmds::math::Vector3d closestComponentVector(const gmds::math::Vector3d& AV,
                                                          gmds::math::Chart& AChart);



        //TODO
        void checkContact();

        int getColor(gmds::TCellID AID, gmds::math::Vector3d AV);





        std::vector<int> createCorner(int AZone,std::tuple<int,int,int,int,int,int,int,int> ABlock, std::vector<int> ADual);

        int findOutFace(const gmds::math::Point&               AFromPnt,
                        const gmds::math::Vector3d&             AFromDir,
                        int                                     AFromTet,
                        const std::vector<gmds::Face>&          AFaces,
                        const std::vector<gmds::math::Triangle>&ATri,
                        gmds::math::Vector3d&                   AToDir);

        void classifyEdges();

        void uncolorTet(gmds::TCellID ATetID, int ASheetID);

        gmds::Mesh* getBlocks();


    protected:

        gmds::Mesh* m_mesh;

        gmds::Mesh* m_imesh;

        gmds::Mesh* m_hmesh;

        //std::vector<db::DualSurfaceCreator> dual_sheet;
        //std::vector<db::BoundarySurfaceCreator> dual_bounder;

        db::DualSurfaceCreator* m_surfaceCreator;
        db::BoundarySurfaceCreator* m_boundaryCreator;

        std::vector<db::DualSheet> dual_sheets;
        std::map<int, int> m_boundarySheets_to_surf;

        //A map used to store the dual zones of the mesh, <id,block node>
        std::map<int,gmds::Node> m_dual_zones;
        std::map<int,std::vector<int>> m_zone_border;

        gmds::Variable<gmds::math::Chart>* m_tetra_chart;
        gmds::Variable<gmds::math::Chart>* m_vertex_chart;

        gmds::Variable<std::vector<db::DualSheetCreator::intersectInfo>>* var_info;
        gmds::Variable<int>* m_propagation_round;


        gmds::Variable<int>* m_sing;

        gmds::Variable<int>* m_sheet_X;
        gmds::Variable<int>* m_sheet_Y;
        gmds::Variable<int>* m_sheet_Z;

        gmds::Variable<int>* m_block;

        gmds::cad::FACManager manager;
        gmds::cad::GeomMeshLinker linkerG_T;
        gmds::cad::GeomMeshLinker linkerH_G;

        gmds::Variable<int>* m_block_id;

        int mark_ghost;


        int mark_visited;

        int wave_precedent;
        int face_treated;
        int wave_tet_mark;

        double min_length;

        int m_axis[2];


        std::vector<double> durations;

        gmds::Variable<int>* m_X;
        gmds::Variable<int>* m_Y;
        gmds::Variable<int>* m_Z;

    };
}

#endif //GMDS_DUAL_H
