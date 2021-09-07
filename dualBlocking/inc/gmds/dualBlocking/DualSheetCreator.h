//
// Created by calderans on 11/03/19.
//

#ifndef GMDS_DUALSHEETCREATOR_H
#define GMDS_DUALSHEETCREATOR_H

#include <gmds/ig/Mesh.h>
#include <gmds/math/Chart.h>
#include <gmds/math/Plane.h>
#include <gmds/cad/FACManager.h>
#include <gmds/cad/GeomMeshLinker.h>
#include <map>

namespace db {
    class DualSheetCreator {

    public:

        struct intersectInfo{
            gmds::math::Point p;
            gmds::math::Vector3d v;//normale du plan
            gmds::math::Vector3d d;//direction de la propagation
            gmds::TCellID e;
            int wave;
            gmds::TCellID tet;

            bool operator==(const intersectInfo& a) const
            {
                return (p == a.p && (v.X() == a.v.X() && v.Y() == a.v.Y() && v.Z() == a.v.Z()) && e == a.e && wave == a.wave);
            }

            bool operator!=(const intersectInfo& a) const
            {
                return !(p == a.p && (v.X() == a.v.X() && v.Y() == a.v.Y() && v.Z() == a.v.Z()) && e == a.e && wave == a.wave );
            }

            bool operator<(const intersectInfo& a) const
            {
                return ( e < a.e || wave < a.wave || (v.X() < a.v.X() || v.Y() < a.v.Y() || v.Z() < a.v.Z()));
            }
        };

        /*------------------------------------------------------------------------*/
        /** \brief Setting the Id of the dual sheet.
         *
         * \param AID id of the sheet
         */
        void setTetID(gmds::TCellID AID);

        /*------------------------------------------------------------------------*/
        /** \brief Setting the Id of the first tetra of the dual sheet.
         *
         * \param AID id of the tetra
         */
        void setSheetID(int AID);

        /*------------------------------------------------------------------------*/
        /** \brief Execute the algorithm of dual sheet creation starting in tetra tetID
         * and following the face on the geometric model.
         *
         * \return true if the sheet creation algorithm end or false if the sheet meet a singularity
         */
        virtual bool execute() = 0;

        /*------------------------------------------------------------------------*/
        /** \brief Get all IDs of tetra forming the dual sheet.
         *
         * \return Vector of tetra ID
         */
        std::map<gmds::TCellID,std::vector<gmds::math::Vector3d>> getSurface();

        /*------------------------------------------------------------------------*/
        /** \brief Get the sheet ID
         *
         * \return sheet ID
         */
        int getID();

        void reset();


        /*------------------------------------------------------------------------*/
        /** \brief Correct the manifoldness of the sheet and the geometry issues due to moving nodes.
         */
        void sheet_cleaning();


        virtual gmds::Mesh* buildSurfaceSheet() = 0;


        DualSheetCreator(gmds::Mesh* AMesh,
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

        virtual ~DualSheetCreator();
    protected:
        bool propagation_loop(const gmds::TCellID &ATet_ID, gmds::math::Vector3d start_norm);

        bool propagation_loop(const std::set<gmds::TCellID> &ATet_list);

        /*------------------------------------------------------------------------*/
        /** \brief Find the closest vector of a chart defined at point APoint in the tetra ATet
         */
        gmds::math::Vector3d findNormal(gmds::math::Point& APoint,
                                              const gmds::math::Vector3d& AV,
                                              gmds::Region& ATet);

        /*------------------------------------------------------------------------*/
        /** \brief Find the closest vector of a chart defined at point APoint in the edge AEdge
         */
        gmds::math::Vector3d findNormal(gmds::math::Point& APoint,
                                              const gmds::math::Vector3d& AV,
                                              gmds::Edge& AEdge);

        /*------------------------------------------------------------------------*/
        /** \brief
         */
        std::vector<intersectInfo> intersect(gmds::TCellID ATetra, gmds::math::Plane APlane,
                                             bool& AMoveNode,
                                             gmds::TCellID& AMovedNodeID);

        /*------------------------------------------------------------------------*/
        /** \brief Predicate if an edge and a plane intersect each other
         */
        int intersectEdge(gmds::math::Plane& APl, gmds::math::Point& PI,
                          gmds::Edge& AE);

        /*------------------------------------------------------------------------*/
        /** \brief Find the closest vector of \p AV in the chart \AChart
         */
        gmds::math::Vector3d closestComponentVector(const gmds::math::Vector3d& AV,
                                                          gmds::math::Chart& AChart);


        /*------------------------------------------------------------------------*/
        /** \brief compute the mean vector of the list of vectors AVectors
         */
        gmds::math::Vector3d meanVector(std::vector<gmds::math::Vector3d> AVectors);

        /*------------------------------------------------------------------------*/
        /** \brief compute the mean point of the list of points APoints
         */
        gmds::math::Point meanPoints(std::vector<gmds::math::Point> const &APoints);


        /*------------------------------------------------------------------------*/
        /** \brief compute the mean point of the list of points APoints
         */
        double averageInfo(std::vector<intersectInfo> AInfos, intersectInfo &AOutput);


        /*------------------------------------------------------------------------*/
        /** \brief Correct an intersection point on a vertex by placing it on the more orthogonal incident edge with the plane
         */
        bool onVertexCorrection(gmds::Node &ANode);

        /*------------------------------------------------------------------------*/
        /** \brief Cut the face of the tetra ATetra from the intersection points in AInfos
         */
        std::vector<intersectInfo> face_cut(gmds::TCellID ATetra, std::vector<intersectInfo> const &AInfos,
                                            bool &move_node,
                                            gmds::TCellID &moved_id);

        /*------------------------------------------------------------------------*/
        /** \brief Compute the Runge-Kutta-4 in the face
         */
        intersectInfo RK_4(const gmds::Face &ATriangle, intersectInfo Ainfo, bool &move_node, gmds::TCellID &moved_id);


        bool isAxisSingularity(gmds::TCellID AId, const gmds::math::Vector3d& AV);

        int isSingular(gmds::TCellID AId, gmds::math::Vector3d& AVOut);

        void clear();


    protected:

        //Following members are used for each dual sheet
        //-----------------------------------------------
        //The id of the sheet
        int sheet_ID;
        //The id of the starting tet
        gmds::TCellID tetID;
        //A map used to list all tet id of the surface and their normal vector of the surface in the tet
        std::map<gmds::TCellID,std::vector<gmds::math::Vector3d>> surface;
        //-----------------------------------------------


        //Following members are const for all dual sheets
        //-----------------------------------------------
        //The tet mesh
        gmds::Mesh* m_mesh;
        //The tet used to build the plan surface
        gmds::Mesh m_smesh;
        //Linker and Manager to get the geometrical informations from the mesh
        gmds::cad::GeomMeshLinker &linker;
        gmds::cad::FACManager &manager;
        //The min length of all the edge in the mesh, used as a threshold in the creation algorithm
        double min_length;
        //Mesh variable that store charts at every vertex of the mesh: used in the face cutting algorithm
        gmds::Variable<gmds::math::Chart>* m_vertex_chart; //only read
        //Mesh variable that store if the tet is singular
        gmds::Variable<int>* m_singularities; //only read
        //-----------------------------------------------

        //Following members are reseted at each dual sheet creation
        //-----------------------------------------------
        int propagation_round;

        //Variable storing the intersection info on the edges
        gmds::Variable<std::vector<intersectInfo>>* var_info;

        //Variable on region used to identify the propagation wave
        gmds::Variable<int>* m_propagation_round;

        //mark edges in the previous wave during sheet creation
        int wave_precedent;
        //mark face who have been treated during the surface creation
        int face_treated;
        //mark the tetra who have been treated during the surface creation
        int mark_wave_tet;
        //-----------------------------------------------

        gmds::Variable<int>* wave_id;
        gmds::Variable<int>* edge_Gdim;
        gmds::Variable<int>* node_Gdim;

        gmds::Variable<int>* m_part;



        int test_mark;

    };
}


#endif //GMDS_DUALSHEETCREATOR_H
