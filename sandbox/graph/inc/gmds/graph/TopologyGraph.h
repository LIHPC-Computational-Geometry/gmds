#ifndef GMDS_TOPOLOGYGRAPH_H
#define GMDS_TOPOLOGYGRAPH_H

#include <vector>
#include <gmds/cadfac/FACManager.h>
#include <gmds/math/Chart.h>
#include <gmds/graph/MinCut.h>


namespace gmds {

    namespace graph {
        class TopologyGraph {

        public:
            TopologyGraph(gmds::cad::FACManager* Amanager, gmds::cad::GeomMeshLinker* Alinker, gmds::Mesh* Amesh);

            int getGeomSurfID(gmds::TCellID AID);

            std::vector<gmds::TCellID> getGeomSurf(int AsurfID);

            int pickSurf(gmds::TCellID AID,std::vector<TCellID> &AIDsOutPut);

            /*
             * Run the cutting algorithm, need to pick surfaces before
             */
            int generateCut();

            /*
             * Method to undo action for the hmi
             */
            int undoLastCut();

            int resetCut();

            int undoSelection();

            /*
             * Method that set the value of cut to 1 for all tets, in other words when there is no cut
             */
            void noCut();

            void cutFace(TCellID ANode1, TCellID ANode2);


        protected:
            int geomEdgeReclassification(std::vector<TCellID> AEdges);

            void computeWeight(std::vector<gmds::TCellID> AFaces);

            std::vector<int> findSelectionBoundary();

            void getCut();

            void surfaceReconstruction(std::vector<gmds::TCellID>);

            void surfaceGraph();

            void cutOptim();

            void boundaryMarking();

            void cutDirectionWeight(int ASurfID);

        protected:
            std::vector<std::vector<int>> m_graph;
            gmds::Variable<int>* BND_COLOR;

            gmds::cad::FACManager* m_manager;
            gmds::cad::GeomMeshLinker* m_linker;

            gmds::Variable<int>* m_surface_color;
            gmds::Variable<int>* m_isomorph;
            gmds::Variable<int>* m_isomorph_tet;

            gmds::Variable<gmds::math::Chart>* m_vertex_chart;
            gmds::Variable<double>* m_face_weight;

            gmds::Variable<int>* m_cut_t;
            gmds::Variable<int>* m_cut_f;

            gmds::Mesh* m_mesh;

            std::vector<int> m_selection;
            std::map<int,std::vector<gmds::TCellID>> m_boundaryEdges;
            int mark_face_surf;
            int mark_tet_offset;


            MinCut* m_minCut;

            //Use to store the number of cut done in the model
            //Can be use to know the last cut done
            int cut_id;

            void computeEdgeFrameWeight(Edge e, math::Vector3d AtargetVec);
        };

    }
}
#endif //GMDS_TOPOLOGYGRAPH_H
