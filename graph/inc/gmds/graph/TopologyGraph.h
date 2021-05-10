#ifndef GMDS_TOPOLOGYGRAPH_H
#define GMDS_TOPOLOGYGRAPH_H

#include <vector>
#include <gmds/cad/FACManager.h>



namespace gmds {

    namespace graph {
        class TopologyGraph {

        public:
            TopologyGraph(gmds::cad::FACManager* Amanager, gmds::Mesh* Amesh);

            bool buildGraph();

            void printGraph();

            void writeGML();

            void colorSurfaces();

            bool findIso();

        protected:
            std::vector<std::vector<int>> m_graph;

            gmds::cad::FACManager* m_manager;

            gmds::Variable<int>* m_surface_color;
            gmds::Variable<int>* m_isomorph;

            gmds::Mesh* m_mesh;

        };

    }
}
#endif //GMDS_TOPOLOGYGRAPH_H
