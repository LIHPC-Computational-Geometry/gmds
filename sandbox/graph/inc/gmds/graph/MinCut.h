//
// Created by calderans on 23/07/20.
//

#ifndef GMDS_MINCUT_H
#define GMDS_MINCUT_H

#include <map>
#include <gmds/utils/CommonTypes.h>
#include <gmds/ig/Mesh.h>
#include <gmds/math/Chart.h>
#include <queue>

namespace gmds {

    namespace graph {
        class MinCut {
        public:
            MinCut(Mesh* AMesh);

            void graph_cut(Variable<int>* ATetAssign,
                           Variable<int>* ATetResult,
                           Variable<double> *AFaceWeight);

            void graph_cut(std::map<TCellID,int> ATetList,
                           std::vector<TCellID> AFaceList,
                           Variable<int>* ATetAssign,
                           Variable<int>* ATetResult,
                           Variable<double> *AFaceWeight);

            std::vector<TCellID> shortest_path(TCellID ASource,
                                               std::vector<TCellID> ANodeGraph,
                                               std::map<TCellID,double> AEdgeWeight,
                                               std::vector<TCellID> ATargetNodes);

        protected:
            Mesh* m_mesh;
        };
    }
}

#endif //GMDS_MINCUT_H
