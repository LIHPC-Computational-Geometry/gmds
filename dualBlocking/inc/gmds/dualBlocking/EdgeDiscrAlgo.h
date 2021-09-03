//
// Created by simon on 12/07/2021.
//

#ifndef GMDS_EDGEDISCRALGO_H
#define GMDS_EDGEDISCRALGO_H

#include <glpk.h>
#include <gmds/ig/Mesh.h>
#include <map>

    namespace db {
        class EdgeDiscrAlgo {

        public:
            EdgeDiscrAlgo(gmds::Mesh* AMesh);

            void init();

            void intervalAssignment();

            void boundaryDiscretization();

            void edgesGrouping();
            void edgesGroupingBlock();

            std::map<gmds::TCellID,int> getEdgeDiscretization();


        protected:

            gmds::Mesh* m_blocks;

            //A map that stores the target nb of edge subdivision for each block edge
            std::map<gmds::TCellID,int> targetDiscretization;
            //A map that stores the computed nb of edge subdivision for each block edge
            std::map<gmds::TCellID,int> computedDiscretization;
            //A map that stores the weight for each block edge
            std::map<int,double> edgeWeight;

            //Ugly way to store the 3 edges group for each block
            std::map<int,std::tuple<std::tuple<gmds::TCellID,gmds::TCellID,gmds::TCellID,gmds::TCellID>,std::tuple<gmds::TCellID,gmds::TCellID,gmds::TCellID,gmds::TCellID>,std::tuple<gmds::TCellID,gmds::TCellID,gmds::TCellID,gmds::TCellID>>> edgeGroupBlock;
            //Face method
            std::map<int,std::tuple<std::tuple<gmds::TCellID,gmds::TCellID>,std::tuple<gmds::TCellID,gmds::TCellID>>> edgeGroupFace;

        };
    }


#endif //GMDS_EDGEDISCRALGO_H
