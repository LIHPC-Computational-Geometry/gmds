/*----------------------------------------------------------------------------*/
#include <gmds/utils/LocalCellTopology.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
const size_t LocalHexTopology::F2N[6][4] = {
        {0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4},
        {1, 2, 6, 5},  {2, 3, 7, 6},  {3, 0, 4, 7}
};
const size_t LocalHexTopology::E2N[12][2] = {
        {0, 1},
        {1, 2},
        {2, 3},
        {3, 0},
        {4, 5},
        {5, 6},
        {6, 7},
        {7, 4},
        {0, 4},
        {1, 5},
        {2, 6},
        {3, 7}
};

const size_t LocalHexTopology::OppositeFace[6] = {1,0,4,5,2,3};
const size_t LocalHexTopology::OppositeEdges[12][3] = {
        {2,6,4},    //edges opp to edge 0
        {3,7,5},    //edges opp to edge 1
        {0,4,6},    //edges opp to edge 2
        {1,5,7},    //edges opp to edge 3
        {0,2,6},    //edges opp to edge 4
        {6,3,7},    //edges opp to edge 5
        {2,0,4},    //edges opp to edge 6
        {3,1,5},    //edges opp to edge 7
        {9,10,11},  //edges opp to edge 8
        {10,11,8},  //edges opp to edge 9
        {11,8,9},   //edges opp to edge 10
        {8,9,10},   //edges opp to edge 11
};
