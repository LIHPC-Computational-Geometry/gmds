/*----------------------------------------------------------------------------*/
#ifndef GMDS_LOCALCELLTOPOLOGY_H
#define GMDS_LOCALCELLTOPOLOGY_H

#include <stddef.h>
/*----------------------------------------------------------------------------*/
namespace gmds {

    class LocalHexTopology{
    public:
        static  const size_t F2N[6][4];
        static const size_t E2N[12][2] ;

        /**@brief Give the index of the face opposed to the face given as an input*/
        static const size_t OppositeFace[6];
        /**@brief Give the index of the edges opposite to the edge given as an input */
        static const size_t OppositeEdges[12][3];
    };
/*----------------------------------------------------------------------------*/
}//namespace gmds
/*----------------------------------------------------------------------------*/
#endif //GMDS_LOCALCELLTOPOLOGY_H
/*----------------------------------------------------------------------------*/
