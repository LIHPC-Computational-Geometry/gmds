/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEODHEXMESHER_H
#define GMDS_GEODHEXMESHER_H
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include "LIB_GMDS_GEOD_HONEY_COMB_export.h"

/*----------------------------------------------------------------------------*/
namespace gmds{
    /*------------------------------------------------------------------------*/
    /** @class This class gathers some basic algorithm to build a honeycomb
     *          mesh of a spherical domain
     */
    class LIB_GMDS_GEOD_HONEY_COMB_API GeodHexMesher{
    public:
        typedef enum {
            GEOD_SUCCESS,
            GEOD_FAILURE
        } OpResult;
        /*-------------------------------------------------------------------*/
        /** @brief Constructor
         */
        GeodHexMesher(const double radius=1.0);
        /*-------------------------------------------------------------------*/
        /** @brief Destructor
         */
        virtual ~GeodHexMesher();
        /*-------------------------------------------------------------------*/
        /** @brief Execute the algorithm to get the final mesh
         * @return The result of the operation (success or failure)
         */
        OpResult execute();

    private:
        double m_radius;

    };
    /*------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_GEODHEXMESHER_H
/*----------------------------------------------------------------------------*/
