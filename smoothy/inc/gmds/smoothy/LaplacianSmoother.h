/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/utils/Exception.h>
#include <gmds/utils/CommonTypes.h>
#include <gmds/cad/GeomMeshLinker.h>
#include <gmds/cad/FACManager.h>
#include <gmds/smoothy/AbstractSmoother.h>
#include "GMDSSmoothy_export.h"
/*----------------------------------------------------------------------------*/
#ifndef GMDS_LAPLACIANSMOOTHER_H
#define GMDS_LAPLACIANSMOOTHER_H
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
    namespace smoothy{
/*----------------------------------------------------------------------------*/
/** \class LaplacianSmoother
 *  \brief This class provides a straightforward implementation of a Laplacian
 *         smoother for meshes classifiend on geometric models. Curves, surfaces
 *         and volumes are smoothed individually.
 */
/*----------------------------------------------------------------------------*/
        class GMDSSmoothy_API LaplacianSmoother: public AbstractSmoother{
        public:
            LaplacianSmoother(cad::GeomMeshLinker* ALinker);
            void smoothCurves(const int ANbIterations=1);
            void smoothSurfaces(const int ANbIterations=1);
            void smoothVolumes(const int ANbIterations=1);

        };
/*----------------------------------------------------------------------------*/
    } // end namespace cad
/*----------------------------------------------------------------------------*/
} // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif //GMDS_LAPLACIANSMOOTHER_H
/*----------------------------------------------------------------------------*/
