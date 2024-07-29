/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/utils/Exception.h>
#include <gmds/utils/CommonTypes.h>
#include <gmds/cad/GeomMeshLinker.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/smoothy/AbstractSmoother.h>
#include <gmds/smoothy/SmoothingClassificationService.h>
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
        class GMDSSmoothy_API LaplacianSmoother: public AbstractSmoother, SmoothingClassificationService{
        public:
            LaplacianSmoother(Mesh* AMesh, cad::GeomMeshLinker* ALinker);
	         bool isValid() const override;
	         int smooth() override;
            void smoothCurves();
            void smoothSurfaces();
            void smoothVolumes();

        };
/*----------------------------------------------------------------------------*/
    } // end namespace cad
/*----------------------------------------------------------------------------*/
} // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif //GMDS_LAPLACIANSMOOTHER_H
/*----------------------------------------------------------------------------*/
