/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "GMDSSmoothy_export.h"
#include <gmds/cad/GeomMeshLinker.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/smoothy/AbstractSmoother.h>
#include <gmds/smoothy/SmoothingClassificationService.h>
#include <gmds/utils/CommonTypes.h>
#include <gmds/utils/Exception.h>
/*----------------------------------------------------------------------------*/
#ifndef GMDS_LAPLACIANSMOOTHER_3C_H
#	define GMDS_LAPLACIANSMOOTHER_3C_H
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace smoothy {
/*----------------------------------------------------------------------------*/
/** \class LaplacianSmoother
 *  \brief This class provides a straightforward implementation of a Laplacian
 *         smoother for meshes classifiend on geometric models. Curves, surfaces
 *         and volumes are smoothed individually.
 */
/*----------------------------------------------------------------------------*/
class GMDSSmoothy_API LaplacianSmoother3C : public AbstractSmoother, SmoothingClassificationService
{
 public:
	LaplacianSmoother3C(Mesh *AMesh, cad::GeomMeshLinker *ALinker);
	bool isValid() const override;
	int smooth() override;
	void smoothCurves();
	void smoothSurfaces();
	void smoothVolumes();
};
/*----------------------------------------------------------------------------*/
}     // namespace smoothy
/*----------------------------------------------------------------------------*/
}     // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_LAPLACIANSMOOTHER_3C_H
/*----------------------------------------------------------------------------*/
