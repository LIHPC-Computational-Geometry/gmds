/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "GMDSSmoothy_export.h"
#include <gmds/smoothy/AbstractSmoother.h>
#include <gmds/smoothy/SmoothingClassificationService.h>
/*----------------------------------------------------------------------------*/
#ifndef GMDS_ANGLE_BASED_QUAD_SMOOTHER_H
#	define GMDS_ANGLE_BASED_QUAD_SMOOTHER_H
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace smoothy {
/*----------------------------------------------------------------------------*/
/** \class LaplacianSmoother
 *  \brief This class provides an angle-based smoother fo smoother
 *  for 2D quad  meshes classified on geometric models.
 */
/*----------------------------------------------------------------------------*/
class GMDSSmoothy_API AngleBasedQuadSmoother : public AbstractSmoother, SmoothingClassificationService
{
 public:
	AngleBasedQuadSmoother(Mesh* AMesh, cad::GeomMeshLinker *ALinker);

	int smooth() override;
	/**@brief validation of input
	 * @return check the linker (cf AbstractSmoother class) and check
	 *         that we have only quad faces in the mesh.
	 */
	bool isValid() const override;

 private:
	/**
	 * Smooth curves only
	 */
	void smoothCurves();
	/**
	 * Smooth surfaces only
	 */
	void smoothSurfaces();
};
/*----------------------------------------------------------------------------*/
}     // namespace smoothy
/*----------------------------------------------------------------------------*/
}     // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_ANGLE_BASED_QUAD_SMOOTHER_H
/*----------------------------------------------------------------------------*/
