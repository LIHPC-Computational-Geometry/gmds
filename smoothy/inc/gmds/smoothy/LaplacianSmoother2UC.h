/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "GMDSSmoothy_export.h"
#include <gmds/cad/GeomMeshLinker.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/smoothy/AbstractSmoother.h>
#include <gmds/utils/CommonTypes.h>
#include <gmds/utils/Exception.h>
/*----------------------------------------------------------------------------*/
#ifndef GMDS_LAPLACIANSMOOTHER_2UC_H
#	define GMDS_LAPLACIANSMOOTHER_2UC_H
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace smoothy {
/*----------------------------------------------------------------------------*/
/** \class LaplacianSmoother
 *  \brief This class provides a straightforward implementation of a Laplacian
 *         smoother for unclassified 2D meshes.
 */
/*----------------------------------------------------------------------------*/
class GMDSSmoothy_API LaplacianSmoother2UC : public AbstractSmoother
{
 public:
	/**@brief constructor
	 * @param AMesh the mesh we work on
	 */
	LaplacianSmoother2UC(Mesh *AMesh);
	/** @brief check if the mesh is valid for performing the algorithm.
	 * More specifically the mesh must not have regions and the relation N2F
	 * must be available.
	 *
	 * @return true if the mesh is valid, false otherwise.
	 */
	bool isValid() const override;
	/** @brief perform the smoothing algorithm
	 * @return 1 if the mesh is smoothed
	 */
	int smooth() override;
};
/*----------------------------------------------------------------------------*/
}     // namespace smoothy
/*----------------------------------------------------------------------------*/
}     // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_LAPLACIANSMOOTHER_2UC_H
/*----------------------------------------------------------------------------*/
