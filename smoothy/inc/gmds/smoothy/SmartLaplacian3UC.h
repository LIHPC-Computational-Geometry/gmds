/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "AbstractSmoother.h"
#include "GMDSSmoothy_export.h"
#include <gmds/utils/CommonTypes.h>
#include <gmds/utils/Exception.h>
/*----------------------------------------------------------------------------*/
#ifndef GMDS_SMARTLAPLACIAN_3UC_H
#	define GMDS_SMARTLAPLACIAN_3UC_H
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
class GMDSSmoothy_API SmartLaplacian3UC : public AbstractSmoother
{
 public:
	/**@brief constructor
	 * @param AMesh the mesh we work on
	 */
	SmartLaplacian3UC(Mesh *AMesh);
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

 private:
	void getAdjacentNodes(const Region& ARegion, const TCellID &ANode1, TCellID &ANode2, TCellID &ANode3, TCellID &ANode4);

	void getAdjacentQuality(const Node& ANode, std::vector<double>& AQualities);
};
/*----------------------------------------------------------------------------*/
}     // namespace smoothy
/*----------------------------------------------------------------------------*/
}     // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_SMARTLAPLACIAN_3UC_H
/*----------------------------------------------------------------------------*/
