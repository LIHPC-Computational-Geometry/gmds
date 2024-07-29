/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "GMDSSmoothy_export.h"
#include <gmds/cad/GeomMeshLinker.h>
#include <gmds/utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
#ifndef GMDS_SMOOTHING_CLASSIFICATION_SERVICE_H
#	define GMDS_SMOOTHING_CLASSIFICATION_SERVICE_H
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace smoothy {
/*----------------------------------------------------------------------------*/
/** \class SmoothingClassificationService
 *  \brief This class provides attributes and methods that are shared by all
 *         smoothing algorithms that smooth classified meshes, i.e. meshes that
 *         are deeply linked to a geometric model.
 */
/*----------------------------------------------------------------------------*/
class GMDSSmoothy_API SmoothingClassificationService
{
 public:
	/**@brief constructor
	 * @param ALinker the linker that has the knowledge and connection
	 *                between geometry and mesh
	 */
	SmoothingClassificationService(cad::GeomMeshLinker *ALinker);

	/**@brief validation of the mesh model
	 * @return true if the model fits the algorithms requirement. False otherwise
	 */
	bool isValidForClassification() const;

 protected:
	/**@brief initialize all the required technical attributes$
	 */
	void init();

	/** linker*/
	cad::GeomMeshLinker *m_linker;
	/** for each curve id, we store the mesh node ids*/
	std::map<int, std::vector<TCellID>> m_c2n;
	/** for each surface id, we store the mesh node ids*/
	std::map<int, std::vector<TCellID>> m_s2n;
	/** for each volume id, we store the mesh node ids*/
	std::map<int, std::vector<TCellID>> m_v2n;
	/** local n2n connection build for the purposes of smoothing*/
	std::map<TCellID, std::vector<TCellID>> m_n2n;
};
/*----------------------------------------------------------------------------*/
}     // namespace smoothy
/*----------------------------------------------------------------------------*/
}     // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_SMOOTHING_CLASSIFICATION_SERVICE_H
/*----------------------------------------------------------------------------*/
