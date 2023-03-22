#ifndef GMDS_ELGMORPHING_ELGMORPHING_H
#define GMDS_ELGMORPHING_ELGMORPHING_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_ELGMORPHING_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/math/Point.h>
#include <string>
//#include <map>
//#include <fstream>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace elgmorphing {
/*----------------------------------------------------------------------------*/
/** \class ElgMorphing
 *  \brief  Class
 */
class LIB_GMDS_ELGMORPHING_API ElgMorphing{

 public:
	/*-------------------------------------------------------------------------*/
	/** @enum  Status code for executing algorithms
	 */
	typedef enum {
		FAIL,
		SUCCESS
	} STATUS;
	/*-------------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param
	 */
	explicit ElgMorphing();
	/*-------------------------------------------------------------------------*/
	/** @brief Default destructor.
         *  @param
	 */
	virtual ~ElgMorphing() =default;
	/*-------------------------------------------------------------------------*/
	/** \brief
	 */
	ElgMorphing::STATUS execute();
	/*-------------------------------------------------------------------------*/
	/** \brief Boundingbox computation
	 */
	void computeBoundingBox(gmds::Mesh* AMesh,
	                        std::string AGroupName,
									gmds::TCoord minXYZ[3],
	                        gmds::TCoord maxXYZ[3]) const;
	/*-------------------------------------------------------------------------*/
	/** \brief Boundingbox transform.
	 */
	gmds::math::Point bbtransform2d(gmds::math::Point APt,
	                              gmds::TCoord minXYZ_orig[3], gmds::TCoord maxXYZ_orig[3],
	                              gmds::TCoord minXYZ_dest[3], gmds::TCoord maxXYZ_dest[3]) const;
	/*-------------------------------------------------------------------------*/

 protected:
	/** pointer to mesh we work on */
	Mesh* m_mesh;
};
/*----------------------------------------------------------------------------*/
}  // namespace elgmorphing
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/
#endif  // #define GMDS_ELGMORPHING_ELGMORPHING_H

