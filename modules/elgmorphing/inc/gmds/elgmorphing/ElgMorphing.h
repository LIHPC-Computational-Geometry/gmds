#ifndef GMDS_ELGMORPHING_ELGMORPHING_H
#define GMDS_ELGMORPHING_ELGMORPHING_H

/*----------------------------------------------------------------------------*/
#include "GMDSelgmorphing_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/math/Point.h>
#include <set>
#include <string>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace elgmorphing {
/*----------------------------------------------------------------------------*/
/** \class ElgMorphing
 *  \brief  Class
 */
class GMDSelgmorphing_API ElgMorphing{

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
	/** \brief
	 */
	void setMeshOrig(gmds::Mesh* AMesh) {mesh_orig_ = AMesh;}
	/*-------------------------------------------------------------------------*/
	/** \brief
	 */
	void setMeshDest(gmds::Mesh* AMesh) {mesh_orig_ = AMesh;}

	/**@brief Identify the nodes that are part of the cells groups
	 * which names begin with @Astring
	    	* @param AMesh
         * @param AString the prefix we look for in the groups names
         *
         * @return a container of those nodes
	 */
	std::set<TCellID> compoundGrpNodes(gmds::Mesh* AMesh,
												  std::string AString);

	/*-------------------------------------------------------------------------*/
	/**@brief Identify the nodes that are part of the cells groups
	 * which names begin with @Astring, but only those on the boundary
	    	* @param AMesh
         * @param AString the prefix we look for in the groups names
         *
         * @return a container of those nodes
	 */
	std::set<TCellID> compoundGrpBnd(gmds::Mesh* AMesh,
												std::string AString);

	/*-------------------------------------------------------------------------*/
	/** \brief Boundingbox computation
	 */
	void computeBoundingBox(gmds::Mesh* AMesh,
	                        std::string AGroupName,
									gmds::TCoord minXYZ[3],
	                        gmds::TCoord maxXYZ[3]) const;
	void computeBoundingBox(const gmds::Mesh* AMesh,
	                        const std::set<TCellID>& AIDs,
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
	Mesh* mesh_orig_;
	Mesh* mesh_dest_;

	TInt mark_fixed_nodes_;
};
/*----------------------------------------------------------------------------*/
}  // namespace elgmorphing
/*----------------------------------------------------------------------------*/
}  // namespace gmds
/*----------------------------------------------------------------------------*/
#endif  // #define GMDS_ELGMORPHING_ELGMORPHING_H

