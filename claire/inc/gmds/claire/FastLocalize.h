/*----------------------------------------------------------------------------*/
#ifndef GMDS_CLAIRE_FASTLOCALIZE_H
#define GMDS_CLAIRE_FASTLOCALIZE_H
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#	include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
#include <map>
#include <ANN/ANN.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API FastLocalize {
 public:

	/*-------------------------------------------------------------------*/
	/** @brief Constructor.
		*  @param AMesh the mesh where we work on
	 */
	FastLocalize(Mesh *AMesh);
	/*-------------------------------------------------------------------*/
	/** @brief Destructor.
	 */
	virtual ~FastLocalize();

	/*-------------------------------------------------------------------*/
	/** @brief Given an input point \p APoint return the cell data that
	 *  indicate what is the closest cell of \APoint (it can be a node or
	 *  a face)
	*  @param APoint point to localize
	*  @return the cell data indicating where \p APoint is in this->m_mesh
	 */
	Cell::Data find(const math::Point& APoint);
	/*-------------------------------------------------------------------*/
	/** @brief Given an input point \p APoint return the TCellID that
	 *  indicate what is the closest region of \APoint
	*  @param APoint point to localize
	*  @return the TCellID indicating where \p APoint is in this->m_mesh
	 */
	TCellID findTetra(const math::Point& APoint);
	/*-------------------------------------------------------------------*/

 private:
	void buildANNTree();
 private:
	/** mesh we work on */
	Mesh *m_mesh;
	std::map<int,TCellID > m_ann2gmds_id;
	ANNkd_tree*	 m_kdTree;
	ANNidxArray	 m_nnIdx;					// near neighbor indices
	ANNdistArray m_dists;					// near neighbor distances
	ANNpoint 	 m_queryPt;  // query point
	ANNpointArray m_dataPts;  // data points
};
/*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_CLAIRE_FASTLOCALIZE_H
/*----------------------------------------------------------------------------*/
