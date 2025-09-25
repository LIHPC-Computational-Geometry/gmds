/*----------------------------------------------------------------------------*/
#ifndef GMDS_MORPHMESH_FASTLOCALIZE_H
#define GMDS_MORPHMESH_FASTLOCALIZE_H
/*----------------------------------------------------------------------------*/
#include "GMDSmorphMesh_export.h"
#	include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
#include <map>
#include <ANN/ANN.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
class GMDSmorphMesh_API FastLocalize {
 public:

	/*-------------------------------------------------------------------*/
	/** @brief Constructor.
		*  @param AMesh the mesh where we work on
	 */
	FastLocalize(const std::vector<Node> &ANodes);
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
	TCellID find(const math::Point& APoint);
	/*-------------------------------------------------------------------*/

 private:
	void buildANNTree();
 private:
	/** mesh we work on */
	std::vector<Node> m_nodes;
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
#endif     // GMDS_MORPHMESH_FASTLOCALIZE_H
/*----------------------------------------------------------------------------*/
