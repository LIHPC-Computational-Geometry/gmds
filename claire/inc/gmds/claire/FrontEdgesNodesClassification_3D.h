//
// Created by rochec on 20/12/2022.
//

#ifndef GMDS_FRONTEDGESNODESCLASSIFICATION_3D_H
#define GMDS_FRONTEDGESNODESCLASSIFICATION_3D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/claire/Front_3D.h>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API FrontEdgesNodesClassification_3D
{

 public:
	/*--------------------------------------------------------------------*/
	/** @enum  Status code for executing algorithms
	 */
	typedef enum {
		FAIL,
		SUCCESS
	} STATUS;

	/*-------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param[in] AMeshT the triangular mesh where we work on
         *  @param[in] AFront front 3D
         *  @param[in] A_EdgesClassification variable for the edges classification
         *  @param[in] A_NodesClassification variable for the nodes classification
         *
	 */
	FrontEdgesNodesClassification_3D(Mesh *AMesh, Front_3D *AFront, Variable<int>* A_EdgesClassification, Variable<int>* A_NodesClassification);

	/*-------------------------------------------------------------------*/
	/** @brief Destructor.
         *  @param
	 */
	~FrontEdgesNodesClassification_3D();
	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/
	/** @brief
	 	* \param[in]
		*
		* \return  return the mark on semi edges
	 */
	int getMarkSemiEdges();
	/*-------------------------------------------------------------------*/
	/** @brief
	 	* \param[in]
		*
		* \return  return the mark on semi nodes
	 */
	int getMarkSemiNodes();
	/*-------------------------------------------------------------------*/
	/** @brief
	 	* \param[in]
		*
		* \return  return the global feature edges
	 */
	std::vector<std::vector<TCellID>> getGlobalFeatureEdge();
	/*-------------------------------------------------------------------*/

 private:
	/*-------------------------------------------------------------------*/
	/** @brief
	 	* \param[in] e_id the id of the edge considered
		*
		* \return  return the edge classification
	 */
	int SingleEdgeClassification(TCellID e_id);
	/*-------------------------------------------------------------------*/
	/** @brief Classification of all the edges on the front
	 	* \param[in]
		*
		* \return fill the edge classification variable
	 */
	void FrontEdgesClassification();
	/*-------------------------------------------------------------------*/
	/** @brief Classification of all the nodes on the front
	 	* \param[in]
		*
		* \return fill the node classification variable
	 */
	void FrontNodesClassification();
	/*-------------------------------------------------------------------*/
	/** @brief Init the mark on semi edges, and the one on semi nodes
	 	* \param[in]
		*
		* \return
	 */
	void MarkSemiEdgesandNodes();
	/*-------------------------------------------------------------------*/
	/** @brief From a feature node and an adjacent feature edge, compute
	 	* the global feature edge
	 	* \param[in] n_id id of the starting feature node
	 	* \param[in] e_id id of the starting feature edge
		*
		* \return IDs of the edges on the global feature edge
	 */
	std::vector<TCellID> ComputeOneGlobalFeatureEdge(TCellID n_id, TCellID e_id);
	/*-------------------------------------------------------------------*/
	/** @brief Compute all the global feature edge
	 	* \param[in]
		*
		* \return
	 */
	std::vector<std::vector<TCellID>> ComputeAllGlobalFeatureEdge();
	/*-------------------------------------------------------------------*/
 private:
	/** the quad mesh we work on */
	Mesh *m_mesh;
	/** the front */
	Front_3D *m_Front;
	/** Edges classification */
	Variable<int>* m_EdgesClassification;
	/** Nodes classification */
	Variable<int>* m_NodesClassification;
	/** Mark on semi edges */
	int m_mark_semiEdges;
	/** Mark on semi nodes */
	int m_mark_semiNodes;
	/** */
	std::vector<std::vector<TCellID>> m_global_feature_edges;
};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_FRONTEDGESNODESCLASSIFICATION_3D_H
