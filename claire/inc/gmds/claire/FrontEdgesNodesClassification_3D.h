//
// Created by rochec on 20/12/2022.
//

#ifndef GMDS_FRONTEDGESNODESCLASSIFICATION_3D_H
#define GMDS_FRONTEDGESNODESCLASSIFICATION_3D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/claire/Front_3D.h>
#include <gmds/claire/FastLocalize.h>
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
	/*--------------------------------------------------------------------*/
	struct Global_Feature_Edge
	{
		TCellID     					Start_n_id;
		TCellID 							End_n_id;
		std::vector<TCellID> 		edges_id;
	};
	/*-------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param[in] AMeshT the triangular mesh where we work on
         *  @param[in] AFront front 3D
         *  @param[in] A_EdgesClassification variable for the edges classification
         *  @param[in] A_NbrFeatureEdgesAroundNode variable for the nodes classification
         *
	 */
	FrontEdgesNodesClassification_3D(Mesh *AMesh, Front_3D *AFront, Variable<int>* A_EdgesClassification, Variable<int>*A_NbrFeatureEdgesAroundNode,
	                                 Mesh *AMesh_T, FastLocalize *Afl, Variable<math::Vector3d>* A_VectorField);

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
		* \return  return the mark
	 */
	TInt getMarkEdgesTemplates();
	/*-------------------------------------------------------------------*/
	/** @brief
	 	* \param[in]
		*
		* \return  return the mark
	 */
	TInt getMarkNodesTemplates();
	/*-------------------------------------------------------------------*/
	/** @brief
	 	* \param[in]
		*
		* \return  return the variable on nodes of the nodes classification
	 */
	Variable<int>* getVarNodesClassification();
	/*-------------------------------------------------------------------*/
	/** @brief
	 	* \param[in]
		*
		* \return  return the variable on nodes of the nodes classification
	 */
	Variable<int>* getVarNbrFeatureEdgesAroundNode();
	/*-------------------------------------------------------------------*/
	/** @brief Return the classification of the front node n
	 	* \param[in] n_id the id of the node considered
		*
		* \return  return the node classification (0: regular, 1: 3 CORNER, ...)
	 */
	int singleNodeClassification(TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 	* \param[in] n_id node id
		*
		* \return  return the mark
	 */
	int getNbrSingularEdgesAroundNode(TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 	* \param[in] e_id edge id
		*
		* \return  return the mark
	 */
	int getEdgeType(TCellID e_id);
	/*-------------------------------------------------------------------*/
 private:
	/*-------------------------------------------------------------------*/
	/** @brief Return the classification of a single edge e
	 	* \param[in] e_id the id of the edge considered
		*
		* \return  return the edge classification (0: side, 1: corner, 2:end, 3:reversal)
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
	/** @brief Return true if the configuration around the node n is valid
	 * to apply one of the implemented templates on nodes.
	 	* \param[in] n_id id of the node
		*
		* \return true if the node is valid for a template, false otherwise.
	 */
	bool isValidNodeForTemplate(TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief From a feature node and an adjacent feature edge, compute
	 	* the global feature edge
	 	* \param[in] n_id id of the starting feature node
	 	* \param[in] e_id id of the starting feature edge
		*
		* \return IDs of the edges on the global feature edge
	 */
	Global_Feature_Edge ComputeOneGFE(TCellID n_id, TCellID e_id);
	/*-------------------------------------------------------------------*/
	/** @brief Compute all the global feature edge
	 	* \param[in]
		*
		* \return
	 */
	std::vector<Global_Feature_Edge> ComputeAllGFE();
	/*-------------------------------------------------------------------*/
	/** @brief Return true if the path is valid to apply templates.
	 	* \param[in] GFE the "path"
		*
		* \return true if the path is valid for templates, false otherwise.
	 */
	bool isThisPathValidForTemplates(Global_Feature_Edge& GFE);
	/*-------------------------------------------------------------------*/
	/** @brief Return true if the path is a valid loop path to apply templates.
	 	* \param[in] GFE the "path"
		*
		* \return true if the path is valid for templates, false otherwise.
	 */
	bool isThisValidLoopPath(Global_Feature_Edge& GFE);
	/*-------------------------------------------------------------------*/
	/** @brief
	 	* \param[in]
		*
		* \return
	 */
	void ComputeValid_GFE();
	/*-------------------------------------------------------------------*/
	/** @brief
	 	* \param[in]
		*
		* \return
	 */
	void ComputeValidLoop_GFE();
	/*-------------------------------------------------------------------*/
	/** @brief Return true if the node is valid for being internal to a
	 * path of edges of type type_path.
	 * \param[in] n_id id of the node
	 * \param[in] type_path type of the path (1: CORNER, 2: END, 3: REVERSAL)
	 *
	 * \return
	 */
	bool isValidNodeForInsidePath(TCellID n_id, int type_path);
	/*-------------------------------------------------------------------*/
	/** @brief Return true if the node is valid for being a limit to a
	 * path of edges of type type_path.
	 * \param[in] n_id id of the node
	 * \param[in] type_path type of the path (1: CORNER, 2: END, 3: REVERSAL)
	 *
	 * \return
	 */
	bool isValidNodeForPathLimit(TCellID n_id, int type_path);
	/*-------------------------------------------------------------------*/
	/** @brief Clean the edges classification of the front
	 	* \param[in]
		*
		* \return
	 */
	void FrontEdgesClassificationCleaner();
	/*-------------------------------------------------------------------*/
	/** @brief
	 	* \param[in]
		*
		* \return
	 */
	void SemiEdgesCleaner();
	/*-------------------------------------------------------------------*/
 private:
	/** the quad mesh we work on */
	Mesh *m_mesh;
	/** the tet mesh we work on */
	Mesh *m_mesh_T;
	/** k-d tree */
	FastLocalize *m_fl;
	/** Vector Field for extrusion */
	Variable<math::Vector3d>* m_VectorField;
	/** the front */
	Front_3D *m_Front;
	/** Edges classification */
	Variable<int>* m_EdgesClassification;
	/** Edges classification */
	Variable<int>* m_NodesClassification;
	/** Nodes classification */
	Variable<int>*m_NbrFeatureEdgesAroundNode;
	/** Mark edges for templates */
	TInt m_mark_EdgesForTemplates;
	/** Mark nodes for templates */
	TInt m_mark_NodesForTemplates;
	/** */
	std::vector<std::vector<TCellID>> m_global_feature_edges;
	/** */
	std::vector<Global_Feature_Edge> m_All_global_feature_edges;
};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_FRONTEDGESNODESCLASSIFICATION_3D_H
