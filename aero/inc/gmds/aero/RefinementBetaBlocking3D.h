//
// Created by rochec on 08/12/23.
//

#ifndef GMDS_REFINEMENTBETABLOCKING3D_H
#define GMDS_REFINEMENTBETABLOCKING3D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_AERO_export.h"
#include "gmds/ig/Mesh.h"
#include <gmds/aero/Blocking3D.h>
#include <gmds/aero/Params.h>
#include <string>
#include <map>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_AERO_API RefinementBetaBlocking3D
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
         *  @param[in] ABlocking3D the blocking 3D
         *  @param[in] Aparams_aero parameters for aero algorithm
         *
	 */
	RefinementBetaBlocking3D(Blocking3D* ABlocking3D,
	                         Variable<int>* Avar_LayerID,
	                         double Asize_first_edge);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/
 private:
	/*-------------------------------------------------------------------*/
	/** @brief Refine a sheet
	 	*
	 	* \param[in] sheet vector of the edge ids of the sheet to refine
	 	*
	 	* \return
	 */
	void sheetRefinement(std::vector<std::pair<TCellID,TCellID>>* sheet);
	/*-------------------------------------------------------------------*/
	/** @brief
	 	*
	 	* \param[in]
	 	*
	 	* \return
	 */
	void computeSheetsToRefine();
	/*-------------------------------------------------------------------*/
	/** @brief
	 	*
	 	* \param[in] b_id Block ID
		* \param[in] bf_id Block Face ID
	 	*
	 	* \return
	 */
	std::vector<std::pair<TCellID,TCellID>> computeOneSheet(TCellID b_id, TCellID bf_id);
	/*-------------------------------------------------------------------*/
	/** @brief
	 	*
	 	* \param[in] b_id Block ID
	 	* \param[in] bf_id Block Face ID
	 	*
	 	* \return
	 */
	std::vector<std::pair<TCellID,TCellID>> compute4adjacentFacesandEdgestoFaceinBlock(TCellID b_id, TCellID bf_id);
	/*-------------------------------------------------------------------*/

 private:
	/** Blocking */
	Blocking3D *m_blocking3D;
	/** Variable layer ID */
	Variable<int>* m_var_LayerID;
	/** Params pour l'aéro */
	double m_size_first_edge;
	/** Ids of the sheets and vector of pair (block,face) to refine the sheet */
	std::vector<std::vector<std::pair<TCellID,TCellID>>> m_sheets;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_REFINEMENTBETABLOCKING3D_H
