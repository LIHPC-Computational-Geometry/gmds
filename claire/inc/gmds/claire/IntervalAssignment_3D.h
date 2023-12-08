//
// Created by rochec on 23/03/2023.
//

#ifndef GMDS_INTERVALASSIGNMENT_3D_H
#define GMDS_INTERVALASSIGNMENT_3D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include "gmds/ig/Mesh.h"
#include <gmds/claire/Blocking3D.h>
#include <gmds/claire/Params.h>
#include <string>
#include <map>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API IntervalAssignment_3D
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
         *  @param[in] ACtrlPnts3D the control points of the curved blocking 3D @p ABlocking3D
         *  @param[in] Aparams_aero parameters for aero algorithm
         *  @param[in] AedgesDiscretization variable to store the discretization of each edge of the blocking @p ABlocking3D
         *
	 */
	IntervalAssignment_3D(Blocking3D* ABlocking3D,
	                      Blocking3D* ACtrlPnts3D,
	                      ParamsAero& Aparams_aero,
	                      Variable<int>* AedgesDiscretization);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/

 private:
	/*-------------------------------------------------------------------*/
	/** @brief Compute the map of chords of the 3D blocking.
	 	* \param[in] e_id starting edge id
		*
		* \return a map with (int, std::vector<TCellID>) where the int is the id
	 	* of the chord and the vector is filled with the ids of the chord edges
	 */
	std::vector<TCellID> ComputeSingleSheet(TCellID e_id);
	/*-------------------------------------------------------------------*/
	/** @brief Compute the map of sheets of the 3D blocking.
	 	* \param[in]
		*
		* \return  a map with (int, std::vector<TCellID>) where the int is the id
	 	* of the sheet and the vector is filled with the ids of the sheet edges
	 */
	std::map<int, std::vector<TCellID>> ComputeSheets();
	/*-------------------------------------------------------------------*/
	/** @brief Get the constraint of the edge e_id
	 	* \param[in] e_id the edge we want to know the constraint
	 	* \param[out] N_ideal the ideal number of cells on the edge
	 	* \param[out] hardConstraint is the constraint on the edge hard or not
		*
		* \return
	 */
	void EdgeConstraint(TCellID e_id, int &N_ideal, bool &hardConstraint);
	/*-------------------------------------------------------------------*/
	/** @brief Get the number of cells for the chord
	 	* \param[in] chord the chord
		*
		* \return the number of cells for all the edges of the chord
	 */
	int ComputeSheetDiscretization(const std::vector<TCellID>& sheet);
	/*-------------------------------------------------------------------*/
	/** @brief Get the length of the curved Bezier block edge
	 	*
	 	* \param[in] Ae the block edge
		*
		* \return
	 */
	double BezierEdgeLength(const Edge Ae);
	/*-------------------------------------------------------------------*/

 private:
	/** 3D Blocking */
	Blocking3D *m_Blocking3D;
	/** Control Points of the 3D curved Blocking m_Blocking3D */
	Blocking3D *m_CtrlPts3D;
	/** Params pour l'aéro */
	ParamsAero m_params_aero;
	/** Variable on the edges, filled with the edge discretization */
	Variable<int>* m_edgesDiscretization;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_INTERVALASSIGNMENT_3D_H
