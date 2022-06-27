//
// Created by rochec on 27/06/2022.
//

#ifndef GMDS_INTERVALASSIGNMENT_2D_H
#define GMDS_INTERVALASSIGNMENT_2D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include "gmds/ig/Mesh.h"
#include <gmds/ig/Blocking2D.h>
#include <gmds/claire/Params.h>
#include <string>
#include <map>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API IntervalAssignment_2D
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
         *  @param[in] ABlocking2D the blocking 2D
         *  @param[in] Aparams_aero parameters for aero algorithm
         *
	 */
	IntervalAssignment_2D(Blocking2D* ABlocking2D, ParamsAero Aparams_aero);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/

 private:
	/*-------------------------------------------------------------------*/
	/** @brief Compute the map of chords of the blocking.
	 	* \param[in] ABlocking the blocking
		*
		* \return  a map with (int, std::vector<TCellID>) where the int is the id
	 	* of the chord and the vector is filled with the ids of the chord edges
	 */
	std::map<int, std::vector<TCellID>> ComputeChords();
	/*-------------------------------------------------------------------*/
	/** @brief Get the opposite edges
	 	* \param[in] e_id the edge we want the opposites
		*
		* \return  std::vector<TCellID> vector of one or two opposite edges
	 	* of e_id
	 */
	std::vector<TCellID> ComputeOppositeEdges(TCellID e_id);
	/*-------------------------------------------------------------------*/
	/** @brief Get the opposite edges
	 	* \param[in] e_id the edge we want to know the constraint
	 	* \param[out] N_ideal the ideal number of cells on the edge
	 	* \param[out] hardConstraint is the constraint on the edge hard or not
		*
		* \return
	 */
	void EdgeConstraint(TCellID e_id, int N_ideal, bool hardConstraint);
	/*-------------------------------------------------------------------*/

 private:
	/** Blocking */
	Blocking2D *m_blocking;
	/** Params pour l'aéro */
	ParamsAero m_params_aero;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_INTERVALASSIGNMENT_2D_H
