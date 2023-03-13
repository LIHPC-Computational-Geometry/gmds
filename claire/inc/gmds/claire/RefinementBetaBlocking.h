//
// Created by rochec on 08/08/2022.
//

#ifndef GMDS_REFINEMENTBETABLOCKING_H
#define GMDS_REFINEMENTBETABLOCKING_H

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
class LIB_GMDS_CLAIRE_API RefinementBetaBlocking
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
	RefinementBetaBlocking(Blocking2D* ABlocking2D, ParamsAero& Aparams_aero);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/

 private:
	/*-------------------------------------------------------------------*/
	/** @brief Refine a chord
	 	* \param[in] ind_chord index of the chord to refine
		*
		* \return
	 */
	void ChordRefinement(int ind_chord);
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

 private:
	/** Blocking */
	Blocking2D *m_blocking;
	/** Params pour l'aéro */
	ParamsAero m_params_aero;
	/** Ids of the chord and vector of blocks */
	std::map<int, std::vector<TCellID>> m_map_chords;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_REFINEMENTBETABLOCKING_H
