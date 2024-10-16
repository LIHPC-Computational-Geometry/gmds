//
// Created by rochec on 09/05/2022.
//

#ifndef GMDS_SMOOTHINGPAVING_2D_H
#define GMDS_SMOOTHINGPAVING_2D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_AERO_export.h"
#include "gmds/ig/Mesh.h"
#include <gmds/aero/Front.h>
#include <gmds/aero/AeroBoundaries_2D.h>
#include <gmds/aero/Params.h>
#include <string>
#include <map>
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_AERO_API SmoothingPaving_2D
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
         *  @param[in] AMesh the mesh where we work on
         *
	 */
	SmoothingPaving_2D(Mesh *AMesh, Front AFront);

	/*-------------------------------------------------------------------*/
	/** @brief Execute the algorithm
	 */
	STATUS execute();
	/*-------------------------------------------------------------------*/

 private:
	/*-------------------------------------------------------------------*/
	/** @brief Met à jour la map des anciennes coordonnées
	 */
	void UpdateOldCoords();
	/*-------------------------------------------------------------------*/
	/** @brief Initialise la map des longueurs idéales ld pour les noeuds
	 * du front
	 */
	void InitializeMapLd();
	/*-------------------------------------------------------------------*/
	/** @brief Lisse les noeuds du front
	 */
	void FrontNodeSmoothing();
	/*-------------------------------------------------------------------*/
	/** @brief Calcule Delta a
	 		*  @param[in] n_id noeud du front concerné
	 */
	math::Vector3d ComputeDa(TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief Calcule Delta b
	 		*  @param[in] n_id noeud du front concerné
	 */
	math::Vector3d ComputeDb(TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief Calcule Delta c
	 		*  @param[in] n_id noeud du front concerné
	 */
	math::Vector3d ComputeDc(TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief Lisse les noeuds internes (tous sauf paroi et front)
	 */
	void InteriorNodeSmoothing();
	/*-------------------------------------------------------------------*/

 private:
	/** mesh we work on */
	Mesh *m_mesh;
	/** Front */
	Front m_Front;
	/** Ancienne coordonnées */
	std::map<TCellID,math::Point> m_map_old_coords;
	/** Map des tailles suivant "y" */
	std::map<TCellID,double> m_map_ld;
	/** Marque les noeuds du front */
	int m_markFrontNodes;

};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_SMOOTHINGPAVING_2D_H
