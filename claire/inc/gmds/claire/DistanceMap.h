//
// Created by rochec on 14/01/2022.
//

#ifndef GMDS_DISTANCEMAP_H
#define GMDS_DISTANCEMAP_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
#include <string>
#include <map>
/*----------------------------------------------------------------------------*/
#include <fstream>
namespace gmds {
/*----------------------------------------------------------------------------*/
class LIB_GMDS_CLAIRE_API DistanceMap {
 public:
	/*--------------------------------------------------------------------*/
	/** @enum  Status code for executing algorithms
	 */
	typedef enum {
		FAIL,
		SUCCESS
	} STATUS;

	/*-------------------------------------------------------------------*/

	DistanceMap();

 public:

	/*-------------------------------------------------------------------*/
	/** @brief Ajoute un élément dans la map
	 */
	void add(double v0, TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief Récupère le premier élément de la map et l'enlève de celle ci
	 */
	void getAndRemoveFirst(double &AMinDist, TCellID &AMinId);
	/*-------------------------------------------------------------------*/
	/** @brief Récupère un élément dans la map et l'enlève de celle ci
	 */
	void remove(double v0, TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief Met à jour la clé d'un id dans la map
	 */
	void update(double v0_old, double v0_new, TCellID n_id);
	/*-------------------------------------------------------------------*/
	/** @brief Vérifie que chaque id est bien présent une seule fois
	 */
	bool check();
	/*-------------------------------------------------------------------*/
	friend  std::ostream& operator<<(std::ostream&, const DistanceMap&);

	std::vector<TCellID> operator()(const double AI) {return m_map[AI];}
 private:
	/** Tas des couples (distance provisoire, liste d'ids) */
	std::map<double, std::vector<TCellID>> m_map;


};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_DISTANCEMAP_H
