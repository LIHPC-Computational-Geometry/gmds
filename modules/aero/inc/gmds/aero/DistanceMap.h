//
// Created by rochec on 14/01/2022.
//

#ifndef GMDS_DISTANCEMAP_H
#define GMDS_DISTANCEMAP_H

/*----------------------------------------------------------------------------*/
#include "GMDSAero_export.h"
#include <gmds/ig/Mesh.h>
#include <iostream>
#include <string>
#include <map>
#include <fstream>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
class GMDSAero_API DistanceMap {
 public:
	/*--------------------------------------------------------------------*/
	/** @enum  Status code for executing algorithms
	 */
	typedef enum {
		FAIL,
		SUCCESS
	} STATUS;

	/*-------------------------------------------------------------------*/
	/** \brief Default constructor
	 */
	DistanceMap();
	/*-------------------------------------------------------------------*/

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
	/** @brief Donne le nombre d'éléments dans le vecteur associé à la clé
	 */
	void getNbrIds(double v0, int &nbr);
	/*-------------------------------------------------------------------*/
	/** @brief Donne le nombre d'éléments dans le tas
	 */
	void getNbrIdsTotal(int &nbr);
	/*-------------------------------------------------------------------*/
	/** @brief Dit si la map est vide ou non
	 */
	bool isEmpty();
	/*-------------------------------------------------------------------*/
	/** @brief Surcharge de l'opérateur << pour l'affichaque des maps
	 */
	friend GMDSAero_API std::ostream& operator<<(std::ostream&, const DistanceMap&);
	/*-------------------------------------------------------------------*/
	/** @brief Surcharge de l'opérateur ()
	 */
	std::vector<TCellID> operator()(const double AI) {return m_map[AI];}
	/*-------------------------------------------------------------------*/
 private:
	/** Tas des couples (distance provisoire, liste d'ids) */
	std::map<double, std::vector<TCellID>> m_map;


};
/*----------------------------------------------------------------------------*/
}     // namespace gmds

#endif     // GMDS_DISTANCEMAP_H
