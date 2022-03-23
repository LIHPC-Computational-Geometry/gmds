//
// Created by rochec on 23/03/2022.
//

#ifndef GMDS_ABSTRACTAEROBOUNDARIES_H
#define GMDS_ABSTRACTAEROBOUNDARIES_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
#include <string>
#include <map>
#include <fstream>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*---------------------------------------------------------------------------*/
// Frame File Headers
/*----------------------------------------------------------------------------*/
/** \class  AbstractAeroBoundaries
 *  \brief  Classe abstraite pour les informations relatives aux frontières
 *  dans le cadre de l'aéro.
 */
class LIB_GMDS_CLAIRE_API AbstractAeroBoundaries{

 public:
	/*--------------------------------------------------------------------*/
	/** @enum  Status code for executing algorithms
	 */
	typedef enum {
		FAIL,
		SUCCESS
	} STATUS;
	/*--------------------------------------------------------------------*/
	/** @brief Constructor.
         *  @param
	 */
	AbstractAeroBoundaries(Mesh *AMesh);
	/*-------------------------------------------------------------------*/
	/** \brief Function to be called for mesh generation
	 */
	virtual AbstractAeroBoundaries::STATUS execute()=0;
	/*-------------------------------------------------------------------*/

 protected:
	/** pointer to mesh we work on */
	Mesh* m_mesh;
	/** Nombre de bords distincts de bords */
	int m_nbrBords;
	/** Nombre de bords distincts sur la paroi */
	int m_nbrBordsParoi;
	/** Couleur qui correspond au bord extérieur */
	int m_bigest_color;
	/** Colorie chaque bord d'une couleur != */
	Variable<int>* m_var_color_bords ;
	/** Marque sur les noeuds de la paroi */
	int m_markNodesParoi;
	/** Marque sur les noeuds de la frontière amont */
	int m_markNodesAmont;

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_ABSTRACTAEROBOUNDARIES_H
