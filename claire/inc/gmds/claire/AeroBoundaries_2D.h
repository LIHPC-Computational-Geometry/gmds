//
// Created by rochec on 23/03/2022.
//

#ifndef GMDS_AEROBOUNDARIES_2D_H
#define GMDS_AEROBOUNDARIES_2D_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/claire/AbstractAeroBoundaries.h>
/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
/** \class  AeroBoundaries_2D
 *  \brief  Caractéristiques des frontières du maillage 2D.
 */
class LIB_GMDS_CLAIRE_API AeroBoundaries_2D: public AbstractAeroBoundaries {
 public:

	/*------------------------------------------------------------------------*/
	/** \brief Default constructor
	 */
	explicit AeroBoundaries_2D(Mesh *AMesh);

	/*------------------------------------------------------------------------*/
	/** @brief Retourne un vecteur ordonné des id des noeuds du bord
	 * correspondant à la couleur color.
	 */
	std::vector<TCellID> BndNodesOrdered(int color);
	/*------------------------------------------------------------------------*/
	/** @brief Retourne la longueur d'un bord (périmètre) à partir de la couleur
	 * color du bord.
	 */
	double ComputeBoundaryLength(int color);
	/*-----------------------------------------------------------------------*/
	/** @brief Retourne un vecteur des id des arêtes du bord
	 * correspondant à la couleur color.
	 */
	std::vector<TCellID> BndEdges(int color);
	/*------------------------------------------------------------------------*/

 protected:
	/*-------------------------------------------------------------------*/
	/** \brief Marque les noeuds sur les bords dans la marque
	 * m_markBoundaryNodes
	 */
	void MarkBoundariesNodes() override;
	/*-------------------------------------------------------------------*/
	/** \brief Identifie de quelle couleur est le front Amont dans la
	 * variable m_var_color_bords.
	 */
	void WhichColorIsAmont() override;
	/*-------------------------------------------------------------------*/

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_AEROBOUNDARIES_2D_H
