//
// Created by rochec on 10/02/2022.
//

#ifndef GMDS_LEVELSETEXTENDED_H
#define GMDS_LEVELSETEXTENDED_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/claire/AbstractLevelSet.h>
/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
/** \class  LevelSetExtended
 *  \brief
 */
class LIB_GMDS_CLAIRE_API LevelSetExtended: public AbstractLevelSet {
 public:

	/*------------------------------------------------------------------------*/
	/** \brief Default constructor
	 */
	LevelSetExtended(Mesh *AMesh, int AmarkFrontNodes, Variable<double>* Adistance);

	/*------------------------------------------------------------------------*/

 protected:
	/** \brief Construit un vecteur des noeuds du voisinage
	 *
	 * @param n
	 * @return
	 */
	virtual std::vector<Node> getNeighbors(Node n);

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_LEVELSETEXTENDED_H
