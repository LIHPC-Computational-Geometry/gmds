//
// Created by rochec on 10/02/2022.
//

#ifndef GMDS_LEVELSETELOI_H
#define GMDS_LEVELSETELOI_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_AERO_export.h"
#include <gmds/aero/AbstractLevelSet.h>
/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
/** \class  LevelSetEloi
 *  \brief
 */
class LIB_GMDS_AERO_API LevelSetEloi: public AbstractLevelSet {
 public:

	/*------------------------------------------------------------------------*/
	/** \brief Default constructor
	*	 @param AMesh pointeur vers un maillage
	*	 @param AmarkFrontNodes marque sur les noeuds sur front
	*	 @param Adistance variable distance
	 */
	LevelSetEloi(Mesh *AMesh, TInt AmarkFrontNodes, Variable<double> *Adistance);

	/*------------------------------------------------------------------------*/

 protected:
	/** \brief Construit un vecteur des noeuds du voisinage
	 *
	 * @param n
	 * @return a vector of nodes
	 */
	std::vector<Node> getNeighbors(Node n) override;

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_LEVELSETELOI_H
