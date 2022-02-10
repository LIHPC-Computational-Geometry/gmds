//
// Created by rochec on 10/02/2022.
//

#ifndef GMDS_LEVELSETELOI_H
#define GMDS_LEVELSETELOI_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/claire/AbstractLevelSet.h>
/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
/** \class  LevelSetEloi
 *  \brief
 */
class LIB_GMDS_CLAIRE_API LevelSetEloi: public AbstractLevelSet {
 public:

	/*------------------------------------------------------------------------*/
	/** \brief Default constructor
	 */
	LevelSetEloi(Mesh *AMesh, int AmarkFrontNodes, Variable<double>* Adistance);

	/*------------------------------------------------------------------------*/

 protected:
	/**
	 *
	 * @param n
	 * @return
	 */
	virtual std::vector<Node> getNeighbors(Node n);

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_LEVELSETELOI_H
