//
// Created by rochec on 09/02/2022.
//

#ifndef GMDS_ABSTRACTAEROPIPELINE_H
#define GMDS_ABSTRACTAEROPIPELINE_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/claire/Params.h>
#include <string>
#include <map>
#include <fstream>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*---------------------------------------------------------------------------*/
// Frame File Headers
/*----------------------------------------------------------------------------*/
/** \class  AbstractAeroPipeline
 *  \brief  Classe abstraite pour l'algorithme de maillage
 *  			pour l'aéro. Les classes filles sont pour le
 *  			cas 2D et le cas 3D.
 */
class LIB_GMDS_CLAIRE_API AbstractAeroPipeline{

 public:
	/** @brief Constructor.
         *  @param
	 */
	AbstractAeroPipeline();

	/*-------------------------------------------------------------------*/
	/** \brief Function to be called for mesh generation
	 */
	virtual void execute()=0;
	/*-------------------------------------------------------------------*/

 protected:
	/** mesh we work on */
	Mesh* m_mesh;
	/** parameters for the algorithm */
	ParamsAero m_params;
	/** Marque sur les noeuds de la paroi */
	int m_markFrontNodesParoi;
	/** Marque sur les noeuds de la frontière ext */
	int m_markFrontNodesExt;

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_ABSTRACTAEROPIPELINE_H
