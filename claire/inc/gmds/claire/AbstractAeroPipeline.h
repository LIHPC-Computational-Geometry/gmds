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
	//AbstractAeroPipeline(ParamsAero Aparams, gmds::Mesh && Am);
	AbstractAeroPipeline(ParamsAero Aparams);
	/*-------------------------------------------------------------------*/
	/** \brief Function to be called for mesh generation
	 */
	virtual void execute()=0;
	/*-------------------------------------------------------------------*/
	/** \brief Function pour savoir si le pipeline s'est terminé
	 */
	bool getIsOver();
	/*-------------------------------------------------------------------*/

 protected:
	/** pointer to mesh we work on */
	Mesh* m_mesh;
	/** pointer to mesh we generate */
	Mesh* m_meshGen;
	/** parameters for the algorithm */
	ParamsAero m_params;
	/** Marque sur les noeuds de la paroi */
	int m_markFrontNodesParoi;
	/** Marque sur les noeuds de la frontière ext */
	int m_markFrontNodesExt;
	/** Variable sur le nouveau maillage, indique à quelle couche appartient un noeud */
	Variable<int>* m_couche_id;
	/** Est-ce que le code s'est bien terminé */
	bool m_isOver ;

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_ABSTRACTAEROPIPELINE_H
