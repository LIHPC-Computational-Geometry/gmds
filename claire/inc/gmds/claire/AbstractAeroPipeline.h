//
// Created by rochec on 09/02/2022.
//

#ifndef GMDS_ABSTRACTAEROPIPELINE_H
#define GMDS_ABSTRACTAEROPIPELINE_H

/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
#include <gmds/claire/Params.h>
#include <gmds/claire/AbstractAeroBoundaries.h>
#include <gmds/cadfac/FACManager.h>
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
	explicit AbstractAeroPipeline(std::string &Aparams, std::string &Aworking_dir);
	/*--------------------------------------------------------------------*/
	/** @brief Destructor.
         *  @param
	 */
	~AbstractAeroPipeline();
	/*-------------------------------------------------------------------*/
	/** \brief Function to be called for mesh generation
	 */
	virtual AbstractAeroPipeline::STATUS execute()=0;
	/*-------------------------------------------------------------------*/

 protected:
	/** pointer to mesh we work on */
	Mesh* m_meshTet;
	/** pointer to mesh we generate */
	Mesh* m_meshHex;
	/** parameters for the algorithm */
	ParamsAero m_params;
	/** Variable sur le nouveau maillage, indique à quelle couche appartient un noeud */
	Variable<int>* m_couche_id;
	/** Manager */
	cad::FACManager* m_manager;
	/** Linker maillage T entrée à la géométrie */
	cad::GeomMeshLinker* m_linker_TG;
	/** Linker maillage Q/H sortie à la géométrie */
	cad::GeomMeshLinker* m_linker_HG;

};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif     // GMDS_ABSTRACTAEROPIPELINE_H
