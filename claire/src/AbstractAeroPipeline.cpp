//
// Created by rochec on 09/02/2022.
//

/*------------------------------------------------------------------------*/
#include <gmds/claire/AbstractAeroPipeline.h>
/*------------------------------------------------------------------------*/
using namespace gmds;
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
AbstractAeroPipeline::AbstractAeroPipeline(ParamsAero Aparams) :
  m_params(Aparams),
  m_meshTet(nullptr),
  m_meshHex(nullptr),
  m_manager(new cad::FACManager()),
  m_linker_TG(new cad::GeomMeshLinker()),
  m_linker_HG(new cad::GeomMeshLinker())
{

}
/*------------------------------------------------------------------------*/


/*------------------------------------------------------------------------*/
AbstractAeroPipeline::~AbstractAeroPipeline(){
	if(m_manager){
		delete m_manager;
	}
	if(m_linker_TG){
		delete m_linker_TG;
	}
	if(m_linker_HG){
		delete m_linker_HG;
	}
	if(m_meshTet){
		delete m_meshTet;
	}
	if(m_meshHex){
		delete m_meshHex;
	}
}
/*------------------------------------------------------------------------*/