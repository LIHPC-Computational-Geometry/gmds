/*----------------------------------------------------------------------------*/
#ifndef GMDS_BLOCKCLASSIFICATOR_H
#define GMDS_BLOCKCLASSIFICATOR_H
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/

#include "LIB_GMDS_BLOCK_MESHER_export.h"
#include "gmds/cadfac/FACManager.h"
#include "gmds/cad/GeomMeshLinker.h"
#include "gmds/ig/Mesh.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
	class LIB_GMDS_BLOCK_MESHER_API BlockClassificator{
	 public:
	   BlockClassificator(Mesh* ABlocks, cad::GeomMeshLinker* ALinker, cad::FACManager* AManager);

	   ~BlockClassificator();

	   void blockCreation(int AMTrix[5][5][5]);

	   void blockClassification();

	   Mesh *getBlocks();


	 private:
	   /* the block structure*/
	   Mesh* m_blocks;
	   /* the linker between blocks and geometry*/
	   cad::GeomMeshLinker* m_linker;
	   /* the geometry model*/
	   cad::FACManager* m_manager;

	   void getEntityToClassify(math::Point APoint, int &ADim, int &AID);
	   void classifyBlockEntity(TCellID ACellID, int ACellDim, int ADim, int AID);
   };
}

#endif     // GMDS_BLOCKCLASSIFICATOR_H
