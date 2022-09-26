//
// Created by calderans on 05/07/22.
//

#ifndef GMDS_BLOCKMESHER2D_H
#define GMDS_BLOCKMESHER2D_H

/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_BLOCK_MESHER_export.h"
#include "gmds/cad/FACManager.h"
#include "gmds/cad/GeomMeshLinker.h"
#include "gmds/ig/Mesh.h"
#include "gmds/utils/Array.h"
/**
 * Classe temporaire créée pour les besoins de Valentin. Elle permet de reclassifier
 * une structure de blocs sur une géométrie et de la mailler par la suite. L'idée est
 * de passer en argument de la classe un maillage de bloc quad optimisé par
 * Valentin qui est semi classifié et correspondant à la géométrie donnée dans
 * AManager. On va déduire la classification des nouveaux sommet et arêtes de blocs
 * créés et mailler les blocs.
 *
 */


namespace gmds {
class LIB_GMDS_BLOCK_MESHER_API BlockMesher2D{

 public:
	BlockMesher2D(Mesh* ABlocks, cad::GeomMeshLinker* ALinker, cad::FACManager* AManager);

	virtual ~BlockMesher2D();

	void updateClassification();

	void projectNodes();

	bool executeMeshing();

	Mesh* getMesh();


 private:

	void classifyBoundary(std::vector<Edge> &AEdges);

 private:
	/* the block structure*/
	Mesh* m_blocks;

	/* the linker between blocks and geometry*/
	cad::GeomMeshLinker* m_linker;
	cad::GeomMeshLinker* m_linkerMesh;
	/* the geometry model*/
	cad::FACManager* m_manager;

	/* the generated quad mesh*/
	Mesh* m_mesh;

	/** variables that store the grids nodes for blocks entities*/
	std::map<TCellID ,TCellID> m_blockV_to_meshN;
	std::map<TCellID , std::vector<TCellID> > m_blockE_to_meshN;
	std::map<TCellID , std::vector<std::vector<TCellID> > > m_blockF_to_meshN;
	Variable<std::vector<std::vector<TCellID> >  >* m_face_grids;
	Variable<std::vector<TCellID> >* m_edge_grids;
	bool getEdgeFrom(const TCellID AN0, const TCellID AN1, const std::vector<Edge> &AEdges, Edge &AResult, bool &AInverted);
	bool meshFaces();
	bool meshEdges();
	void meshVertices();
};
}

#endif     // GMDS_BLOCKMESHER2D_H
