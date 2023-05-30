//
// Created by bourmaudp on 04/01/23.
//

#ifndef GMDS_LINKERBLOCKINGGEOM_H
#define GMDS_LINKERBLOCKINGGEOM_H

#include "gmds/cad/GeomMeshLinker.h"
#include "gmds/ig/Mesh.h"
#include <gmds/ig/Node.h>

#include <gmds/blockMesher/BlockMesher.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/cad/GeomPoint.h>
#include <gmds/cad/GeomCurve.h>
#include <gmds/cad/GeomSurface.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/smoothy/LaplacianSmoother.h>


#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
#include <gmds/ig/MeshDoctor.h>

namespace gmds{
void reverseLink(Mesh* ABlocks, Mesh* AGeom, cad::GeomMeshLinker* ALinker,
                 std::map<TCellID ,TCellID > mapNodes2Points,std::map<TCellID ,TCellID > mapEdges2Curves,
                 std::map<TCellID ,TCellID > mapFaces2Surfaces,std::map<TCellID ,TCellID > mapRegions2Volumes);




class LinkerBlockingGeom{
 public :
	LinkerBlockingGeom(Mesh* ABlocks, cad::FACManager* AGeom);

	virtual ~LinkerBlockingGeom();

	void execute(cad::GeomMeshLinker* ALinker);

 private :
	/** the block structure*/
	Mesh* m_blocks;
	/** the geometry mesh*/
	cad::FACManager* m_geom;

};
}
#endif     // GMDS_LINKERBLOCKINGGEOM_H
