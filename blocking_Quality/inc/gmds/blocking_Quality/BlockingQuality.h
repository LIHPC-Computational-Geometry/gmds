//
// Created by bourmaudp on 02/12/22.
//
/*----------------------------------------------------------------------------------------*/
#ifndef GMDS_BLOCKINGQUALITY_H
#define GMDS_BLOCKINGQUALITY_H
/*----------------------------------------------------------------------------------------*/
#include "gmds/cad/GeomMeshLinker.h"
#include <gmds/ig/Mesh.h>
#include <gmds/igalgo/BoundaryExtractor3D.h>
#include <gmds/math/Matrix.h>
#include <gmds/quality/HexQuality.h>
/*----------------------------------------------------------------------------------------*/
namespace gmds {
/*--------------------------------------------------------------------------------------	--*/

double blockQuality(Region ARegion);

double allBlocksQuality(Mesh* AMesh);

double layeringNodes(Mesh* AMesh, const Mesh* AImprintMesh);

double internEdgeValence(Mesh* AMesh);

double boundaryEdgeValence(Mesh* AMesh);

double globalEdgeValence(Mesh* AMesh);

double ratioGlobalEdgeValence(Mesh* AMesh);

int blockingSimplicity(Mesh* AMesh);

double blockingQuality(Mesh* AMesh, Mesh* AGeometry);

std::unordered_map<TCellID ,std::vector<TCellID>> mapNodes2Points(Mesh* ABlocks, Mesh* AGeom, cad::GeomMeshLinker* ALinker);
std::unordered_map<TCellID ,std::vector<TCellID>> mapEdges2Curves(Mesh* ABlocks, Mesh* AGeom, cad::GeomMeshLinker* ALinker);
std::unordered_map<TCellID ,std::vector<TCellID>> mapFaces2Surfaces(Mesh* ABlocks, Mesh* AGeom, cad::GeomMeshLinker* ALinker);


bool linkNodes2Points(Mesh* ABlocks, Mesh* AGeom, cad::GeomMeshLinker* ALinker);
bool linkEdges2Curves(Mesh* ABlocks, Mesh* AGeom, cad::GeomMeshLinker* ALinker);
bool linkFaces2Surfaces(Mesh* ABlocks, Mesh* AGeom, cad::GeomMeshLinker* ALinker);


std::vector<Node> noLinkNodes2Points(Mesh* ABlocks, Mesh* AGeom, cad::GeomMeshLinker* ALinker);



}
#endif     // GMDS_BLOCKINGQUALITY_H
