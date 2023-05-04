//
// Created by bourmaudp on 21/02/23.
//

#ifndef GMDS_VALIDBLOCKING_H
#define GMDS_VALIDBLOCKING_H

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

namespace gmds {

class ValidBlocking
{
 public:
	ValidBlocking(Mesh *ABlocks, cad::FACManager *AGeom, cad::GeomMeshLinker *ALinker, std::map<std::vector<TCellID>,int> *elementsNoClassified);

	virtual ~ValidBlocking();

	bool execute();

 private:
	bool checkValidNodes();
	bool checkValidityNode(TInt ANodeId);

	bool checkValidEdges();
	bool checkValidityEdge(TInt AEdgeId);
	std::vector<std::pair<int,TCellID>> elementsOnCurve(cad::GeomCurve *ACurve);
	bool checkSuiteOfElements(std::vector<std::pair<int,TCellID>> AElementsOnCurve);
	bool checkExtremityCurve(std::vector<Edge> AEdgesVector, std::vector<Node> ANodesVector, cad::GeomCurve *ACurve);

	bool checkValidFaces();
	bool checkValidityFace(TInt AFaceId);
	bool checkValidityElementsFace(TInt AFaceId);
	bool checkFacesLinked();



	//Part Geom Check
	bool checkGeomValidity();
	//Points Check and return list of points no classified
	bool checkPoints();
	bool checkAPoint(cad::GeomPoint *APoint);
	std::map<int,TCellID> returnNoValidsPoints();

	//Curves Check and return list of curves no classified
	bool checkCurves();
	bool checkACurve(cad::GeomCurve *ACurve);
	std::map<int,TCellID> returnNoValidsCurves();

	//Surfaces Check and return list of surfaces no classified
	bool checkSurfaces();
	bool checkASurface(cad::GeomSurface *ASurface);
	std::map<int,TCellID> returnNoValidsSurfaces();

	//Volumes Check and return list of volumes no classified
	bool checkVolumes();
	bool checkAVolume(cad::GeomVolume *AVolume);
	std::map<int,TCellID> returnNoValidsVolumes();


	/** the block structure*/
	Mesh *m_blocks;
	/** the geometry mesh*/
	cad::FACManager *m_geom;
	/** the linker*/
	cad::GeomMeshLinker *a_linker;
	/** list of elements no classified*/
	std::map<std::vector<TCellID>,int> *elements_No_Classified;
	;

};
}



#endif     // GMDS_VALIDBLOCKING_H
