//
// Created by Paul Bourmaud on 12/04/2022.
//

#ifndef GMDS_Tools_H
#define GMDS_Tools_H
/*----------------------------------------------------------------------------*/
#include "gmds/ig/Mesh.h"
#include "GMDSIgAlgo_export.h"
#include "gmds/paul/Grid.h"

/*----------------------------------------------------------------------------*/


namespace gmds{
	//Tools
	class Tools{
	 public:
	   Tools(GridBuilderAround *AGrid);

	   int getValueActivateFace(const int faceIDChecked);

	   bool checkExistEdge(const int i1, const int i2);

	   bool checkCommonFace(const int i1, const int i2);

	   bool checkFollowIdNode(const int i1, const int i2);

	   std::vector<Node> getListNodesOfFace(const int faceID);

	   std::vector<Face> getListFacesOfNode(const int nodeID);

	   std::vector<Face> getFaceCommon(const int i1, const int i2);

	   void getIdPreviousNode(const int idNode, const int idFaceNode);

	   void getIdNextNode(const int idNode, const int idFaceNode);

	   gmds::GridBuilderAround g_grid;

   };
}

#endif     // GMDS_Tools_H
