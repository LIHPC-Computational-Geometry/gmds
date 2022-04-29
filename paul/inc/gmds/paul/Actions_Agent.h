//
// Created by Paul Bourmaud on 07/04/2022.
//
#ifndef GMDS_ACTIONS_AGENT_H
#define GMDS_ACTIONS_AGENT_H
/*----------------------------------------------------------------------------*/
#include "gmds/ig/Mesh.h"
#include "GMDSIgAlgo_export.h"
#include "gmds/paul/Grid.h"
#include "gmds/paul/Tools.h"
/*----------------------------------------------------------------------------*/

namespace gmds {
	// Action Delete Face (swap value activate variable)
	class Actions
   {
	 public:
	   Actions(GridBuilderAround* AGrid);
	   /** @brief action delete face
	    *
	    * @param faceID id of the face to delete
	    */
	   void executeDeleteFace(const int faceID);

	   /** \brief action for cute a face
	    *
	    * @param firstNodeID first node
	    * @param secondNodeID second node
	    */
	   void executeCutEdge(Node firstNodeID,Node secondNodeID);

	   void executeGlideNode(gmds::Node node);

	   gmds::GridBuilderAround g_grid;
	   gmds::Tools tool;
   };
}
#endif     // GMDS_ACTIONS_AGENT_H
