//
// Created by Paul Bourmaud on 07/04/2022.
//
#ifndef GMDS_ACTIONS_AGENT_H
#define GMDS_ACTIONS_AGENT_H
/*----------------------------------------------------------------------------*/
#include "gmds/ig/Mesh.h"
#include "GMDSIgAlgo_export.h"
#include "gmds/paul/Grid.h"
/*----------------------------------------------------------------------------*/

namespace gmds {
	// Action Delete Face (swap value activate variable)
	class Actions
   {
	 public:
	   Actions(GridBuilderAround* AGrid);
	   /** @brief action delete face
	    */
	   void executeDeleteFace(const int faceID);

	   void executeCutEdge(gmds::Edge edge);

	   void executeGlideNode(gmds::Node node);

	   gmds::GridBuilderAround g_grid;
   };
}
#endif     // GMDS_ACTIONS_AGENT_H
