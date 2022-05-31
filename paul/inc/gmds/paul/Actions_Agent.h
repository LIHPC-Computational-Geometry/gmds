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

	class Actions
   {
	 public:
	   Actions(GridBuilderAround* AGrid);
	   /** @brief action delete face
	    *
	    * @param faceID id of the face to delete
	    */
	   void executeDeleteFace(int faceID);

	   /** \brief action for cute a face
	    *
	    * @param firstNodeID first node
	    * @param secondNodeID second node
	    */
	   void executeCutEdge(Node firstNodeID,Node secondNodeID);
		/** \brief cut the face with a direction indicator
		 *
		 * @param AFace the face cut
		 * @param direction if direction = 0 cut horizontal, if direction = 1 cut vertical
		 */
	   void executeCutFace(Face AFace, int direction);

	   /** \brief glide a node to the boundary of the target mesh
	    *
	    * @param node The node to drag
	    * @param AMesh the target mesh
	    */

	   void executeGlideNode(Node node,Mesh *AMesh);

	   void executeGlideMaxNodeFace(Face AFace);
	   void executeGlideMinNodeFace(Face AFace);

	   void executeAction(int actionSelect,Face AFace);

	   gmds::GridBuilderAround g_grid;
	   gmds::Tools tool;
   };
}
#endif     // GMDS_ACTIONS_AGENT_H
