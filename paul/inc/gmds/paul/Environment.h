//
// Created by bourmaudp on 10/05/22.
//

#ifndef GMDS_ENVIRONMENT_H
#define GMDS_ENVIRONMENT_H

/*----------------------------------------------------------------------------*/
#include "gmds/igalgo/VolFracComputation.h"
#include "gmds/ig/Mesh.h"
#include "gmds/paul/Grid.h"
#include "gmds/paul/Actions_Agent.h"
/*----------------------------------------------------------------------------*/

namespace gmds {
	class Environment
   {
	 public :
	   Environment(GridBuilderAround* AGrid,Mesh* AMeshTarget,Actions* AAction);

	   //virtual ~Environment();

	   /** \brief calc the global Intersection over Union
	    *
	    * @return the IoU global
	    */
	   double globalIoU();

	   void deleteEnv();

	   /** \brief calc the local Intersection over Union for a specify face
	    *
	    * @param AFaceID the face
	    * @return the local IoU
	    */
	   double localIoU(const Face AFace);

	   double reward(const Face AFace);


		/** \brief execute the action
		 *
		 * @param AFace the face select by the environment
		 * @param numberAction 0 = executeDeleteFace ; 1 = executeCutFace vertical; 2 = executeCutFace horizontal; 3 = Glide point max range; 4 = glide point min range
		 */
	   void executeAction(Face AFace, int numberAction);

	   /** \brief select a activate face with the littlest local IoU
	    *
	    * @return the selected face
	    */
	   Face faceSelect();

	   void calcVolFrac();

	   gmds::GridBuilderAround g_grid;
	   gmds::Mesh m_mesh;

	   gmds::Actions action;
   };
}
#endif     // GMDS_ENVIRONMENT_H
