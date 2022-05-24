//
// Created by bourmaudp on 10/05/22.
//

#ifndef GMDS_ENVIRONMENT_H
#define GMDS_ENVIRONMENT_H

/*----------------------------------------------------------------------------*/
#include "gmds/igalgo/VolFracComputation.h"
#include "gmds/ig/Mesh.h"
#include "gmds/paul/Grid.h"
/*----------------------------------------------------------------------------*/

namespace gmds {
	class Environment
   {
	 public :
	   Environment(GridBuilderAround* AGrid,Mesh* AMesh);

	   double globalIoU();

	   double localIoU(const Face AFaceID);

	   Face faceSelect();

	   gmds::GridBuilderAround g_grid;
	   gmds::Mesh m_mesh;
   };
}
#endif     // GMDS_ENVIRONMENT_H
