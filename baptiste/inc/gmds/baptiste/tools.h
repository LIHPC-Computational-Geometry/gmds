#ifndef GMDS_TOOLS_H
#define GMDS_TOOLS_H
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_BAPTISTE_export.h"
#include <gmds/ig/Mesh.h>
#include "gmds/baptiste/RLBlockSet.h"
#include "GMDSIgAlgo_export.h"
/*----------------------------------------------------------------------------*/

namespace gmds
{
	std::vector<double> LinearSpacedArray(double a, double b, std::size_t N);

   void cloneMesh(const Mesh &originalMesh, Mesh &newMesh);

   void cloneBlockSet(const RLBlockSet &originalBlockSet, RLBlockSet &newBlockSet);

   Mesh readMesh(std::string filename);

   class Tools
   {
	 public:
	   Tools();

	   Mesh m_mesh;
	   void readMesh(std::string filename);
	   void computeVolFrac(Mesh AMesh, Mesh targetMesh);
   };

}
#endif     // GMDS_TOOLS_H