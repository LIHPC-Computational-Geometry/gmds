#ifndef GMDS_TOOLS_H
#define GMDS_TOOLS_H
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_BAPTISTE_export.h"
#include <gmds/ig/Mesh.h>
#include "GMDSIgAlgo_export.h"
/*----------------------------------------------------------------------------*/

namespace gmds
{
	std::vector<double> LinearSpacedArray(double a, double b, std::size_t N);

   // Mesh readMesh(std::string filename);



   Variable<double>* getVolFrac(Mesh mesh);

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