#ifndef GMDS_CAPSULE_H
#define GMDS_CAPSULE_H
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_BAPTISTE_export.h"
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/

namespace gmds
{
	class Capsule
   {
	 public:
	   Capsule();

	   ~Capsule();

	   Mesh m_mesh;

	   void readMesh(std::string filename);

	   void saveMesh(std::string filename);
   };
}
#endif     // GMDS_CAPSULE_H