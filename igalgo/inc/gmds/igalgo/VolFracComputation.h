/*----------------------------------------------------------------------------*/
#ifndef GMDS_VOL_FRAC_COMPUTATION_H
#define GMDS_VOL_FRAC_COMPUTATION_H
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include "GMDSIgAlgo_export.h"
/*----------------------------------------------------------------------------*/
//#include <map>
//#include <vector>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** @brief Computes the volume fractions of mesh AImprintMesh
      *
      * @param Amesh the mesh to work on
      * @param AImprintMesh, the mesh to imprint onto the first mesh
      * @param AVolFrac the variable that will carry the volume fractions 
 */
void GMDSIgAlgo_API volfraccomputation_2d(Mesh* AMesh, const Mesh* AImprintMesh, gmds::Variable<double>* AVolFrac);
/*----------------------------------------------------------------------------*/
//class GMDSIgAlgo_API GridBuilder
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_VOL_FRAC_COMPUTATION_H
