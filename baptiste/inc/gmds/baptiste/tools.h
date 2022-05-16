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
}
#endif     // GMDS_TOOLS_H