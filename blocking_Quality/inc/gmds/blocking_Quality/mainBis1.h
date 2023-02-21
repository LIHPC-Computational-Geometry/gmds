//
// Created by bourmaudp on 20/12/22.
//

#ifndef GMDS_MAINBIS1_H
#define GMDS_MAINBIS1_H
/*----------------------------------------------------------------------------------------*/

#include <cstdlib>
#include <iostream>
#include <ctime>
#include <string>


#include <gmds/ig/Mesh.h>
#include <gmds/ig/Node.h>
#include <gmds/quality/QuadQuality.h>
#include <gmds/quality/HexQuality.h>
#include <gmds/blocking_Quality/BlockingQuality.h>
#include <gmds/igalgo/GridBuilder.h>
#include <gmds/blockMesher/BlockMesher.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/blockMesher/BlockMesher.h>
#include <gmds/cadfac/FACManager.h>
#include <gmds/igalgo/BoundaryOperator.h>
#include <gmds/smoothy/LaplacianSmoother.h>

#include <gmds/io/IGMeshIOService.h>
#include <gmds/io/VTKWriter.h>
#include <gmds/io/VTKReader.h>
/*----------------------------------------------------------------------------------------*/
namespace gmds {


void fonctionMain1();
void fonctionMain2();

void test();

}
/*--------------------------------------------------------------------------------------	--*/

#endif     // GMDS_MAINBIS1_H
