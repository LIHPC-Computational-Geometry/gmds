/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    InterfaceNodesPosSmoothVF.h
 *  \author  legoff
 *  \date    10/29/2018
 */
/*----------------------------------------------------------------------------*/
#ifndef ORTOOLS_INTERFACENODESPOSSMOOTHVF_H_
#define ORTOOLS_INTERFACENODESPOSSMOOTHVF_H_
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <KM/DS/Mesh.h>
#include <KM/Utils/Graph.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/DATACMPNT/FracPres.h"
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
using namespace elg3d;
namespace milp {
/*----------------------------------------------------------------------------*/


    void
    InterfaceNodesPosSmoothVF_assignORTOOLS_xD(const kmds::GrowingView<kmds::TCellID> *ACellsIDs_source,
                                           const kmds::GrowingView<kmds::TCellID> *ACellsIDs_target,
                                           const FracPres *Afp_source,
                                           MaterialAssignment *Ama_target,
                                           const kmds::Variable<kmds::TCellID> *AVarSubCells2oldCells,
                                           const kmds::Variable<kmds::TCellID> *AVarOldCells2firstSubCells,
                                           const kmds::Variable<bool> *AVarMixedCells_source,
                                           const kmds::Variable<kmds::TCellID> *AVarCells2vertices,
                                           const kmds::Graph *AGraph,
                                           const int ANbSubPixels);



/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ORTOOLS_INTERFACENODESPOSSMOOTHVF_H_ */
