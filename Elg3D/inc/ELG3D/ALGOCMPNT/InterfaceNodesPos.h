/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    InterfaceNodesPos.h
 *  \author  legoff
 *  \date    03/26/2018
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_INTERFACENODESPOS_H_
#define ELG3D_INTERFACENODESPOS_H_
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
#include <KM/DS/Mesh.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include "ELG3D/ALGOCMPNT/MaterialGradientComputation.h"
#include "ELG3D/DATACMPNT/FracPres.h"
#include "ELG3D/DATACMPNT/MaterialAssignment.h"
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/

    const kmds::TFloat InterfaceNodesPos_ISOVALUEFP = 1./2.;


/*----------------------------------------------------------------------------*/

    /*------------------------------------------------------------------------*/
    /** \brief  Computes the Pcross position between two cells where the iso-contour
     *          for one material should be.
     *
     *  \param[in]  Acid0 the id of the first cell
     *  \param[in]  Acid1 the id of the second cell
     *  \param[in]  AMat the material id
     *  \param[in]  Afp the fracpres
     *  \param[in]  AVarCenter a variable carried by the cells storing the centroids
     *  \param[out]  APcross the computed position
     *
     *  \return true if the position was computed as a linear interpolation, false
     *          if the iso-value InterfaceNodesPos_ISOVALUEFP was outside the range
     *          [fracpres(AMat, Acid0) .. fracpres(AMat, Acid1)] in which case
     *          APcross is either one end-points or right at the middle in case of equality
     */
    bool InterfaceNodesPos_computePcross_oneface_onemat_3D(const kmds::TCellID Acid0,
                                         const kmds::TCellID Acid1,
                                         const int AMat,
                                         const elg3d::FracPres* Afp,
                                         const kmds::Variable<gmds::math::Point>* AVarCenter,
                                         gmds::math::Point& APcross);


    /*------------------------------------------------------------------------*/
    /** \brief  Computes a new position for the interface nodes; also computes
     *          the normal for each material
     *
     *  \param[in]  ASelection a container of the nodes on the interfaces
     *  \param[in]  AMesh the mesh
     *  \param[in]  Ama the material assignment for the cells
     *  \param[in]  Afp the fracpres
     *  \param[in]  AVarCenter a variable carried by the cells storing the centroids
    // *  \param[out]  APcross the computed position
     *
     */
    void InterfaceNodesPos_computeNodesNewPos_2D(const kmds::GrowingView<kmds::TCellID>* ASelectionInterfaceNodes,
                                                 kmds::Mesh* AMesh,
                                                 const kmds::Connectivity* c_N2F,
                                                 const kmds::Connectivity* c_F2C,
                                                 const MaterialAssignment* Ama,
                                                 const elg3d::FracPres* Afp,
                                                 const kmds::Variable<gmds::math::Point>* AVarCenter,
                                                 const kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads,
                                                 kmds::Variable<gmds::math::Point>* AVarNewPos);

    /*------------------------------------------------------------------------*/
    /** \brief  Computes a new position for the interface nodes; also computes
     *          the normal for each material
     *
     *  \param[in]  ASelection a container of the nodes on the interfaces
     *  \param[in]  AMesh the mesh
     *  \param[in]  Ama the material assignment for the cells
     *  \param[in]  Afp the fracpres
     *  \param[in]  AVarCenter a variable carried by the cells storing the centroids
     *  \param[out]  APcross the computed position
     *
     */
    void InterfaceNodesPos_computeNodesNewPos_3D(const kmds::GrowingView<kmds::TCellID>* ASelectionInterfaceNodes,
                                                 kmds::Mesh* AMesh,
                                                 const kmds::Connectivity* c_N2F,
                                                 const kmds::Connectivity* c_F2C,
                                                 const MaterialAssignment* Ama,
                                                 const elg3d::FracPres* Afp,
                                                 const kmds::Variable<gmds::math::Point>* AVarMidpoints,
                                                 const kmds::Variable<MaterialGradientComputation_midcellGradients>* AVarGrads,
                                                 kmds::Variable<gmds::math::Point>* AVarNewPos);

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_INTERFACENODESPOS_H_ */
