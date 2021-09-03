/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * This software is a computer program whose purpose is to provide a set of
 * functionnalities to represent and handle any type of meshes (2D, 3D,
 * triangles, tetrahedra, quad, hexa, polygons, polyhedra, etc.) and write
 * meshing algorithms. So it gathers many mathematical objects like points,
 * segment, quaternions, etc. and basic algorithms useful to build more evolved
 * ones.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and, more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
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
