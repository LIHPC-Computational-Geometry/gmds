/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    BoundingBoxGeomAssociation.h
 *  \author  legoff
 *  \date    04/23/2018
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_BOUNDINGBOXGEOMASSOCIATION_H_
#define ELG3D_BOUNDINGBOXGEOMASSOCIATION_H_
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
#include <gmds/cad/FACManager.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/

    /*------------------------------------------------------------------------*/
    /** \brief  Creates the faceted geom model from the bounding box of the mesh and
     *          associates the nodes to the corresponding geom entities.
     *
     *          WARNING : works only on rectangular domains
     *
     *  \param[in]  AMesh the mesh
     *  \param[out]  AAFacGeomManager the faceted geom manager
     *  \param[out]  AVarNodeGeomAssociation
     */
    void BoundingBoxGeomAssociation_init_2D(kmds::Mesh* AMesh,
                                            gmds::cad::FACManager* AFacGeomManager,
                                            kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation);

    /*------------------------------------------------------------------------*/
    /** \brief  Creates the faceted geom model from the bounding box of the mesh and
     *          associates the nodes to the corresponding geom entities.
     *
     *          WARNING : works only on rectangular domains
     *
     *  \param[in]  AMesh the mesh
     *  \param[out]  AAFacGeomManager the faceted geom manager
     *  \param[out]  AVarNodeGeomAssociation
     */
    void BoundingBoxGeomAssociation_init_3D(kmds::Mesh* AMesh,
                                            gmds::cad::FACManager* AFacGeomManager,
                                            kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation);

    void BoundingBoxGeomAssociation_assoc_3D(kmds::Mesh* AMesh,
                                             const double minXYZ[3],
                                             const double maxXYZ[3],
                                             const gmds::cad::FACManager* AFacGeomManager,
                                             kmds::Variable<std::uintptr_t>* AVarNodeGeomAssociation);

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_BOUNDINGBOXGEOMASSOCIATION_H_ */
