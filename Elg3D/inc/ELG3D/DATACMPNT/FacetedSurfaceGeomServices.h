/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    FacetedSurfaceGeomServices.h
 *  \author  legoff
 *  \date    04/09/2019
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_FACETEDSURFACEGEOMSERVICES_H_
#define ELG3D_FACETEDSURFACEGEOMSERVICES_H_
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <map>
#include <set>
#include <string>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
#include <Kokkos_Core.hpp>

#include <gts.h>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <KM/Utils/KTypes.h>

#include <gmds/cadfac/FACSurface.h>
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/
/** \class FacetedGeomServices
 *  \brief This is a dummy class.
 */
/*----------------------------------------------------------------------------*/
    class FacetedSurfaceGeomServices
    {
    public:
        /*------------------------------------------------------------------------*/
        /** \brief  Constructor
         */
        FacetedSurfaceGeomServices();

        /*------------------------------------------------------------------------*/
        /** \brief Copy constructor
         */
        FacetedSurfaceGeomServices(const FacetedSurfaceGeomServices&) =delete;

        /*------------------------------------------------------------------------*/
        /** \brief  Destructor
         */
        ~FacetedSurfaceGeomServices() =default;

        /*------------------------------------------------------------------------*/
        /** \brief  Overloaded operator=
         */
        FacetedSurfaceGeomServices& operator=(const FacetedSurfaceGeomServices&) =delete;


        /*------------------------------------------------------------------------*/
        /** \brief  Builds the gts axis-aligned bounding box trees, one per surface.
         *          Needs to be called before actually using the structure.
         *
         *  \return
         */
        void buildAABBSurfacesTriangulationTrees(
                std::vector<gmds::cad::FACSurface*>& ASurfaces);

        /*------------------------------------------------------------------------*/
        /** \brief  Projects unto the surface
         *
         *  \return
         */
        void project(const gmds::cad::GeomSurface* ASurf,
                     gmds::math::Point& APoint) const;

        void
        project(gmds::math::Point& APoint, GNode* ATree) const;

        void
        project(gmds::math::Point& APoint) const;


        // projection function unto a triangle
        // given to gts to be used internally
        static GtsPoint* FacetedGeomServices_triangle_project(
                GtsPoint* AP,
                gpointer bounded)
        {
            gmds::math::Point p(AP->x,AP->y,AP->z);
            gmds::math::Triangle* tri = (gmds::math::Triangle*) bounded;
            gmds::math::Point newP = tri->project(p);
            GtsPoint* gtsp = gts_point_new(gts_point_class (),
                                           newP.X(),newP.Y(),newP.Z());

            return gtsp;
        }


    private:

        // Axis-Aligned Bounding Box tree for the surfaces triangles
        GNode* m_aabbSurfacesTrianglesTree;

        // Axis-Aligned Bounding Box tree for the surfaces triangles;
        // one tree for each surface
        std::map<gmds::cad::GeomSurface* const, GNode*> m_aabbSurfacesTrianglesTrees;

        std::map<gmds::cad::GeomSurface* const, std::vector<gmds::math::Triangle> > m_surfacesTriangulation;

        //        std::map<gmds::cad::GeomSurface*,std::vector<gmds::math::Point> > m_surfacesBoundingBox;


    };
/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_FACETEDSURFACEGEOMSERVICES_H_ */
