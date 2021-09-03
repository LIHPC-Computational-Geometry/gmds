/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    FacetedCurveGeomServices.h
 *  \author  legoff
 *  \date    07/10/2019
 */
/*----------------------------------------------------------------------------*/
#ifndef ELG3D_FACETEDCURVEGEOMSERVICES_H_
#define ELG3D_FACETEDCURVEGEOMSERVICES_H_
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

#include <gmds/cad/FACCurve.h>
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
    class FacetedCurveGeomServices
    {
    public:
        /*------------------------------------------------------------------------*/
        /** \brief  Constructor
         */
        FacetedCurveGeomServices();

        /*------------------------------------------------------------------------*/
        /** \brief Copy constructor
         */
        FacetedCurveGeomServices(const FacetedCurveGeomServices&) =delete;

        /*------------------------------------------------------------------------*/
        /** \brief  Destructor
         */
        ~FacetedCurveGeomServices() =default;

        /*------------------------------------------------------------------------*/
        /** \brief  Overloaded operator=
         */
        FacetedCurveGeomServices& operator=(const FacetedCurveGeomServices&) =delete;


        /*------------------------------------------------------------------------*/
        /** \brief  Builds the gts axis-aligned bounding box trees, one per curve.
         *          Needs to be called before actually using the structure.
         *
         *  \return
         */
        void buildAABBCurvesTriangulationTrees(
                std::vector<gmds::cad::FACCurve*>& ACurves);

        /*------------------------------------------------------------------------*/
        /** \brief  Projects unto the curve
         *
         *  \return
         */
        void project(const gmds::cad::GeomCurve* ACurv,
                     gmds::math::Point& APoint) const;

        void
        project(gmds::math::Point& APoint, GNode* ATree) const;


        // projection function unto a segment
        // given to gts to be used internally
        static GtsPoint* FacetedGeomServices_segment_project(
                GtsPoint* AP,
                gpointer bounded)
        {
            gmds::math::Point p(AP->x,AP->y,AP->z);
            gmds::math::Segment* seg = (gmds::math::Segment*) bounded;
            gmds::math::Point newP = seg->project(p);
            GtsPoint* gtsp = gts_point_new(gts_point_class (),
                                           newP.X(),newP.Y(),newP.Z());

            return gtsp;
        }


    private:

        // Axis-Aligned Bounding Box tree for the curves segments
        GNode* m_aabbCurvesSegmentsTree;

        // Axis-Aligned Bounding Box tree for the curves segments;
        // one tree for each curve
        std::map<gmds::cad::GeomCurve* const, GNode*> m_aabbCurvesSegmentsTrees;

        std::map<gmds::cad::GeomCurve* const, std::vector<gmds::math::Segment> > m_curvesTriangulation;

        //        std::map<gmds::cad::GeomCurve*,std::vector<gmds::math::Point> > m_curvesBoundingBox;


    };
/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
#endif /* ELG3D_FACETEDCURVEGEOMSERVICES_H_ */
