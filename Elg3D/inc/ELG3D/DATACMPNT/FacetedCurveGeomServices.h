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
