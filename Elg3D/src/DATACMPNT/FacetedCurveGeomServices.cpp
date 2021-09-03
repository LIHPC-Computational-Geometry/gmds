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
/** \file    FacetedCurveGeomServices.cpp
 *  \author  legoff
 *  \date    04/09/2019
 */
/*----------------------------------------------------------------------------*/
#include "ELG3D/DATACMPNT/FacetedCurveGeomServices.h"
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// Elg3D File Headers
/*----------------------------------------------------------------------------*/
#include <KM/Utils/Exception.h>
/*----------------------------------------------------------------------------*/
namespace elg3d {
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
    FacetedCurveGeomServices::FacetedCurveGeomServices()
    {

    }

/*----------------------------------------------------------------------------*/
    void
    FacetedCurveGeomServices::buildAABBCurvesTriangulationTrees(std::vector<gmds::cad::FACCurve*>& ACurves)
    {
        for(size_t iCurve=0; iCurve<ACurves.size(); iCurve++)
        {
            std::vector<gmds::math::Segment> segments;
            ACurves[iCurve]->getTriangulation(segments);
            m_curvesTriangulation[ACurves[iCurve]] = segments;
        }


        GSList* boxList = NULL;

        for(unsigned int iCurve=0; iCurve<ACurves.size(); iCurve++) {

            std::vector<gmds::math::Segment>& curveTriangulation = m_curvesTriangulation[ACurves[iCurve]];

            GSList* boxList_local = NULL;

            for(unsigned int iSegment=0; iSegment<curveTriangulation.size(); iSegment++) {
                double minXYZ[3];
                double maxXYZ[3];

                curveTriangulation[iSegment].computeBoundingBox(minXYZ,maxXYZ);

                gpointer pointer = &(m_curvesTriangulation[ACurves[iCurve]][iSegment]);
                GtsBBox* bbox = gts_bbox_new(
                        gts_bbox_class (),
                        pointer,
                        minXYZ[0],minXYZ[1],minXYZ[2],
                        maxXYZ[0],maxXYZ[1],maxXYZ[2]);

                boxList = g_slist_prepend(boxList,bbox);
                boxList_local = g_slist_prepend(boxList_local,bbox);
            }
            GNode* boxTree_local = gts_bb_tree_new(boxList_local);
            if(boxTree_local == NULL) {
                throw kmds::KException("FacetedMeshIntersectionService::buildAABBCurvesTriangulationTree : failed to build local tree");
            }
            m_aabbCurvesSegmentsTrees[ACurves[iCurve]] = boxTree_local;
        }


        std::cout<<"boxList "<<g_slist_length(boxList)<<std::endl;
        GNode* boxTree = gts_bb_tree_new(boxList);

        if(boxTree == NULL) {
            throw kmds::KException("FacetedMeshIntersectionService::buildAABBCurvesTriangulationTree : failed to build tree");
        }
        m_aabbCurvesSegmentsTree = boxTree;
    }

    /*----------------------------------------------------------------------------*/
    void
    FacetedCurveGeomServices::project(gmds::math::Point& APoint, GNode* ATree) const
    {
        GtsPoint* p = gts_point_new(gts_point_class (),
                                    APoint.X(),APoint.Y(),APoint.Z());

        gdouble* distance = NULL;
        GtsPoint* newP = gts_bb_tree_point_closest(
                ATree,
                p,
                FacetedGeomServices_segment_project,
                distance);

        APoint.setXYZ(newP->x,newP->y,newP->z);
    }
/*----------------------------------------------------------------------------*/
    void
    FacetedCurveGeomServices::project(
            const gmds::cad::GeomCurve* ACurv,
            gmds::math::Point& APoint) const
    {

        if((this->m_aabbCurvesSegmentsTrees).find(const_cast<gmds::cad::GeomCurve*> (ACurv)) == (this->m_aabbCurvesSegmentsTrees).end()) {
            throw kmds::KException("FacetedCurveGeomServices::project : curve not found in aabbCurvesSegmentsTrees_. Probably missing init.");
        }

        this->project(APoint,(this->m_aabbCurvesSegmentsTrees).find(const_cast<gmds::cad::GeomCurve*> (ACurv))->second);

//        this->project(APoint, m_aabbCurvesSegmentsTree);
    }

/*----------------------------------------------------------------------------*/
}  // end namespace elg3d
/*----------------------------------------------------------------------------*/
