/*----------------------------------------------------------------------------*/
/** \file    FacetedVolume.t.h
 *  \author  F. LEDOUX
 *  \date    29/06/2011
 */
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/cadfac/FACVolume.h>
#include <algorithm>
#include <set>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
    namespace cad{
/*----------------------------------------------------------------------------*/
        int FACVolume::m_next_id=1;
/*----------------------------------------------------------------------------*/
        void FACVolume::resetIdCounter() {m_next_id=1;}
/*----------------------------------------------------------------------------*/
        FACVolume::FACVolume(const std::string& AName)
                :GeomVolume(AName),m_id(m_next_id++)
        {}
/*----------------------------------------------------------------------------*/
        FACVolume::~FACVolume()
        {}
/*----------------------------------------------------------------------------*/
        TCoord FACVolume::computeArea() const {
            throw GMDSException("FACVolume::computeArea: not yet implemented");
        }
/*----------------------------------------------------------------------------*/
        void FACVolume::computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const
        {
            if(m_adjacent_surfaces.empty()){
                throw GMDSException("The bounding box of a volume can only be computed if it is enclosed by surfaces!");

            }
            //we have at least one surface
            m_adjacent_surfaces[0]->computeBoundingBox(minXYZ,maxXYZ);
            for(auto s_id=1; s_id<m_adjacent_surfaces.size(); s_id++){
                GeomSurface* s = m_adjacent_surfaces[s_id];
                TCoord  min_s[3], max_s[3];
                s->computeBoundingBox(min_s,max_s);
		          minXYZ[0] = std::min(minXYZ[0],min_s[0]);
		          minXYZ[1] = std::min(minXYZ[1],min_s[1]);
		          minXYZ[2] = std::min(minXYZ[2],min_s[2]);
		          maxXYZ[0] = std::max(maxXYZ[0],max_s[0]);
		          maxXYZ[1] = std::max(maxXYZ[1],max_s[1]);
		          maxXYZ[2] = std::max(maxXYZ[2],max_s[2]);
            }
        }
/*----------------------------------------------------------------------------*/
        std::vector<GeomPoint*>& FACVolume::points() {
            return m_adjacent_points;
        }
/*----------------------------------------------------------------------------*/
        std::vector<GeomCurve*>& FACVolume::curves() {
            return m_adjacent_curves;
        }
/*----------------------------------------------------------------------------*/
        std::vector<GeomSurface*>& FACVolume::surfaces() {
            return m_adjacent_surfaces;
        }
/*----------------------------------------------------------------------------*/
    } // namespace cad
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
