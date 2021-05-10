/*----------------------------------------------------------------------------*/
/** \file    FacetedVolume.t.h
 *  \author  F. LEDOUX
 *  \date    29/06/2011
 */
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/cad/FACVolume.h>
#include <set>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
    namespace cad{
/*----------------------------------------------------------------------------*/
        int FACVolume::m_next_id=1;
/*----------------------------------------------------------------------------*/
        FACVolume::FACVolume(const std::string& AName)
                :GeomVolume(AName)
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
                for(auto i=0;i<3;i++){
                    if(min_s[i]<minXYZ[i]){
                        minXYZ[i]=min_s[i];
                    }
                    if(max_s[i]>maxXYZ[i]){
                        maxXYZ[i]=max_s[i];
                    }

                }
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
