/*----------------------------------------------------------------------------*/
//
// Created by totoro on 2019-04-13.
//
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/utils/Exception.h>
#include <gmds/utils/CommonTypes.h>
#include <gmds/cad/GeomMeshLinker.h>
#include <gmds/cad/FACManager.h>
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOMSMOOTHER_H
#define GMDS_GEOMSMOOTHER_H
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
    namespace cad{
/*----------------------------------------------------------------------------*/
/** \class GeomPoint
 *  \brief This class describe the services that are required by the
 *  	   mesh to the geometrical model. As a consequence, this interface only
 *  	   contains query methods.
 */
/*----------------------------------------------------------------------------*/
        class EXPORT_GMDS GeomSmoother {
        public:
            GeomSmoother(GeomMeshLinker* ALinker);
            void smoothCurves(const int ANbIterations=1);
            void smoothSurfaces(const int ANbIterations=1);
            void smoothVolumes(const int ANbIterations=1);

            bool isValid() const;
        private:
            void init();
        private:
            GeomMeshLinker* m_linker;
            /** for each curve id, we store the mesh node ids*/
            std::map<int, std::vector<TCellID> > m_c2n;
            /** for each surface id, we store the mesh node ids*/
            std::map<int, std::vector<TCellID> > m_s2n;
            /** for each volume id, we store the mesh node ids*/
            std::map<int, std::vector<TCellID> > m_v2n;
            /** local n2n connection build for the purposes of smoothing*/
            std::map<TCellID, std::vector<TCellID> > m_n2n;
        };
/*----------------------------------------------------------------------------*/
    } // end namespace cad
/*----------------------------------------------------------------------------*/
} // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif //GMDS_GEOMSMOOTHER_H
/*----------------------------------------------------------------------------*/
