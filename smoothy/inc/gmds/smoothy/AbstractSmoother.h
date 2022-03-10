/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/cad/GeomMeshLinker.h>
#include "GMDSSmoothy_export.h"
/*----------------------------------------------------------------------------*/
#ifndef GMDS_ABSTRACT_SMOOTHER_H
#define GMDS_ABSTRACT_SMOOTHER_H
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
    namespace smoothy{
/*----------------------------------------------------------------------------*/
/** \class AbstractSmoother
 *  \brief This class provides attributes and methods that are shared by
 *         smoothing algorithms as child classes. It is abstract and so cannot
 *         be instanciated.
 */
/*----------------------------------------------------------------------------*/
        class GMDSSmoothy_API AbstractSmoother {
        public:
            /**@brief validation of the linker
             * @return true if the linker is valid, i.e. the underlying mesh
             *         has the right connection. False otherwise
             */
            virtual bool isValid() const;
        protected:
            /**@brief constructor
             * @param ALinker the linker that has the knowledge and connection
             *                between geometry and mesh
             */
            AbstractSmoother(cad::GeomMeshLinker* ALinker);

            /**@brief initialize all the required technical attributes$
             */
            void init();

            /** linker*/
            cad::GeomMeshLinker* m_linker;
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
#endif //GMDS_ABSTRACT_SMOOTHER_H
/*----------------------------------------------------------------------------*/
