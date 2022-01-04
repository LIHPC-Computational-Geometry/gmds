/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include <gmds/smoothy/AbstractSmoother.h>
#include "GMDSSmoothy_export.h"
/*----------------------------------------------------------------------------*/
#ifndef GMDS_ANGLE_BASED_QUAD_SMOOTHER_H
#define GMDS_ANGLE_BASED_QUAD_SMOOTHER_H
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
    namespace smoothy{
/*----------------------------------------------------------------------------*/
/** \class LaplacianSmoother
 *  \brief This class provides a straightforward implementation of a Laplacian
 *         smoother for meshes classifiend on geometric models. Curves, surfaces
 *         and volumes are smoothed individually.
 */
/*----------------------------------------------------------------------------*/
        class GMDSSmoothy_API AngleBasedQuadSmoother: public AbstractSmoother{
        public:
            AngleBasedQuadSmoother(cad::GeomMeshLinker* ALinker);
            void smooth(const int ANbIterations=1);


            /**@brief validation of input
             * @return check the linker (cf AbstractSmoother class) and check
             *         that we have only quad faces in the mesh.
             */
            virtual bool isValid() const;

        private:
            void smoothCurves();
            void smoothSurfaces();
        };
/*----------------------------------------------------------------------------*/
    } // end namespace cad
/*----------------------------------------------------------------------------*/
} // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif //GMDS_ANGLE_BASED_QUAD_SMOOTHER_H
/*----------------------------------------------------------------------------*/
