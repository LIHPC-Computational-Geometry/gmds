/*----------------------------------------------------------------------------*/
#ifndef GMDS_GRID_BUILDER_H
#define GMDS_GRID_BUILDER_H
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include "GMDSIgAlgo_export.h"
/*----------------------------------------------------------------------------*/
#include <map>
#include <vector>
/*----------------------------------------------------------------------------*/
namespace gmds{

    /*----------------------------------------------------------------------------*/
    /** @class  GridBuilder
     *  @brief  Class that provides a way to create a structured grid
     */
    class GMDSIgAlgo_API GridBuilder
    {
    public:

        /*------------------------------------------------------------------------*/
        /** @brief Constructor.
         *
         * @param Amesh the mesh to work on
         * @param ADim the grid dimension: 2 or 3 (default)
         */
        GridBuilder(Mesh* AMesh, const TInt ADim=3);

        /*------------------------------------------------------------------------*/
        /** @brief  Destructor.	*/
        virtual ~GridBuilder();

        /*------------------------------------------------------------------------*/
        /** @brief Check if the mesh fits algorithm requirements, which are:
         *         - 3D mesh with R, N and R2N
         *         - 2D mesh with F, N and F2N
         */
        bool isValid() const;
        /*------------------------------------------------------------------------*/
        /** @brief  Performs the grid building (two last parameters are only used
         *          for 3D meshes)
         * @param AXNb      NB cells in X direction
         * @param AXStep    Length of each step in in X direction
         * @param AYNb      NB cells in Y direction
         * @param AYStep    Length of each step in in Y direction
         * @param AZNb      NB cells in Z direction
         * @param AZStep    Length of each step in in Z direction
         */
        void execute(const TInt  AXNb, const TCoord AXStep,
                     const TInt  AYNb, const TCoord AYStep,
                     const TInt  AZNb=1, const TCoord AZStep=1.0);

    private:
        void build2D(const TInt  AXNb, const TCoord AXStep,
                     const TInt  AYNb, const TCoord AYStep
        );
        void build3D(const TInt  AXNb, const TCoord AXStep,
                     const TInt  AYNb, const TCoord AYStep,
                     const TInt  AZNb, const TCoord AZStep
        );

    private:
        /** a mesh */
        Mesh* m_mesh;

        /** Grid dimension*/
        TInt  m_dim;
    };
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/

#endif //GMDS_GRID_BUILDER_H
