/*----------------------------------------------------------------------------*/
#ifndef GMDS_SURFACE_REORIENT_H
#define GMDS_SURFACE_REORIENT_H
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include "GMDSIgAlgo_export.h"
/*----------------------------------------------------------------------------*/
#include <map>
#include <vector>
/*----------------------------------------------------------------------------*/
namespace gmds{

    /*----------------------------------------------------------------------------*/
    /** @class  SurfaceReorient
     *  @brief  Class that provides an algorithm to reorient surfaces in a consistent
     *  manner
     */
    class GMDSIgAlgo_API SurfaceReorient
    {
    public:

        /*------------------------------------------------------------------------*/
        /** @brief Constructor.
         *
         * @param Amesh the mesh to work on
         * @param ADim the grid dimension: 2 or 3 (default). Dimension of 2 means
         * that we know we work in the 2D plane. In this case the orientation is
         * consistent with FEM requirements.
         */
	   SurfaceReorient(Mesh* AMesh, const TInt ADim=3);

        /*------------------------------------------------------------------------*/
        /** @brief  Destructor.	*/
        virtual ~SurfaceReorient();

        /*------------------------------------------------------------------------*/
        /** @brief Check if the mesh fits algorithm requirements, which are:
         *         - 3D mesh with R, N and R2N
         *         - 2D mesh with F, N and F2N
         */
        bool isValid() const;
        /*------------------------------------------------------------------------*/
        /** @brief  Performs the reorientation algorithm.
         * @return the number of faces that have been reoriented
         */
        int execute();

	   private:
	     int orient2d();
	     int orient3d();
	     bool orient2d(gmds::Face& AF);
	     static TCoord isLeft(Node& AN1, Node& AN2, Node& AN3);
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
