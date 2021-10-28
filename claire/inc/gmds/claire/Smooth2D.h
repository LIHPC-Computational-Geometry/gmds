/*----------------------------------------------------------------------------*/
#ifndef GMDS_CLAIRE_SMOOTH2D_H
#define GMDS_CLAIRE_SMOOTH2D_H
/*----------------------------------------------------------------------------*/
#include "LIB_GMDS_CLAIRE_export.h"
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
#include <string>
#include <map>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
    class LIB_GMDS_CLAIRE_API Smooth2D {
    public:
        /*--------------------------------------------------------------------*/
        /** @enum  Status code for executing algorithms
         */
        typedef enum {
            FAIL,
            SUCCESS
        } STATUS;

        /*-------------------------------------------------------------------*/
        /** @brief Constructor.
         *  @param AMesh the mesh where we work on
         *  @param AVarBnd node variable (value 1 means constrained node)
         *  @param ANbIterations nb max iterations
         */
        Smooth2D(Mesh *AMesh,
                 const Variable<int> *AVarBnd,
                 const int ANbIterations = 100);

        /*-------------------------------------------------------------------*/
        /** @brief Set the max number of iterations
         *  @param[in] ANbIterations
         */
        void setNbIterations(const int ANbIterations);

        /*-------------------------------------------------------------------*/
        /** @brief Execute the algorithm
         */
        STATUS execute();
    private:
        /*-------------------------------------------------------------------*/
        /** @brief Goes through all the free nodes in the mesh and build
         *         corresponding stencils
         */
        void buildStencils();
	     //void PerturbationMaillage(const Variable<int>* var_bnd, const double dx, const double dy);
	     math::Point FindMidBranche(const math::Point A, const math::Point B, const math::Point C);
	     //bool IntersectionSegments2D(const math::Point A, const math::Point B, const math::Point C, const math::Point D);
    private:
        /** mesh we work on */
        Mesh *m_mesh;
        /** Variable to store which nodes are constrained*/
        const Variable<int> *m_node_constrained;
        /** nb max iterations */
        int m_nb_max_iterations;

        /** free nodes in the mesh */
        std::vector<TCellID> m_free_nodes;
        typedef struct {
            unsigned int val[3][3];
        } stencil;
        std::map<TCellID, stencil> m_stencil;
    };
/*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_CLAIRE_SMOOTH2D_H
/*----------------------------------------------------------------------------*/
