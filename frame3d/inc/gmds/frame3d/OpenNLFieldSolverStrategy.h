/*----------------------------------------------------------------------------*/
#ifndef SH_OPENNL_FIELD_SOLVER_STRATEGY_H_
#define SH_OPENNL_FIELD_SOLVER_STRATEGY_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
/*---------------------------------------------------------------------------*/
// Frame File Headers
#include <gmds/frame3d/FieldSolverStrategyItf.h>
#include "LIB_GMDS_FRAME_3D_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*----------------------------------------------------------------------------*/
    /** \class  OpenNLFieldSolverStrategy
     *  \brief  Implementation of the FieldSolverStrategyItf interface used by an
     *          FieldGenerator object to build a frame field.
     */
    class LIB_GMDS_FRAME_3D_API OpenNLFieldSolverStrategy: public FieldSolverStrategyItf {
    public:
        
        /*------------------------------------------------------------------------*/
        /** \brief Default constructor
         * \param AWithCurvature add an extra control on the curvature
         *        preservation
         */
        OpenNLFieldSolverStrategy(const bool AWithCurvature = false);
        
        /*------------------------------------------------------------------------*/
        /** \brief Function to be called for solving the system
         */
        virtual void solve();
        virtual void addSmoothingTerms();
        virtual void addBoundaryConstraints();
        virtual void addLockedTerms();
        virtual void initializeAssembly();
        virtual void finalizeAssembly();
        virtual void addLocalOptimConstraints();
        virtual void clean();
        virtual void setX();
        virtual void getFeasibleSolution();

    private:
        bool m_with_curvature;
                
    };
    /*------------------------------------------------------------------------*/
}//namespace fhedo
/*----------------------------------------------------------------------------*/
#endif /* SH_OPENNL_FIELD_SOLVER_STRATEGY_H_ */
/*----------------------------------------------------------------------------*/
