/*----------------------------------------------------------------------------*/
#ifndef SH_EIGEN_FIELD_SOLVER_STRATEGY_H_
#define SH_EIGEN_FIELD_SOLVER_STRATEGY_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
/*---------------------------------------------------------------------------*/
// Frame File Headers
#include <gmds/frame3d/FieldSolverStrategyItf.h>
#include "LIB_GMDS_FRAME_3D_export.h"
/*----------------------------------------------------------------------------*/
namespace  gmds {
/*----------------------------------------------------------------------------*/
/** \class  EigenFieldSolverStrategy
 *  \brief  Implementation of the FieldSolverStrategyItf interface used by an
 *          FieldGenerator object to build a frame field.
 */
class LIB_GMDS_FRAME_3D_API EigenFieldSolverStrategy: public FieldSolverStrategyItf {
public:
    
    /*------------------------------------------------------------------------*/
    /** \brief Default constructor
     */
    EigenFieldSolverStrategy();
    
    /*------------------------------------------------------------------------*/
    /** \brief Function to be called for solving the system
     */
    virtual void solve();
    virtual void addSmoothingTerms();
    virtual void addBoundaryConstraints();
    virtual void addLockedTerms();
    virtual void initializeAssembly();
    virtual void finalizeAssembly();
    virtual void clean();
    virtual void setX();
    virtual void addLocalOptimConstraints();
    virtual void getFeasibleSolution();

private:
    //    int m_equation_numbering;
    //
    //    /*** Matrix A of the ||AX-b||^2 */
    //    Eigen::SparseMatrix<double> m_A;
    //    /*** Vector X of the ||AX-b||^2 */
    //    Eigen::VectorXd m_X;
    //    /*** Vector b of the ||AX-b||^2 */
    //    Eigen::VectorXd m_b;
    

    
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* SH_EIGEN_FIELD_SOLVER_STRATEGY_H_ */
/*----------------------------------------------------------------------------*/
