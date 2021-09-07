/*----------------------------------------------------------------------------*/
#ifndef FEM_SOLVER_STRATEGY_H_
#define FEM_SOLVER_STRATEGY_H_
/*---------------------------------------------------------------------------*/
// Frame File Headers
#include "FieldSolverStrategyItf.h"
#include "LIB_GMDS_FRAME_3D_export.h"
/*---------------------------------------------------------------------------*/
// Eigen File Headers
#include <Eigen/Sparse>
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*----------------------------------------------------------------------------*/
    /** \class  FEMSolverStrategy
     *  \brief  Implementation of the FieldSolverStrategyItf interface used by an
     *          FieldGenerator object to build a frame field.
     */
    class LIB_GMDS_FRAME_3D_API FEMSolverStrategy: public FieldSolverStrategyItf {
    public:

        /*------------------------------------------------------------------------*/
        /** \brief Default constructor
         */
        FEMSolverStrategy();

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
        void buildStiffnessMatrix(Eigen::SparseMatrix<double>& AM,
        const int i);
void buildRightHandDirichlet(Eigen::VectorXd& AB, const int AI);
    private:
            int m_equation_numbering;

            /*** Matrix A of the ||AX-b||^2 */
            Eigen::SparseMatrix<double> m_A[9];
            /*** Vector X of the ||AX-b||^2 */
            Eigen::VectorXd m_X[9];
            /*** Vector b of the ||AX-b||^2 */
            Eigen::VectorXd m_b[9];

        int m_dirichlet_mark[9];

    };
    /*------------------------------------------------------------------------*/
}//namespace gmds
/*----------------------------------------------------------------------------*/
#endif //FEM_SOLVER_STRATEGY_H_
/*----------------------------------------------------------------------------*/
