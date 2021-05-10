/*----------------------------------------------------------------------------*/
#ifndef FUNCTIONAL_APPROACH_CROSS_2D_H_
#define FUNCTIONAL_APPROACH_CROSS_2D_H_
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <map>
#include <set>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/math/Cross2D.h>
/*---------------------------------------------------------------------------*/
// Usage of the Eigen Template library
#include <Eigen/Sparse>
/*---------------------------------------------------------------------------*/
class CrossFieldGeneration2D;
/*----------------------------------------------------------------------------*/
/** \class  LaplaceCross2D
 *  \brief  Compute a 2D cross field on a 2D triangular mesh by solving a
 *          simple Laplace equation via FEM but with a functional approach
 *          where frames are represented into a Fourier function basis.
 */
class EXPORT_GMDS FunctionalApproachCross2D
{
public:
    
    /*------------------------------------------------------------------------*/
    /** \brief Constructor.
     *
     *  \param AMesh  the mesh where we work on
     *  \param AField the cross field to propagate, it is already
     *                initialized on AFrom
     *  \param AFrom  the nodes we start from
     *  \param ATo    the nodes we propagate through
     *  \param AFaces the faces we work on
     *
     *   Warning: AFrom and ATo must be disjoint!!!
     */
    FunctionalApproachCross2D(CrossFieldGeneration2D* AParent,
                              gmds::Mesh* AMesh,
                              gmds::Variable<gmds::math::Cross2D>* AField,
                              std::vector<gmds::Node>& AFromNodes,
                              std::vector<gmds::Node>& AToNodes,
                              std::vector<gmds::Face>& AFaces );
    
    /*------------------------------------------------------------------------*/
    /** \brief Function to be called for executing a frame field generation
     *	     algorithm. This function follows the template method DP. It
     *	     means that it drives the global execution algorithm, each
     *	     specific part being delegated to child classes.
     */
    void execute();
    
private:
    void initialization(Eigen::VectorXd& AX, const int AMaxNbIterations=1);
    void smoothing(Eigen::VectorXd& AX, const int AMaxNbIterations=1000,
                   const int AMaxInnerLoopIterations=10);
    double getFieldCurvature();
private:
    /*** calling class for debug purpose */
    CrossFieldGeneration2D* m_parent;
    
    /*** Background triangular mesh */
    gmds::Mesh* m_mesh;
    
    /*** Cross field that we build on the mesh (node association) */
    gmds::Variable<gmds::math::Cross2D>* m_field;
    
    /*** Set of nodes we start from */
    std::vector<gmds::Node> m_from_nodes;
    /*** Set of nodes we start from */
    std::vector<gmds::Node> m_to_nodes;
    
    /*** Set of nodes we work on */
    std::vector<gmds::Node> m_nodes;
    
    /*** Set of edges we work on */
    std::vector<gmds::Edge> m_edges;
    
    /*** We store the correspondance from node ids to local index in m_nodes*/
    std::map<gmds::TCellID, int> m_id_to_local_index;
    std::map<gmds::TCellID, int> m_id_to_I;
    
    /*** Mark used to find the nodes of m_from_nodes */
    int  m_mark_dirichlet;
};
/*----------------------------------------------------------------------------*/
#endif //FUNCTIONAL_APPROACH_CROSS_2D_H_
/*----------------------------------------------------------------------------*/
