/*----------------------------------------------------------------------------*/
#ifndef LAPLACE_CROSS_2D_H_
#define LAPLACE_CROSS_2D_H_
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <map>
#include <set>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/math/Cross2D.h>
#include "LIB_GMDS_FRAME_export.h"
//#include <GMDS/Algo/DistanceFieldBuilder2D.h> //?doesnt exist anymore, adapt it!
/*---------------------------------------------------------------------------*/
// Usage of the Eigen Template library
#include <Eigen/Sparse>
/*----------------------------------------------------------------------------*/
/** \class  LaplaceCross2D
 *  \brief  Compute a 2D cross field on a 2D triangular mesh by solving a
 *          simple Laplace equation via FEM
 */
class LIB_GMDS_FRAME_API LaplaceCross2D
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
  LaplaceCross2D(gmds::Mesh* AMesh,
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
  /*------------------------------------------------------------------------*/
  /** \brief fill in the matrix A in Ax=b with the stiffness contribution
   */
  void buildStiffnessMatrix(Eigen::SparseMatrix<double>& AMat);

  /*------------------------------------------------------------------------*/
  /** \brief fill in the right-hand term b in Ax=b with Dirichlet conditions
   */
  void buildRightHandDirichlet(Eigen::VectorXd& ABCos, Eigen::VectorXd& ABSin);

  /*------------------------------------------------------------------------*/
  /** \brief Return the maximum number of adjacency for a node we work on
   */
  int getMaxAdjacency() const;

  /*---------------------------------------------------------------------------*/
  std::vector<gmds::Edge>  getEdgesOnCurve(const gmds::Node& ANode) const;
  gmds::Node getNeighboorOn(const gmds::Node& ANode, const gmds::Edge& AEdge) const;
 private:
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

  /*** Set of faces we work on */
  std::vector<gmds::Face> m_faces;

  /*** We store the correspondance from node ids to local index in m_nodes*/
  std::map<gmds::TCellID, int> m_id_to_local_index;
  std::map<gmds::TCellID, int> m_id_to_I;

  /*** Mark used to find the nodes of m_from_nodes */
  int  m_mark_dirichlet;
};
/*----------------------------------------------------------------------------*/
#endif //LAPLACE_CROSS_2D_H_
/*----------------------------------------------------------------------------*/
