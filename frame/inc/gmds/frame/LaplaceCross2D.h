/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux (2015)
 *
 * franck.ledoux@cea.fr
 *
 * The FRAME software is a computer program whose purpose is to provide a set
 * of algorithms to build 2D and 3D meshes using frame field concept. The main
 * focus of these algorithms is quadrilateral and hexahedral meshing.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability. 
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or 
 * data to be ensured and,  more generally, to use and operate it in the 
 * same conditions as regards security. 
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
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
