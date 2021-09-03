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
#ifndef CROSS_FIELD_GENERATION_2D_H_
#define CROSS_FIELD_GENERATION_2D_H_
/*----------------------------------------------------------------------------*/
#include <vector>
#include <set>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/math/Cross2D.h>
#include "LIB_GMDS_FRAME_export.h"
/*----------------------------------------------------------------------------*/
/** \class  CrossFieldGeneration2D
 *  \brief  Computes a 2D cross field onto a 2D triangular mesh by solving a
 *          simple Laplace equation via FEM
 */
class LIB_GMDS_FRAME_API CrossFieldGeneration2D {

  public:
  /*------------------------------------------------------------------------*/
  /** \brief Strategies to compute the cross field
   */
    typedef enum {
        laplace_solve,
        level_sets,
        functional_approach,
        stab_balls
    } Strategy;
  /*------------------------------------------------------------------------*/
  /** \brief Constructor.
   *
   *  \param AMesh the mesh where we work on
   */
  CrossFieldGeneration2D(gmds::Mesh* AMesh);
     
  /*------------------------------------------------------------------------*/
  /** \brief specify the prefix of output files for debug
   *
   *  \param AName
   */
  void setDebugPrefix(const std::string& AName);
     
  /*------------------------------------------------------------------------*/
  /** \brief Function to be called for executing a frame field generation 
   *		   algorithm. This function follows the template method DP. It
   *		   means that it drives the global execution algorithm, each 
   *		   specific part being delegated to child classes.
   */
  void execute(const Strategy AStrategy);

    void writeForDebug(const std::string AFileName="");

 private:
    
    /*------------------------------------------------------------------------*/
    /** \brief Build the frame field by solving a global Laplace equation on the
     *   domain
     */
    void buildFieldViaLaplaceEDP();
    
    
    /*------------------------------------------------------------------------*/
    /** \brief Build the frame field by solving a global Laplace equation on the
     *   domain using a functional approach
     */
    void buildFieldViaFunctionalApproach();
    
  /*------------------------------------------------------------------------*/
  /** \brief Build the frame field by applying a advancing front algorithm,
   *   which uses a distance field to drive the insertion process.
   */
  void buildFieldViaLevelSets();

  /*------------------------------------------------------------------------*/
  /** \brief Initialize all the marks used for the algorithm
   */
  void initMarks();

  /*------------------------------------------------------------------------*/
  /** \brief Clean all the marks used for the algorithm
   */
  void cleanMarks();

  /*------------------------------------------------------------------------*/
  /** \brief Mark all the boundary cells of the mesh
   */
  void markBoundaryCells();

  /*------------------------------------------------------------------------*/
  /** \brief Compute the crosses that are put on the node classified on curves
   */
   void initCrossesOnCurves();
  /*------------------------------------------------------------------------*/
  /** \brief Compute the crosses that are put on the node classified on points
   */
   void initCrossesOnPoints();
   
  /*------------------------------------------------------------------------*/
  /** \brief Get the edges classified on a curve and incident to ANode
   *
   * \param ANode the node we want to get incident edges 
   *
   * \return the edges classified on a curve and incident to ANode
   */
   std::vector<gmds::Edge> getEdgesOnCurve(const gmds::Node& ANode) const;

  /*------------------------------------------------------------------------*/
  /** \brief Get the node connected to ANonde by AEdge
   *
   * \param ANode a node
   * \param AEdge an edge
   *
   * \return the node of AEdge, which is opposite to ANode
   */
   gmds::Node getNeighboorOn(const gmds::Node& ANode, 
			     const gmds::Edge& AEdge) const;
                             
/*------------------------------------------------------------------------*/
  /** \brief Computes the deviation per face as the mean between the angle deviation of the refence vectors at the corners of the face
   *
   * \param  
   *
   */
   void computeReferenceVectorDeviationPerFace();

   /*------------------------------------------------------------------------*/
   void smooth(const int AMark);
   void smoothAll();
   void colorSimplices();

 private:
  /** Background triangular mesh */
  gmds::Mesh* m_mesh; 

  /** Cross field that we build on the mesh (node association) */
  gmds::Variable<gmds::math::Cross2D>* m_cross_field_2D; 


  std::string m_debug_output;

  /** mark for nodes on curves */
  int m_markNodeOnCurv;
  /** mark for nodes on points */
  int m_markNodeOnPnt;
  /** mark for edges on curves */
  int m_markEdgeOnCurv;
  /** mark for faces we work on */
  int m_markFace;
  /** mark for isolated nodes (connected to nothing) */
  int m_markIsolated;

  /** all the nodes on a curve (but not on a point) */
  std::vector<gmds::Node> m_curve_nodes;
  /** all the nodes on a surface (not curve, not point) */
  std::vector<gmds::Node> m_surf_nodes; 

};
/*----------------------------------------------------------------------------*/
#endif //CROSS_FIELD_GENERATION_2D_H_
/*----------------------------------------------------------------------------*/
