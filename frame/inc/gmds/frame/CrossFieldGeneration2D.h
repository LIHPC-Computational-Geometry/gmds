/*----------------------------------------------------------------------------*/
#ifndef CROSS_FIELD_GENERATION_2D_H_
#define CROSS_FIELD_GENERATION_2D_H_
/*----------------------------------------------------------------------------*/
#include <vector>
#include <set>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/math/Cross2D.h>
/*----------------------------------------------------------------------------*/
/** \class  CrossFieldGeneration2D
 *  \brief  Computes a 2D cross field onto a 2D triangular mesh by solving a
 *          simple Laplace equation via FEM
 */
class EXPORT_GMDS CrossFieldGeneration2D {

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
