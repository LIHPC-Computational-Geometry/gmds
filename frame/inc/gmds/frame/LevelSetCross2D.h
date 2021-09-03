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
#ifndef LEVEL_SET_CROSS_2D_H_
#define LEVEL_SET_CROSS_2D_H_
/*----------------------------------------------------------------------------*/
#include <map>
#include <list>
#include <set>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/math/Cross2D.h>
/*----------------------------------------------------------------------------*/
/** \class  LevelSetCross2D
 *  \brief  Propagate a cross field from crosses defined on the boundary 
 *          into the mesh
 */
class LIB_GMDS_FRAME_API LevelSetCross2D
{
 public:

  /*------------------------------------------------------------------------*/
  /** \brief Constructor.
   *
   *  \param AMesh     the mesh where we work on
   *  \param AField    the cross field to propagate
   *  \param AFrom     the nodes we start from
   *  \param ATo       the nodes we propagate through
   *  \param AFaceMark a mark indicating the face we can work on
   */
  LevelSetCross2D(gmds::Mesh* AMesh,
		  gmds::Variable<gmds::math::Cross2D>* AField,
		  std::vector<gmds::Node>& AFromNodes,
		  std::vector<gmds::Node>& AToNodes,
		  const int AFaceMark);
  /*------------------------------------------------------------------------*/
  /** \brief Function to be called for executing the cross field generation 
   *	     algorithm. 
   */
  void execute();
     
 private:  


  /*------------------------------------------------------------------------*/
  /** \brief Compute the distance field that will be used to define the
   *         insertion order of the advancing front process
   */
  void computeDistanceFields();

  /*------------------------------------------------------------------------*/
  /** \brief Return all the nodes sharing a face with ANode
   *  
   * \param ANode the node we want the neighborhood
   *
   * \return the nodes sharing a face with ANode 
   */
   std::vector<gmds::Node> getAdjacentNodes(gmds::Node& ANode) const;


   /*------------------------------------------------------------------------*/
   /** \brief Extrapolate a cross value in ATo from the crosses already 
    *         defined in AFrom2
    *  
    * \param ATo    a node
    * \param AFrom  a vector of nodes
    */
   void extrapolateCross(gmds::Node& ATo, std::vector<gmds::Node>& AFrom);

   /*------------------------------------------------------------------------*/
   /** \brief Iterative process to assign a cross to every node of m_mesh
    */
   void propagateCrosses();


 private:

   /*** Background triangular mesh */
   gmds::Mesh* m_mesh; 
  
   /*** Cross field that we build on the mesh (node association) */
   gmds::Variable<gmds::math::Cross2D>* m_field;

   /*** Set of nodes we start from */
   std::vector<gmds::Node> m_from_nodes;

   /*** Set of nodes we must defined a cross on */
   std::vector<gmds::Node> m_to_nodes;

   /*** Boolean mark for faces we can work on */
   int m_mark_face;

   /*** Distance field that is going to be computed on m_mesh */
   gmds::Variable<double>* m_distance_field; 
	  
   /*** Distance field that is going to be computed on m_mesh */
   std::vector<gmds::TCellID> m_insertion_order;

   /*** Boolean mark for nodes where a cross is already defined */
   int m_mark_alive;

};
/*----------------------------------------------------------------------------*/
#endif //ADVANCING_FRONT_FIELD_GEN_2D_H_
/*----------------------------------------------------------------------------*/
