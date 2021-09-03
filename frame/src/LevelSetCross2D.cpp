/*----------------------------------------------------------------------------*/
#include <iostream>
/*---------------------------------------------------------------------------*/
#include <gmds/igalgo/BoundaryOperator.h>
#include <GMDS/Algo/DistanceFieldBuilder2D.h>
#include "gmds/math/Numerics.h"
/*---------------------------------------------------------------------------*/
#include <gmds/io/VTKWriter.h>
#include <sstream>
/*---------------------------------------------------------------------------*/
#include "LevelSetCross2D.h"
/*---------------------------------------------------------------------------*/
using namespace gmds;
/*---------------------------------------------------------------------------*/
LevelSetCross2D::LevelSetCross2D(gmds::IGMesh* AMesh,
				 gmds::Variable<gmds::math::Cross2D>* AField,
				 std::vector<gmds::Node>& AFromNodes,
				 std::vector<gmds::Node>& AToNodes,
				 const int AMarkFace)
  : m_mesh(AMesh), m_field(AField), m_from_nodes(AFromNodes),
    m_to_nodes(AToNodes), m_mark_face(AMarkFace)
{ 
  m_mark_alive = m_mesh->getNewMark<Node>();
  
  for(unsigned int i=0; i <m_from_nodes.size(); i++)
    m_mesh->mark(m_from_nodes[i],m_mark_alive);
  
  m_distance_field = 0;
}
/*---------------------------------------------------------------------------*/
void LevelSetCross2D::execute()
{
  //==================================================================
  //STEP 1- Computes the distance field(s)
  //==================================================================
  std::cout << "======================================"<< std::endl;
  std::cout << "Computation of distance fields" << std::endl;
  computeDistanceFields();
  std::cout << "    DONE" << std::endl;
  


  //==================================================================
  // Now, we assign a cross to every  nodes !!!
  // We follow the order of insertion done to compute the distance 
  // field.
  //==================================================================
  // Step 2 - Cross field propagation
  //==================================================================
  // A node of the narrow band must be adjacent to 2 alive nodes, 
  // while being not alive.
  std::cout << "======================================"<< std::endl;
  std::cout << "Adv.-front cross field generation" << std::flush;
  propagateCrosses();
  std::cout << "    DONE" << std::endl;

  //==================================================================
  // CLEANING - Boolean marks are cleaned
  //==================================================================
  m_mesh->unmarkAll<Node>(m_mark_alive);
  m_mesh->freeMark<Node>(m_mark_alive);
}
     
/*---------------------------------------------------------------------------*/
void LevelSetCross2D::computeDistanceFields()
{
 
  gmds::DistanceFieldBuilder2D distanceBuilder(m_mesh);
  std::cout << "Computation of a 2D distance field" << std::endl;
  if (!distanceBuilder.isValid())
    {
      std::cout << "Invalid model for distance computation" << std::endl;
      throw GMDSException("Invalid model for distance computation");
    }
   m_distance_field = distanceBuilder.computeDistance(m_from_nodes, 
						     m_to_nodes,
						     m_mark_face);
  
  // We store the node ids in the order they have been inserted to 
  // compute the distance field
  m_insertion_order = distanceBuilder.getInsertionOrder();
  std::cout<<"ordered list: "<<m_insertion_order.size()<<std::endl;
  std::cout<<"to list     : "<<m_to_nodes.size()<<std::endl;
  std::cout<<"all         : "<<m_mesh->getNbNodes()<<std::endl;

}
/*---------------------------------------------------------------------------*/
void LevelSetCross2D::propagateCrosses()
{
  for(unsigned int i=0;i<m_to_nodes.size(); i++) {
    Node n_i = m_to_nodes[i];
    (*m_field)[n_i.getID()] =  math::Cross2D(0);
  }

  for(unsigned int i=0;i<m_insertion_order.size(); i++) {
    Node trial_node = m_mesh->get<Node>(m_insertion_order[i]);
    m_mesh->mark(trial_node, m_mark_alive);
    
    //now among the adjacent nodes, we get the ones thar are alive
    std::vector<Node> trial_adj_nodes = getAdjacentNodes(trial_node);
    std::vector<Node> alive_adj_nodes;

    for (unsigned int j = 0; j < trial_adj_nodes.size(); j++) {
      Node n_j = trial_adj_nodes[j];
      if (m_mesh->isMarked(n_j, m_mark_alive))
	alive_adj_nodes.push_back(n_j);
    }
    
    // We can now extrapolate a quaternion in this node... 
    extrapolateCross(trial_node, alive_adj_nodes);
  }
}
/*----------------------------------------------------------------------------*/
std::vector<Node> LevelSetCross2D::getAdjacentNodes(Node& ANode) const
{
  std::set<Node> adj_nodes;
  std::vector<Face>  adj_faces = ANode.get<Face>();

  for (unsigned int j = 0; j < adj_faces.size(); j++)
    {
      Face f_j = adj_faces[j];
      std::vector<Node> nodes_fj = f_j.get<Node>();
      for (unsigned int j = 0; j < nodes_fj.size(); j++)
	{
	  Node n_j = nodes_fj[j];
	  if (n_j.getID() != ANode.getID())
	    {
	      adj_nodes.insert(n_j);
	    }
	}
    }
  std::vector<Node> adj;
  adj.insert(adj.end(), adj_nodes.begin(), adj_nodes.end());

  return adj;
}
	 

/*----------------------------------------------------------------------------*/
void LevelSetCross2D::extrapolateCross(gmds::Node& ATo, 
				       std::vector<gmds::Node>& AFrom)
{
  math::Point p0 = ATo.getPoint();

  //===================================================================
  //We get crosses and  coeff (inverse of distance)
  //===================================================================
  std::vector<math::Cross2D> from_crosses;
  std::vector<TCoord> from_weights;
  for (unsigned int i = 0; i < AFrom.size(); i++) {
    Node current = AFrom[i];
    math::Cross2D current_c = (*m_field)[current.getID()];
    from_crosses.push_back(current_c);
    //coeef
    math::Point p1 = current.getPoint();
    math::Vector v01(p0, p1);
    from_weights.push_back(1.0/v01.norm());
  }
  //===================================================================
  //2D cross extrapolated
  //===================================================================
  
  math::Cross2D c_mean = math::Cross2D::mean(from_crosses,from_weights);
  
  //===================================================================
  //We assign the 2D cross on the mesh node
  //===================================================================
  (*m_field)[ATo.getID()] = c_mean;
 
}
/*----------------------------------------------------------------------------*/
