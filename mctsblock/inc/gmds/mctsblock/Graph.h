/*----------------------------------------------------------------------------*/
#ifndef GMDS_MCTS_GRAPH_H
#define GMDS_MCTS_GRAPH_H
/*----------------------------------------------------------------------------*/
#include <LIB_GMDS_MCTSBLOCK_export.h>
#include <gmds/utils/CommonTypes.h>
#include <gmds/utils/Exception.h>
/*----------------------------------------------------------------------------*/
#include <limits.h>
#include <set>
#include <map>
#include <memory>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace mctsblock {
/*----------------------------------------------------------------------------*/
// A structure to represent a
// node in adjacency list
struct LIB_GMDS_MCTSBLOCK_API AdjListNode
{
  TCellID dest;
  double weight;
  std::shared_ptr<AdjListNode> next;
  AdjListNode(TCellID ANodeID=0, double AWeight=0, std::shared_ptr<AdjListNode> AN=nullptr);
};
/*----------------------------------------------------------------------------*/
// A structure to represent
// an adjacency list
struct LIB_GMDS_MCTSBLOCK_API AdjList
{
  // Pointer to head node of list
  std::shared_ptr<AdjListNode> head;
  AdjList( std::shared_ptr<AdjListNode> h=nullptr);
};
/*----------------------------------------------------------------------------*/
// A structure to represent a graph.
// A graph is an array of adjacency lists.
// Size of array will be V (number of
// vertices in graph)
class LIB_GMDS_MCTSBLOCK_API Graph
{
public:
  Graph(const std::set<TCellID>& ANodeIDs);
  //Be careful, input ids are in the global input space numbering, not the local one
  void addEdge(TCellID ASrcNode, TCellID ADestNode, double AWeight);

  // The main function that calculates
  // distances of shortest paths from src to all
  // vertices. It is a O(ELogV) function
  void computeDijkstra(TCellID ASrcNode);
  //Be careful, input ids are in the global input space numbering, not the local one
 	void setWeight(const TCellID AN1, const TCellID AN2, const double AW);
	/**
	 *
	 * @return for each node of the graph D, you get the path from the source node S to D.
	 * The first item of the vector is D, and the last one S.
	 */
  std::map<TCellID , std::vector<TCellID> > getShortestPath();
  std::map<TCellID , double > getShortestPathWeights();
private:
  std::vector<std::pair<TCellID , double> > getAdjNodes(const TCellID  ANodeId);
  void buildShortestPaths(const TCellID ASrc);
  // A utility function to create
  // a new adjacency list node
  void  newAdjListNode(TCellID ADest, double AWeight, std::shared_ptr<AdjList> ANode);
  void  setWeightOneWay(const TCellID ASrc, const TCellID ADest, const double AWeight);

private:
  int m_nb_vertices;
  std::set<TCellID> m_input_vertices;
  std::vector<TCellID> m_graph_to_input_vertices;
  std::map<TCellID , TCellID> m_input_to_graph_vertices;
  std::vector<std::shared_ptr<AdjList>> m_array;

  std::vector<double> m_dist;

  std::map<TCellID , std::vector<TCellID> > m_shortest_path_nodes;
  std::map<TCellID , double > m_shortest_path_weight;
};
/*----------------------------------------------------------------------------*/
}     // namespace mctsblock
/*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/
#endif     // GMDS_MCTS_GRAPH_H


