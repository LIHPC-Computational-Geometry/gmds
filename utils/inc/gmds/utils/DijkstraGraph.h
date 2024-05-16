/*----------------------------------------------------------------------------*/
#ifndef GMDS_DIJKSTRAGRAPH_H
#define GMDS_DIJKSTRAGRAPH_H
/*----------------------------------------------------------------------------*/
#include<bits/stdc++.h>

#include <gmds/utils/CommonTypes.h>


# define INF 0x3f3f3f3f
/*----------------------------------------------------------------------------*/
namespace gmds{

class DijkstraGraph
{
	int V;    // No. of vertices

	int idV = 0;

	// In a weighted graph, we need to store vertex
	// and weight pair for every edge
	std::list<std::pair<int, int>> *adj;
	std::map<TCellID , int> v2n;
	std::map<int, TCellID> n2v;
	std::map<TCellID, std::pair<int, int>> graphEdge2blockingEdge;
	std::map<TCellID, double> weightEdge; //TODO maybe not used

 public:
	DijkstraGraph(std::set<int> AV);  // Constructor

	void setWeightEdge(TCellID AEdge, double AWeight);

	// function to add an edge to graph
	void addEdge(TCellID u, TCellID v, double w, TCellID idEdgeBlocking);

	// prints shortest path from s
	void shortestPath(int s);
};

class GraphOrigin{
	int V;    // No. of vertices

	// In a weighted graph, we need to store vertex
	// and weight pair for every edge
	std::list< std::pair<int, int> > *adj;

 public:
	GraphOrigin(int V);  // Constructor

	// function to add an edge to graph
	void addEdgeGraph(int u, int v, int w);

	// prints shortest path from s
	void shortestPathGraph(int s);
};
}

#endif     // GMDS_DIJKSTRAGRAPH_H
