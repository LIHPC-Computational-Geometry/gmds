#include "gmds/utils/DijkstraGraph.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
using namespace std;

/*----------------------------------------------------------------------------*/
DijkstraGraph::DijkstraGraph(std::set<int> AV)
{
	this->V = AV.size();
	adj = new std::list<std::pair<int, int>>[V];
}

void DijkstraGraph::setWeightEdge(gmds::TCellID AEdge, double AWeight)
{
	auto idVertEdge = graphEdge2blockingEdge.find(AEdge)->second;
	//int idVertU = v2n.find(idVertEdge.first)->second;
	int idVertU = idVertEdge.first;
	int idVertV = idVertEdge.second;
	//int idVertV = v2n.find(idVertEdge.second)->second;
	list< pair<int, int> >::iterator i;
	for (i = adj[idVertU].begin(); i != adj[idVertU].end(); ++i){
		int v = (*i).first;
		if(v == idVertV){
			(*i).second = AWeight;
		}
	}
	list< pair<int, int> >::iterator j;
	for (j = adj[idVertV].begin(); j != adj[idVertV].end(); ++j){
		int v = (*j).first;
		if(v == idVertU){
			(*j).second = AWeight;
		}
	}

	//adj[idVertU].push_back(std::make_pair(idVertV, AWeight));
	//adj[idVertV].push_back(std::make_pair(idVertU, AWeight));

}

void
DijkstraGraph::addEdge(TCellID u, TCellID v, double w, TCellID idEdgeBlocking)
{
	int idVertU,idVertV;
	if(v2n.find(u)!=v2n.end()){
		if(v2n.find(v) != v2n.end()){
			idVertU = v2n.find(u)->second;
			idVertV = v2n.find(v)->second;
		}
		else{
			idVertU = v2n.find(u)->second;
			idVertV = idV;
			idV++;
		}
	}
	else if(v2n.find(v) != v2n.end()){
			idVertV = v2n.find(v)->second;
			idVertU = idV;
		   idV++;
	}
	else{
		idVertU= idV;
		idV++;
		idVertV= idV;
		idV++;
	}
	adj[idVertU].push_back(std::make_pair(idVertV, w));
	adj[idVertV].push_back(std::make_pair(idVertU, w));

	v2n.insert(pair<TCellID,int>(u,idVertU));
	v2n.insert(pair<TCellID,int>(v,idVertV));

	n2v.insert(pair<int,TCellID>(idVertU,u));
	n2v.insert(pair<int,TCellID>(idVertV,v));

	graphEdge2blockingEdge.insert(pair<TCellID,std::pair<int, int>>(idEdgeBlocking,std::make_pair(idVertU,idVertV)));

	weightEdge.insert(pair<TCellID ,double>(idEdgeBlocking,w));


}
void
DijkstraGraph::shortestPath(int src)
{
	// Create a set to store vertices that are being
	// processed
	set<std::pair<int, int>> setds;

	// Create a vector for distances and initialize all
	// distances as infinite (INF)
	vector<int> dist(V, INF);

	// Insert source itself in Set and initialize its
	// distance as 0.
	setds.insert(make_pair(0, src));
	dist[src] = 0;

	/* Looping till all shortest distance are finalized
	   then setds will become empty    */
	while (!setds.empty()) {
		// The first vertex in Set is the minimum distance
		// vertex, extract it from set.
		pair<int, int> tmp = *(setds.begin());
		setds.erase(setds.begin());

		// vertex label is stored in second of pair (it
		// has to be done this way to keep the vertices
		// sorted distance (distance must be first item
		// in pair)
		int u = tmp.second;

		// 'i' is used to get all adjacent vertices of a vertex
		list<pair<int, int>>::iterator i;
		for (i = adj[u].begin(); i != adj[u].end(); ++i) {
			// Get vertex label and weight of current adjacent
			// of u.
			int v = (*i).first;
			int weight = (*i).second;

			//    If there is shorter path to v through u.
			if (dist[v] > dist[u] + weight) {
				/*  If distance of v is not INF then it must be in
				    our set, so removing it and inserting again
				    with updated less distance.
				    Note : We extract only those vertices from Set
				    for which distance is finalized. So for them,
				    we would never reach here.  */
				if (dist[v] != INF) setds.erase(setds.find(make_pair(dist[v], v)));

				// Updating distance of v
				dist[v] = dist[u] + weight;
				setds.insert(make_pair(dist[v], v));
			}
		}
		printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

		for(auto s : setds){
			std::cout<<"1 : "<<s.first<<" & 2 : "<<s.second<<std::endl;
		}
	}

	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	// Print shortest distances stored in dist[]
	printf("Vertex   Distance from Source\n");
	for (int i = 0; i < V; ++i)
		printf("id node blocking :%u \n%d \t\t %d\n",n2v.find(i)->second, i, dist[i]);
}

/*----------------------------------------------------------------------------*/
// Allocates memory for adjacency list
GraphOrigin::GraphOrigin(int V)
{
	this->V = V;
	adj = new list< pair<int, int> >[V];
}

void GraphOrigin::addEdgeGraph(int u, int v, int w)
{
	adj[u].push_back(make_pair(v, w));
	adj[v].push_back(make_pair(u, w));
}

// Prints shortest paths from src to all other vertices
void GraphOrigin::shortestPathGraph(int src)
{
	// Create a set to store vertices that are being
	// processed
	set< pair<int, int> > setds;

	// Create a vector for distances and initialize all
	// distances as infinite (INF)
	vector<int> dist(V, INF);

	// Insert source itself in Set and initialize its
	// distance as 0.
	setds.insert(make_pair(0, src));
	dist[src] = 0;

	/* Looping till all shortest distance are finalized
	   then setds will become empty    */
	while (!setds.empty())
	{
		// The first vertex in Set is the minimum distance
		// vertex, extract it from set.
		pair<int, int> tmp = *(setds.begin());
		setds.erase(setds.begin());

		// vertex label is stored in second of pair (it
		// has to be done this way to keep the vertices
		// sorted distance (distance must be first item
		// in pair)
		int u = tmp.second;

		// 'i' is used to get all adjacent vertices of a vertex
		list< pair<int, int> >::iterator i;
		for (i = adj[u].begin(); i != adj[u].end(); ++i)
		{
			// Get vertex label and weight of current adjacent
			// of u.
			int v = (*i).first;
			int weight = (*i).second;

			//    If there is shorter path to v through u.
			if (dist[v] > dist[u] + weight)
			{
				/*  If distance of v is not INF then it must be in
				    our set, so removing it and inserting again
				    with updated less distance.
				    Note : We extract only those vertices from Set
				    for which distance is finalized. So for them,
				    we would never reach here.  */
				if (dist[v] != INF)
					setds.erase(setds.find(make_pair(dist[v], v)));

				// Updating distance of v
				dist[v] = dist[u] + weight;
				setds.insert(make_pair(dist[v], v));
			}
		}
	}

	// Print shortest distances stored in dist[]
	printf("Vertex   Distance from Source\n");
	for (int i = 0; i < V; ++i)
		printf("%d \t\t %d\n", i, dist[i]);
}
}
