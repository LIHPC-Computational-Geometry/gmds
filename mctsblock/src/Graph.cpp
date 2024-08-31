/*----------------------------------------------------------------------------*/
#include <gmds/mctsblock/Graph.h>
#include <stdio.h>
#include <stdlib.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::mctsblock;

/*----------------------------------------------------------------------------*/
// Structure to represent a min heap node
struct MinHeapNode
{
	TCellID node_id;
	double dist;

	MinHeapNode(const TCellID AId, const double ADist)
	{
		node_id = AId;
		dist = ADist;
	}
};
/*----------------------------------------------------------------------------*/
// Structure to represent a min heap
struct MinHeap
{
	MinHeap(const int ACapacity);
	// A utility function to create a
	// new Min Heap Node
	void new_node(const TCellID AId, const double ADist);
	bool empty() const;
	std::shared_ptr<MinHeapNode> extractMin();
	// A standard function to
	// heapify at given idx
	// This function also updates
	// position of nodes when they are swapped.
	// Position is needed for decreaseKey()
	void minHeapify(const TCellID AIdx);
	void swapMinHeapNode(std::shared_ptr<MinHeapNode> *a, std::shared_ptr<MinHeapNode> *b);

	// A utility function to check if a given vertex
	// 'ANodeId' is in min heap or not
	bool isInMinHeap(const TCellID ANodeId);
	// Function to decrease key dist value of a given vertex v. This function
	// uses pos[] of min heap to get the current index of node in min heap
	void decreaseKey(const TCellID ANodeId, const double ADist);

	// Number of heap nodes present currently
	int size;
	// Capacity of min heap
	int capacity;
	// This is needed for decreaseKey()
	std::vector<TCellID> pos;
	std::vector<std::shared_ptr<MinHeapNode>> array;
};
/*----------------------------------------------------------------------------*/
void
MinHeap::new_node(const gmds::TCellID AId, const double ADist)
{
	array[AId] = std::make_shared<MinHeapNode>(AId, ADist);
	pos[AId] = AId;
}
/*----------------------------------------------------------------------------*/
// A utility function to create a Min Heap
MinHeap::MinHeap(const int ACapacity)
{
	pos.resize(ACapacity);
	array.resize(ACapacity);
	size = 0;
	capacity = ACapacity;
}
/*----------------------------------------------------------------------------*/
// A utility function to swap two
// nodes of min heap.
// Needed for min heapify
void
MinHeap::swapMinHeapNode(std::shared_ptr<MinHeapNode> *a, std::shared_ptr<MinHeapNode> *b)
{
	auto t = *a;
	*a = *b;
	*b = t;
}
/*----------------------------------------------------------------------------*/
void
MinHeap::minHeapify(const gmds::TCellID AIdx)
{
	int smallest, left, right;
	smallest = AIdx;
	left = 2 * AIdx + 1;
	right = 2 * AIdx + 2;

	if (left < size && array[left]->dist < array[smallest]->dist) smallest = left;

	if (right < size && array[right]->dist < array[smallest]->dist) smallest = right;

	if (smallest != AIdx) {
		// The nodes to be swapped in min heap
		auto smallestNode = array[smallest];
		auto idxNode = array[AIdx];

		// Swap positions
		pos[smallestNode->node_id] = AIdx;
		pos[idxNode->node_id] = smallest;

		// Swap nodes
		swapMinHeapNode(&array[smallest], &array[AIdx]);

		minHeapify(smallest);
	}
}
/*----------------------------------------------------------------------------*/
bool
MinHeap::empty() const
{
	return size == 0;
}
/*----------------------------------------------------------------------------*/
std::shared_ptr<MinHeapNode>
MinHeap::extractMin()
{
	if (empty()) return nullptr;
	// Store the root node
	auto root = array[0];

	// Replace root node with last node
	auto lastNode = array[size - 1];
	array[0] = lastNode;

	// Update position of last node
	pos[root->node_id] = size - 1;
	pos[lastNode->node_id] = 0;

	// Reduce heap size and heapify root
	--size;
	minHeapify(0);
	return root;
}
/*----------------------------------------------------------------------------*/
void
MinHeap::decreaseKey(const gmds::TCellID ANodeId, const double ADist)
{
	// Get the index of v in heap array
	auto i = pos[ANodeId];
	// Get the node and update its dist value
	array[i]->dist = ADist;

	// Travel up while the complete
	// tree is not heapified.
	// This is a O(Logn) loop
	while (i && array[i]->dist < array[(i - 1) / 2]->dist) {
		// Swap this node with its parent
		pos[array[i]->node_id] = (i - 1) / 2;
		pos[array[(i - 1) / 2]->node_id] = i;
		swapMinHeapNode(&array[i], &array[(i - 1) / 2]);

		// move to parent index
		i = (i - 1) / 2;
	}
}
/*----------------------------------------------------------------------------*/
bool
MinHeap::isInMinHeap(const gmds::TCellID ANodeId)
{
	return (pos[ANodeId] < size);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
AdjList::AdjList(std::shared_ptr<AdjListNode> h) : head(h) {}
/*----------------------------------------------------------------------------*/
AdjListNode::AdjListNode(gmds::TCellID ANodeID, double AWeight, std::shared_ptr<AdjListNode> AN)
  : dest(ANodeID), weight(AWeight), next(AN) {}
/*----------------------------------------------------------------------------*/
void
Graph::newAdjListNode(gmds::TCellID ADest, double AWeight, std::shared_ptr<AdjList> ANode)
{
	auto newNode = std::make_shared<AdjListNode>(ADest, AWeight, nullptr);
	newNode->next = ANode->head;
	ANode->head = newNode;
}
/*----------------------------------------------------------------------------*/
void
Graph::setWeightOneWay(const gmds::TCellID ASrc, const gmds::TCellID ADest, const double AWeight)
{
	auto adj_node = m_array[ASrc]->head;
	while (adj_node->dest != ADest)
		adj_node = adj_node->next;
	if (adj_node != nullptr && adj_node->dest == ADest) adj_node->weight = AWeight;
}
/*----------------------------------------------------------------------------*/
void
Graph::setWeight(const gmds::TCellID AN1, const gmds::TCellID AN2, const double AW)
{
	setWeightOneWay(m_input_to_graph_vertices[AN1], m_input_to_graph_vertices[AN2], AW);
	setWeightOneWay(m_input_to_graph_vertices[AN2], m_input_to_graph_vertices[AN1], AW);
}
/*----------------------------------------------------------------------------*/
Graph::Graph(const std::set<TCellID>& ANodeIDs)
{
	m_input_vertices = ANodeIDs;
	m_nb_vertices = m_input_vertices.size();
	m_graph_to_input_vertices.resize(m_nb_vertices);
	auto index=0;
	for (auto input_id: ANodeIDs){
		m_graph_to_input_vertices[index]=input_id;
		m_input_to_graph_vertices[input_id]=index;
		index++;
	}

	m_array.resize(m_nb_vertices);
	for (auto i = 0; i < m_nb_vertices; i++) {
		m_array[i] = std::make_shared<AdjList>();
	}
}
/*----------------------------------------------------------------------------*/
// Adds an edge to an undirected graph
void
Graph::addEdge(gmds::TCellID ASrcNode, gmds::TCellID ADestNode, double AWeight)
{
	// Add an edge from src to dest.
	// A new node is added to the adjacency
	// list of src. The node is
	// added at the beginning
	newAdjListNode(m_input_to_graph_vertices[ADestNode], AWeight, m_array[m_input_to_graph_vertices[ASrcNode]]);

	// Since graph is undirected,
	// add an edge from dest to src also
	newAdjListNode(m_input_to_graph_vertices[ASrcNode], AWeight, m_array[m_input_to_graph_vertices[ADestNode]]);
}
/*----------------------------------------------------------------------------*/
std::vector<std::pair<TCellID, double>>
Graph::getAdjNodes(const gmds::TCellID ANodeId)
{
	std::vector<std::pair<TCellID, double>> neigh;
	auto node = m_array[ANodeId]->head;
	while (node != nullptr) {
		neigh.push_back(std::make_pair(node->dest, node->weight));
		node = node->next;
	}
	return neigh;
}
/*----------------------------------------------------------------------------*/
void
Graph::computeDijkstra(gmds::TCellID ASrcNode)
{
	auto src_node = m_input_to_graph_vertices[ASrcNode];
	// dist values used to pick
	// minimum weight edge in cut
	m_dist.resize(m_nb_vertices);
	// minHeap represents the set of edges E
	MinHeap minHeap(m_nb_vertices);

	// Initialize min heap with all
	// vertices. dist value of all vertices
	for (auto v = 0; v < m_nb_vertices; ++v) {
		m_dist[v] = INT_MAX;
		minHeap.new_node(v, m_dist[v]);
	}

	// Make dist value of src vertex
	// as 0 so that it is extracted first
	minHeap.new_node(src_node, m_dist[src_node]);
	m_dist[src_node] = 0;
	minHeap.decreaseKey(src_node, m_dist[src_node]);

	// Initially size of min heap is equal to V
	minHeap.size = m_nb_vertices;

	// In the following loop, min heap contains all nodes
	// whose shortest distance is not yet finalized.
	while (!minHeap.empty()) {
		// Extract the vertex with minimum distance value
		auto minHeapNode = minHeap.extractMin();

		// Store the extracted vertex number
		int u = minHeapNode->node_id;

		// Traverse through all adjacent  vertices of u (the extracted
		// vertex) and update  their distance values
		auto pCrawl = m_array[u]->head;
		while (pCrawl != NULL) {
			int v = pCrawl->dest;

			// If shortest distance to v is  not finalized yet, and distance to v
			// through u is less than its  previously calculated distance
			if (minHeap.isInMinHeap(v) && m_dist[u] != INT_MAX && pCrawl->weight + m_dist[u] < m_dist[v]) {
				m_dist[v] = m_dist[u] + pCrawl->weight;
				// update distance value in min heap also
				minHeap.decreaseKey(v, m_dist[v]);
			}
			pCrawl = pCrawl->next;
		}
	}
	buildShortestPaths(src_node);
}
/*----------------------------------------------------------------------------*/
void
Graph::buildShortestPaths(const gmds::TCellID ASrc)
{
	// for each node different from src, we build the path
	m_shortest_path_nodes.clear();
	m_shortest_path_weight.clear();
	for (auto i = 0; i < m_nb_vertices; i++) {
		if (i == ASrc) continue;
		std::vector<TCellID> current_path;

		auto cur = i;
		//the second condition is here to ensure that we are not trying to
		//find a path between two unconnected parts of the graph

		while (cur != ASrc && current_path.size()<m_nb_vertices) {
			current_path.push_back(cur);
			auto neigh = getAdjNodes(cur);
			auto cur_dist = m_dist[cur];
			for (auto n_info : neigh) {
				if (n_info.second == cur_dist - m_dist[n_info.first]) cur = n_info.first;
			}
		}
		current_path.push_back(ASrc);
		m_shortest_path_nodes[i] = current_path;
		m_shortest_path_weight[i] = m_dist[i];
	}
}
/*----------------------------------------------------------------------------*/
std::map<TCellID, std::vector<TCellID>>
Graph::getShortestPath()
{
	// the shortest paths stored in m_shortest_path_nodes gives vector of ids
	// that correspond to the graph local numbering. We need to transpose
	// this numbering onto the input nodes
	std::map<TCellID, std::vector<TCellID>> paths;
	for (auto info: m_shortest_path_nodes){
		std::vector<TCellID> ids;
		for (auto i : info.second){
			ids.push_back(m_graph_to_input_vertices[i]);
		}
		paths[m_graph_to_input_vertices[info.first]] = ids;
	}
	return paths;
}
/*----------------------------------------------------------------------------*/
std::map<TCellID, double>
Graph::getShortestPathWeights()
{
	std::map<TCellID, double> paths;
	for (auto info: m_shortest_path_weight){
		paths[m_graph_to_input_vertices[info.first]] = info.second;
	}
	return paths;
}