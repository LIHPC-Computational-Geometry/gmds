//
// Created by chenyt on 27/09/24.
//

#include "gmds/medialaxis/QuantizationGraph.h"
#include "gmds/ig/MeshDoctor.h"
#include <queue>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
QuantizationGraph::QuantizationGraph()
{
	m_mesh_representation = new Mesh(MeshModel(DIM3 | E | N |
	                                           E2N | N2E));
	// Mark vwith 1 visited nodes in the graph building process
	m_mesh_representation->newVariable<int,GMDS_NODE>("alreadyVisited");
	// Mark with 1 pushed nodes in the graph building process
	m_mesh_representation->newVariable<int,GMDS_NODE>("alreadyPushed");
	// Mark with 1 graph edges connecting two half-edges belonging to the same quad
	m_mesh_representation->newVariable<int,GMDS_EDGE>("inQuad");
	// Quantization solution
	m_mesh_representation->newVariable<int,GMDS_NODE>("solution");
}

/*----------------------------------------------------------------------------*/
Node QuantizationGraph::newNode()
{
	Node n = m_mesh_representation->newNode();
	return n;
}

/*----------------------------------------------------------------------------*/
Edge QuantizationGraph::newEdge(gmds::TCellID AId1, gmds::TCellID AId2)
{
	Edge e = m_mesh_representation->newEdge(AId1,AId2);
	return e;
}

/*----------------------------------------------------------------------------*/
void QuantizationGraph::display()
{
	for (auto e_id:m_mesh_representation->edges())
	{
		Edge e = m_mesh_representation->get<Edge>(e_id);
		std::cout<<"Edge "<<e_id<<" : ("<<e.get<Node>()[0].id()<<","<<e.get<Node>()[1].id()<<")"<<std::endl;
	}
}

/*----------------------------------------------------------------------------*/
void QuantizationGraph::markAsVisited(gmds::TCellID AID)
{
	auto visited = m_mesh_representation->getVariable<int,GMDS_NODE>("alreadyVisited");
	visited->set(AID,1);
}

/*----------------------------------------------------------------------------*/
int QuantizationGraph::alreadyVisited(gmds::TCellID AID)
{
	auto visited = m_mesh_representation->getVariable<int,GMDS_NODE>("alreadyVisited");
	return  visited->value(AID);
}

/*----------------------------------------------------------------------------*/
void QuantizationGraph::markAsPushed(gmds::TCellID AID)
{
	auto visited = m_mesh_representation->getVariable<int,GMDS_NODE>("alreadyPushed");
	visited->set(AID,1);
}

/*----------------------------------------------------------------------------*/
int QuantizationGraph::alreadyPushed(gmds::TCellID AID)
{
	auto visited = m_mesh_representation->getVariable<int,GMDS_NODE>("alreadyPushed");
	return  visited->value(AID);
}

/*----------------------------------------------------------------------------*/
void QuantizationGraph::markAsInQuad(TCellID AID)
{
	auto inQuad = m_mesh_representation->getVariable<int,GMDS_EDGE>("inQuad");
	inQuad->set(AID,1);
}

/*----------------------------------------------------------------------------*/
void QuantizationGraph::updateConnectivity()
{
	std::cout<<"> Updating quantization graph connectivity"<<std::endl;
	MeshDoctor doc(m_mesh_representation);
	doc.updateUpwardConnectivity();
}

/*----------------------------------------------------------------------------*/
std::vector<Edge> QuantizationGraph::getLeavingEdges(gmds::Node &AN)
{
	std::vector<Edge> leaving_edges;
	for (auto e:AN.get<Edge>())
	{
		if (e.get<Node>()[0].id() == AN.id())
			leaving_edges.push_back(e);
	}
	return leaving_edges;
}

/*----------------------------------------------------------------------------*/
std::vector<Edge> QuantizationGraph::getEnteringEdges(gmds::Node &AN)
{
	// m2.newNode(1.,2.);
	// m2.newNode(3.,2.);
	// m2.newNode(4.,2.);
	// m2.newNode(0.,3.);

	std::vector<Edge> entering_edges;
	for (auto e:AN.get<Edge>())
	{
		if (e.get<Node>()[1].id() == AN.id())
			entering_edges.push_back(e);
	}
	return entering_edges;
}

/*----------------------------------------------------------------------------*/
void QuantizationGraph::propagateFromRoot(TCellID AID)
{
	std::cout<<"> Propagating solution from root "<<AID<<std::endl;
	Node root = m_mesh_representation->get<Node>(AID);
	// Check if the given node is indeed a root
	if (getEnteringEdges(root).empty())
	{
		auto sol = m_mesh_representation->getVariable<int,GMDS_NODE>("solution");
		// Build all the paths from the root to a leaf or an already seen node
		std::vector<std::vector<int>> paths;
		std::vector<int> path;
		std::vector<int> newPath;
		path.push_back(root.id());
		std::vector<int> alreadySeen(m_mesh_representation->getNbNodes());
		std::vector<int> newAlreadySeen(m_mesh_representation->getNbNodes());
		for (int i = 0; i < m_mesh_representation->getNbNodes(); i++)
			alreadySeen[i] = 0;
		std::queue<std::vector<int>> path_to_continue;
		path_to_continue.push(path);
		std::queue<std::vector<int>> already_seen;
		already_seen.push(alreadySeen);
		while (!path_to_continue.empty())
		{
			path = path_to_continue.front();
			path_to_continue.pop();
			alreadySeen = already_seen.front();
			already_seen.pop();
			// Continue the path
			Node n1,n2,n3;
			n1 = m_mesh_representation->get<Node>(path[path.size()-1]);
			while(!getLeavingEdges(n1).empty() && alreadySeen[n1.id()] == 0)
			{
				alreadySeen[n1.id()] = 1;
				for (auto e: getLeavingEdges(n1))
				{
					n2 = e.get<Node>()[1];
					if (alreadySeen[n2.id()] == 0)
					{
						path.push_back(n2.id());
						break;
					}
				}
				for (auto e: getLeavingEdges(n1))
				{
					n3 = e.get<Node>()[1];
					if (n3.id() != n2.id() && alreadySeen[n3.id()] == 0)
					{
						newPath = path;
						newPath.erase(newPath.begin()+newPath.size()-1);
						newPath.push_back(n3.id());
						newAlreadySeen = alreadySeen;
						path_to_continue.push(newPath);
						already_seen.push(newAlreadySeen);
					}
				}
				n1 = n2;
			}
			paths.push_back(path);
		}
		// Update the solution
		for (auto p:paths)
		{
			Node end = m_mesh_representation->get<Node>(p[p.size()-1]);
			if (getLeavingEdges(end).empty())
			{
				for (auto i:p)
					sol->set(i,sol->value(i)+1);
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
void QuantizationGraph::propagateFromRoots()
{
	for (auto n_id:m_mesh_representation->nodes())
	{
		Node n = m_mesh_representation->get<Node>(n_id);
		if (getEnteringEdges(n).empty())
			propagateFromRoot(n_id);
	}
}

/*----------------------------------------------------------------------------*/
void QuantizationGraph::addOnCycle(gmds::TCellID AID)
{
	std::cout<<"> Adding 1 to the solution to every vertex of the cycle of "<<AID<<std::endl;
	Node n = m_mesh_representation->get<Node>(AID);
	auto sol = m_mesh_representation->getVariable<int,GMDS_NODE>("solution");
	// Build all the paths from the node to a leaf or an already seen node
	std::vector<std::vector<int>> paths;
	std::vector<int> path;
	std::vector<int> newPath;
	path.push_back(n.id());
	std::vector<int> alreadySeen(m_mesh_representation->getNbNodes());
	std::vector<int> newAlreadySeen(m_mesh_representation->getNbNodes());
	for (int i = 0; i < m_mesh_representation->getNbNodes(); i++)
		alreadySeen[i] = 0;
	std::queue<std::vector<int>> path_to_continue;
	path_to_continue.push(path);
	std::queue<std::vector<int>> already_seen;
	already_seen.push(alreadySeen);
	while (!path_to_continue.empty())
	{
		path = path_to_continue.front();
		path_to_continue.pop();
		alreadySeen = already_seen.front();
		already_seen.pop();
		// Continue the path
		Node n1,n2,n3;
		n1 = m_mesh_representation->get<Node>(path[path.size()-1]);
		while(!getLeavingEdges(n1).empty() && alreadySeen[n1.id()] == 0)
		{
			alreadySeen[n1.id()] = 1;
			bool add = false;
			for (auto e: getLeavingEdges(n1))
			{
				n2 = e.get<Node>()[1];
				if (alreadySeen[n2.id()] == 0)
				{
					path.push_back(n2.id());
					break;
				}
			}
			for (auto e: getLeavingEdges(n1))
			{
				n3 = e.get<Node>()[1];
				if (n3.id() != n2.id() && alreadySeen[n3.id()] == 0)
				{
					newPath = path;
					newPath.erase(newPath.begin()+newPath.size()-1);
					newPath.push_back(n3.id());
					newAlreadySeen = alreadySeen;
					path_to_continue.push(newPath);
					already_seen.push(newAlreadySeen);
				}
			}
			n1 = n2;
		}
		paths.push_back(path);
	}
	// Update the solution
	for (auto p:paths)
	{
		Node end = m_mesh_representation->get<Node>(p[p.size()-1]);
		bool cycle = false;
		for (auto e: getLeavingEdges(end))
		{
			if (e.get<Node>()[1].id() == n.id())
			{
				cycle = true;
				break;
			}
		}
		if (cycle)
		{
			for (auto i:p)
			{
				sol->set(i,sol->value(i)+1);
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
void QuantizationGraph::addOnCycles()
{
	auto sol = m_mesh_representation->getVariable<int,GMDS_NODE>("solution");
	for (auto n_id:m_mesh_representation->nodes())
	{
		if (sol->value(n_id) == 0)
			addOnCycle(n_id);
	}
}

/*----------------------------------------------------------------------------*/
void QuantizationGraph::buildCompleteSolution()
{
	std::cout << " " << std::endl;
	std::cout << "====== Building quantization solution =====" << std::endl;
	propagateFromRoots();
	addOnCycles();
	std::cout << "===========================================" << std::endl;
	std::cout << " " << std::endl;
}

/*----------------------------------------------------------------------------*/
void QuantizationGraph::displaySolution()
{
	auto sol = m_mesh_representation->getVariable<int,GMDS_NODE>("solution");
	for (auto n_id:m_mesh_representation->nodes())
		std::cout<<"Length of half edge "<<n_id<<" : "<<sol->value(n_id)<<std::endl;
}

/*----------------------------------------------------------------------------*/
int QuantizationGraph::quantizationSolutionValue(gmds::TCellID AID)
{
	auto sol = m_mesh_representation->getVariable<int,GMDS_NODE>("solution");
	return sol->value(AID);
}

/*----------------------------------------------------------------------------*/
std::vector<Node> QuantizationGraph::shortestPathToASink(gmds::TCellID AID)
{
	// Get the node
	Node n0 = m_mesh_representation->get<Node>(AID);
	// Initialize queues and paths
	std::queue<std::vector<Node>> paths;
	std::queue<std::vector<int>> alreadySeen;
	std::vector<Node> path;
	std::vector<int> already_seen(m_mesh_representation->getNbNodes());
	for (int i = 0; i < m_mesh_representation->getNbNodes(); i++)
		already_seen[i] = 0;
	path.push_back(n0);
	already_seen[n0.id()] = 1;
	paths.push(path);
	alreadySeen.push(already_seen);
	bool Continue = true;
	std::vector<Node> shortestPath;
	std::vector<Edge> leaving_edges;
	Node front;
	std::vector<Node> new_path;
	std::vector<int> new_already_seen(m_mesh_representation->getNbNodes());
	// Build the paths
	while (Continue && !paths.empty())
	{
		// Get the first path to be continued
		path = paths.front();
		paths.pop();
		already_seen = alreadySeen.front();
		alreadySeen.pop();
		// Get the current last node of the path
		front = path[path.size()-1];
		leaving_edges = getLeavingEdges(front);
		if (leaving_edges.empty())
		{
			// It means we have reached a sink
			shortestPath = path;
			Continue = false;
		}
		else
		{
			for (auto e:leaving_edges)
			{
				Node n = e.get<Node>()[1];
				if (already_seen[n.id()] == 0)
				{
					// Then the path continues
					new_path = path;
					new_path.push_back(n);
					paths.push(new_path);
					new_already_seen = already_seen;
					new_already_seen[n.id()] = 1;
					alreadySeen.push(new_already_seen);
				}
			}
		}
	}

	return shortestPath;
}

/*----------------------------------------------------------------------------*/
std::vector<Node> QuantizationGraph::shortestPathFromASource(TCellID AID)
{
	Node n0 = m_mesh_representation->get<Node>(AID);
	std::queue<std::vector<Node>> paths;
	std::queue<std::vector<int>> alreadySeen;
	std::vector<Node> path;
	std::vector<int> already_seen(m_mesh_representation->getNbNodes());
	for (int i = 0; i < m_mesh_representation->getNbNodes(); i++)
		already_seen[i] = 0;
	path.push_back(n0);
	already_seen[n0.id()] = 1;
	paths.push(path);
	alreadySeen.push(already_seen);
	bool Continue = true;
	std::vector<Node> shortestPath;
	std::vector<Edge> entering_edges;
	Node front;
	std::vector<Node> new_path;
	std::vector<int> new_already_seen(m_mesh_representation->getNbNodes());
	while (Continue && !paths.empty())
	{
		path = paths.front();
		paths.pop();
		already_seen = alreadySeen.front();
		alreadySeen.pop();
		front = path[path.size()-1];
		entering_edges = getEnteringEdges(front);
		if (entering_edges.empty())
		{
			// It means we have reached a source
			shortestPath = path;
			Continue = false;
		}
		else
		{
			for (auto e:entering_edges)
			{
				Node n = e.get<Node>()[0];
				if (already_seen[n.id()] == 0)
				{
					// Then the path continues
					new_path = path;
					new_path.push_back(n);
					paths.push(new_path);
					new_already_seen = already_seen;
					new_already_seen[n.id()] = 1;
					alreadySeen.push(new_already_seen);
				}
			}
		}
	}
	std::reverse(shortestPath.begin(),shortestPath.end());
	return shortestPath;
}

/*----------------------------------------------------------------------------*/
std::vector<Node> QuantizationGraph::shortestCycle(TCellID AID)
{
	Node n0 = m_mesh_representation->get<Node>(AID);
	std::queue<std::vector<Node>> paths;
	std::queue<std::vector<int>> alreadySeen;
	std::vector<Node> path;
	std::vector<int> already_seen(m_mesh_representation->getNbNodes());
	for (int i = 0; i < m_mesh_representation->getNbNodes(); i++)
		already_seen[i] = 0;
	path.push_back(n0);
	already_seen[n0.id()] = 1;
	paths.push(path);
	alreadySeen.push(already_seen);
	bool Continue = true;
	std::vector<Node> shortestCycle;
	std::vector<Edge> leaving_edges;
	Node front;
	std::vector<Node> new_path;
	std::vector<int> new_already_seen(m_mesh_representation->getNbNodes());
	while (Continue && !paths.empty())
	{
		path = paths.front();
		paths.pop();
		already_seen = alreadySeen.front();
		alreadySeen.pop();
		front = path[path.size()-1];
		leaving_edges = getLeavingEdges(front);
		for (auto e:leaving_edges)
		{
			Node n = e.get<Node>()[1];
			if (n.id() == n0.id())
			{
				Continue = false;
				shortestCycle = path;
			}
			else
			{
				if (already_seen[n.id()] == 0)
				{
					// Then the path continues
					new_path = path;
					new_path.push_back(n);
					paths.push(new_path);
					new_already_seen = already_seen;
					new_already_seen[n.id()] = 1;
					alreadySeen.push(new_already_seen);
				}
			}
		}
	}
	return shortestCycle;
}

/*----------------------------------------------------------------------------*/
std::vector<Node> QuantizationGraph::shortestElementaryPath(TCellID AID)
{
	std::vector<Node> shortestPath;
	std::vector<Node> node2sink = shortestPathToASink(AID);
	std::vector<Node> source2node = shortestPathFromASource(AID);
	std::vector<Node> cycle = shortestCycle(AID);
	if (node2sink.empty() || source2node.empty())
		shortestPath = cycle;
	else
	{
		if ((node2sink.size()+source2node.size()-1 <= cycle.size()) || cycle.empty())
		{
			shortestPath = source2node;
			for (int i = 1; i < node2sink.size(); i++)
				shortestPath.push_back(node2sink[i]);
		}
		else
			shortestPath = cycle;
	}
	return shortestPath;
}

/*----------------------------------------------------------------------------*/
std::vector<Node> QuantizationGraph::problematicCycle(TCellID AID)
{
	Node n0 = m_mesh_representation->get<Node>(AID);
	std::queue<std::vector<Node>> paths;
	std::queue<std::vector<int>> alreadySeen;
	std::vector<Node> path;
	std::vector<int> already_seen(m_mesh_representation->getNbNodes());
	for (int i = 0; i < m_mesh_representation->getNbNodes(); i++)
		already_seen[i] = 0;
	path.push_back(n0);
	already_seen[n0.id()] = 1;
	paths.push(path);
	alreadySeen.push(already_seen);
	bool Continue = true;
	std::vector<Node> problematicPath;
	std::vector<Edge> leaving_edges;
	Node front,lastNode;
	std::vector<Node> new_path;
	std::vector<int> new_already_seen(m_mesh_representation->getNbNodes());
	while (Continue && !paths.empty())
	{
		path = paths.front();
		paths.pop();
		already_seen = alreadySeen.front();
		alreadySeen.pop();
		front = path[path.size()-1];
		leaving_edges = getLeavingEdges(front);
		for (auto e:leaving_edges)
		{
			Node n = e.get<Node>()[1];
			if (already_seen[n.id()] == 0)
			{
				// Then the path continues
				new_path = path;
				new_path.push_back(n);
				paths.push(new_path);
				new_already_seen = already_seen;
				new_already_seen[n.id()] = 1;
				alreadySeen.push(new_already_seen);
			}
			else
			{
				// Then we found a cycle, the path stops
				Continue = false;
				problematicPath = path;
				lastNode = n;
			}
			
		}
	}
	// Extract the problematic cycle from the problematic path
	std::vector<Node> problematicCycle;
	bool add = false;
	for (auto n:problematicPath)
	{
		if (n.id() == lastNode.id())
			add =true;
		if (add)
			problematicCycle.push_back(n);
	}
	return problematicCycle;
}

/*----------------------------------------------------------------------------*/
Edge QuantizationGraph::getCorrespondingEdge(Node AN1, Node AN2)
{
	Edge e;
	for (auto e1:getLeavingEdges(AN1))
	{
		if (e1.get<Node>()[1].id() == AN2.id())
		{
			e = e1;
			break;
		}
	}
	return e;
}

/*----------------------------------------------------------------------------*/
Edge QuantizationGraph::middleQuadEdge(std::vector<Node> AV)
{
	auto inQuad = m_mesh_representation->getVariable<int,GMDS_EDGE>("inQuad");
	int i = int(double(AV.size())/2.);
	Node n1 = AV[i];
	Node n2;
	if (i < AV.size()-1)
		n2 = AV[i+1];
	else 
		n2 = AV[0];
	Edge e = getCorrespondingEdge(n1,n2);
	if (inQuad->value(e.id()))
		return e;
	else
	{
		n2 = n1;
		if (i > 0)
			n1 = AV[i-1];
		else
			n1 = AV[AV.size()-1];
		e = getCorrespondingEdge(n1,n2);
		return e;
	}
}

/*----------------------------------------------------------------------------*/
bool QuantizationGraph::increaseSolution(TCellID AID)
{
	auto sol = m_mesh_representation->getVariable<int,GMDS_NODE>("solution");
	std::vector<Node> path = shortestElementaryPath(AID);
	if (path.empty())
		return false;
	for (auto n:path)
		sol->set(n.id(),sol->value(n.id())+1);
	return true;
}

/*----------------------------------------------------------------------------*/
std::vector<std::vector<Node>> QuantizationGraph::buildMinimalSolution()
{
	std::vector<std::vector<Node>> problematicCouplesOfNodes;
	std::cout<<"> Building minimal quantization solution"<<std::endl;
	bool success;
	auto sol = m_mesh_representation->getVariable<int,GMDS_NODE>("solution");
	// Ensure non zero conditions
	for (auto n_id:m_non_zero_verticies)
	{
		if (sol->value(n_id) == 0)
		{
			success = increaseSolution(n_id);
			if (!success)
			{
				std::vector<Node> probCycle = problematicCycle(n_id);
				Edge probEdge = middleQuadEdge(probCycle);
				std::vector<Node> probNodes = {probEdge.get<Node>()[0],probEdge.get<Node>()[1]};
				problematicCouplesOfNodes.push_back(probNodes);
			}
		}	
	}
	// Ensure condition "at least one non zero in each group"
	for (auto group:m_non_zero_groups)
	{
		bool ok = false;
		for (auto n_id:group)
		{
			if (sol->value(n_id) > 0)
				ok = true;
		}
		if (!ok)
		{
			success = false;
			for (auto n:group)
			{
				success = increaseSolution(n);
				if (success)
					break;
			}
			if (!success)
			{
				std::vector<Node> probCycle = problematicCycle(group[0]);
				Edge probEdge = middleQuadEdge(probCycle);
				std::vector<Node> probNodes = {probEdge.get<Node>()[0],probEdge.get<Node>()[1]};
				problematicCouplesOfNodes.push_back(probNodes);
			}
		}
	}
	return problematicCouplesOfNodes;
}

/*----------------------------------------------------------------------------*/
void QuantizationGraph::roughlyRepairSolution()
{
	std::cout<<"> Roughly repairing solution"<<std::endl;
	auto sol = m_mesh_representation->getVariable<int,GMDS_NODE>("solution");
	auto inQuad = m_mesh_representation->getVariable<int,GMDS_EDGE>("inQuad");
	for (auto e_id:m_mesh_representation->edges())
	{
		if (inQuad->value(e_id) == 1)
		{
			Edge e = m_mesh_representation->get<Edge>(e_id);
			Node n1 = e.get<Node>()[0];
			Node n2 = e.get<Node>()[1];
			int l1 = sol->value(n1.id());
			int l2 = sol->value(n2.id());
			if (l1 != l2)
			{
				std::cout<<"... changing a length"<<std::endl;
				if (l1 < l2)
					sol->set(n2.id(),l1);
				if (l2 < l1)
					sol->set(n1.id(),l2);
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
std::vector<std::vector<Node>> QuantizationGraph::checkSolutionValidity()
{
	std::vector<std::vector<Node>> problematicCouplesOfNodes;
	std::cout<<"> Checking validity of the quantization solution"<<std::endl;
	bool success;
	auto sol = m_mesh_representation->getVariable<int,GMDS_NODE>("solution");
	// Ensure non zero conditions
	for (auto n_id:m_non_zero_verticies)
	{
		if (sol->value(n_id) == 0)
		{
			std::vector<Node> probCycle = problematicCycle(n_id);
			Edge probEdge = middleQuadEdge(probCycle);
			std::vector<Node> probNodes = {probEdge.get<Node>()[0],probEdge.get<Node>()[1]};
			problematicCouplesOfNodes.push_back(probNodes);
		}	
	}
	// Ensure condition "at least one non zero in each group"
	for (auto group:m_non_zero_groups)
	{
		bool ok = false;
		for (auto n_id:group)
		{
			if (sol->value(n_id) > 0)
				ok = true;
		}
		if (!ok)
		{
			std::vector<Node> probCycle = problematicCycle(group[0]);
			Edge probEdge = middleQuadEdge(probCycle);
			std::vector<Node> probNodes = {probEdge.get<Node>()[0],probEdge.get<Node>()[1]};
			problematicCouplesOfNodes.push_back(probNodes);
		}
	}
	return problematicCouplesOfNodes;
}

/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/