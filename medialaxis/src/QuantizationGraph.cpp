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
	// Quantization solution on nodes of the graph
	m_mesh_representation->newVariable<int,GMDS_NODE>("solution");
	// Geometrical length of the half edge that the node represents
	m_mesh_representation->newVariable<double,GMDS_NODE>("geometrical_length");
	// Mark with 1 nodes where the quantization solution must be 0
	m_mesh_representation->newVariable<int,GMDS_NODE>("forbiden");
	// Quantization solution on edges of the graph (flux going through the edge)
	m_mesh_representation->newVariable<int,GMDS_EDGE>("flux");
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
void QuantizationGraph::setGeometricalLength(TCellID AID, double ALength)
{
	auto gl = m_mesh_representation->getVariable<double,GMDS_NODE>("geometrical_length");
	gl->set(AID,ALength);
}

/*----------------------------------------------------------------------------*/
void QuantizationGraph::markAsZero(TCellID AID)
{
	auto forbiden = m_mesh_representation->getVariable<int,GMDS_NODE>("forbiden");
	forbiden->set(AID,1);
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
Edge QuantizationGraph::getIntEdge(gmds::Node &AN)
{
	auto in_quad = m_mesh_representation->getVariable<int,GMDS_EDGE>("inQuad");
	Edge int_edge;
	for (auto e:AN.get<Edge>())
	{
		if (in_quad->value(e.id()) == 1)
		{
			int_edge = e;
			break;
		}
	}
	return int_edge;
}

/*----------------------------------------------------------------------------*/
std::vector<Edge> QuantizationGraph::getExtEdges(gmds::Node &AN)
{
	auto in_quad = m_mesh_representation->getVariable<int,GMDS_EDGE>("inQuad");
	std::vector<Edge> ext_edges;
	for (auto e:AN.get<Edge>())
	{
		if (in_quad->value(e.id()) == 0)
			ext_edges.push_back(e);
	}
	return ext_edges;
}

/*----------------------------------------------------------------------------*/
Node QuantizationGraph::getOtherNode(Node AN, Edge AE)
{
	if (AN.id() == AE.get<Node>()[0].id())
		return AE.get<Node>()[1];
	else 
		return AE.get<Node>()[0];
}

/*----------------------------------------------------------------------------*/
void QuantizationGraph::propagateFromRoot(TCellID AID)
{
	std::cout<<"> Propagating solution from root "<<AID<<std::endl;
	auto sol = m_mesh_representation->getVariable<int,GMDS_NODE>("solution");
	auto flux = m_mesh_representation->getVariable<int,GMDS_EDGE>("flux");
	auto geo_len = m_mesh_representation->getVariable<double,GMDS_NODE>("geometrical_length");
	auto forbiden = m_mesh_representation->getVariable<int,GMDS_NODE>("forbiden");
	Node root = m_mesh_representation->get<Node>(AID);
	// Check if the given node is indeed a root, and if its length is 0
	if (getExtEdges(root).empty() && sol->value(AID) == 0 && forbiden->value(AID) == 0)
	{
		// Build all the paths from the root to a leaf or an already seen node
		std::vector<std::vector<int>> paths;
		std::vector<int> path;
		std::vector<int> newPath;
		path.push_back(root.id());
		std::vector<int> alreadySeen(m_mesh_representation->getNbNodes());
		std::vector<int> newAlreadySeen(m_mesh_representation->getNbNodes());
		for (int i = 0; i < m_mesh_representation->getNbNodes(); i++)
			alreadySeen[i] = 0;
		alreadySeen[root.id()] = 1;
		int dir = 1; // Convention : when we continue the path trough an interior edges, we say we go in the positive direction.
		std::queue<std::vector<int>> path_to_continue;
		path_to_continue.push(path);
		std::queue<std::vector<int>> already_seen;
		already_seen.push(alreadySeen);
		std::queue<int> direction;
		direction.push(dir);
		while (!path_to_continue.empty())
		{
			path = path_to_continue.front();
			path_to_continue.pop();
			alreadySeen = already_seen.front();
			already_seen.pop();
			dir = direction.front();
			direction.pop();
			Node last_node = m_mesh_representation->get<Node>(path[path.size()-1]);
			if (dir == 1)
			{
				Edge e = getIntEdge(last_node);
				Node nxt = getOtherNode(last_node,e);
				if (alreadySeen[nxt.id()] == 0)
				{
					// We continue the path
					newPath = path;
					newPath.push_back(nxt.id());
					path_to_continue.push(newPath);
					newAlreadySeen = alreadySeen;
					newAlreadySeen[nxt.id()] = 1;
					already_seen.push(newAlreadySeen);
					direction.push(-1);
				}
			}
			if (dir == -1)
			{
				std::vector<Edge> extEdges = getExtEdges(last_node);
				if (extEdges.empty())
				{
					// Then this path is over, we reached an other root
					paths.push_back(path);
				}
				else
				{
					for (auto e:extEdges)
					{
						Node nxt = getOtherNode(last_node,e);
						if (alreadySeen[nxt.id()] == 0)
						{
							// We continue the path
							newPath = path;
							newPath.push_back(nxt.id());
							path_to_continue.push(newPath);
							newAlreadySeen = alreadySeen;
							newAlreadySeen[nxt.id()] = 1;
							already_seen.push(newAlreadySeen);
							direction.push(1);
						}
					}
				}
			}
		}
		// Find which path to increase
		int pos = -1;
		double maxMinLength = 0.;
		int maxNbZeros = 0;
		for (int i = 0; i < paths.size(); i ++)
		{
			double ok = true;
			int NbZeros = 0;
			double minLength = 1e6;
			for (auto id:paths[i])
			{
				if (sol->value(id) == 0)
					NbZeros += 1;
				if (geo_len->value(id)/double(sol->value(id)) < minLength)
					minLength = geo_len->value(id)/double(sol->value(id));
				if (forbiden->value(id) == 1)
					ok = false;
			}
			if (NbZeros > maxNbZeros && ok)
			{
				pos = i;
				maxNbZeros = NbZeros;
				maxMinLength = minLength;
			}
			if (NbZeros == maxNbZeros && minLength > maxMinLength && ok)
			{
				pos = i;
				maxMinLength = minLength;
			}
		}
		// Increase the path with the highest number of zeros
		if (pos >= 0)
		{
			for (int i = 0; i < paths[pos].size(); i++)
			{
				int id = paths[pos][i];
				sol->set(id,sol->value(id)+1);
				if (i < paths[pos].size()-1)
				{
					Node n1 = m_mesh_representation->get<Node>(id);
					Node n2 = m_mesh_representation->get<Node>(paths[pos][i+1]);
					Edge e = getCorrespondingEdge(n1,n2);
					flux->set(e.id(),flux->value(e.id())+1);
				}
			}
		}
		// for (auto i:paths[pos])
		// {
		// 	sol->set(i,sol->value(i)+1);
		// }
		// // Update the solution
		// for (auto p:paths)
		// {
		// 	// Check if increasing this path is useful
		// 	bool increase = false;
		// 	for (auto i:p)
		// 	{
		// 		if (sol->value(i) == 0)
		// 		{
		// 			increase = true;
		// 			break;
		// 		}
		// 	}
		// 	if (increase)
		// 	{
		// 		for (auto i:p)
		// 		{
		// 			sol->set(i,sol->value(i)+1);
		// 		}
		// 	}
		// }
	}
}

/*----------------------------------------------------------------------------*/
void QuantizationGraph::propagateFromRoots()
{
	//auto sol = m_mesh_representation->getVariable<int,GMDS_NODE>("solution");
	for (auto n_id:m_mesh_representation->nodes())
	{
		Node n = m_mesh_representation->get<Node>(n_id);
		if (getExtEdges(n).empty())
			propagateFromRoot(n_id);
	}
}

/*----------------------------------------------------------------------------*/
void QuantizationGraph::addOnCycles()
{
	auto sol = m_mesh_representation->getVariable<int,GMDS_NODE>("solution");
	auto forbiden = m_mesh_representation->getVariable<int,GMDS_NODE>("forbiden");
	for (auto n_id:m_mesh_representation->nodes())
	{
		if (sol->value(n_id) == 0 && forbiden->value(n_id) == 0)
		{
			std::cout<<"> Adding 1 to the solution to every vertex of the cycle of "<<n_id<<std::endl;
			// std::vector<Node> cycle = shortestCycle(n_id);
			// for (auto n:cycle)
			// 	sol->set(n.id(),sol->value(n.id())+1);
			increaseSolution(n_id);
		}
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
int QuantizationGraph::fluxValue(gmds::TCellID AID)
{
	auto flux = m_mesh_representation->getVariable<int,GMDS_EDGE>("flux");
	return flux->value(AID);
}

/*----------------------------------------------------------------------------*/
std::vector<Node> QuantizationGraph::shortestHalfPath(gmds::TCellID AID, int ADirection)
{
	auto forbiden = m_mesh_representation->getVariable<int,GMDS_NODE>("forbiden");
	// Get the node
	Node n0 = m_mesh_representation->get<Node>(AID);
	// Initialize queues and paths
	std::queue<std::vector<Node>> paths;
	std::queue<std::vector<int>> alreadySeen;
	std::queue<int> direction;
	std::vector<Node> path;
	std::vector<int> already_seen(m_mesh_representation->getNbNodes());
	for (int i = 0; i < m_mesh_representation->getNbNodes(); i++)
		already_seen[i] = 0;
	path.push_back(n0);
	already_seen[n0.id()] = 1;
	paths.push(path);
	alreadySeen.push(already_seen);
	int dir = ADirection;
	direction.push(dir);
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
		dir = direction.front();
		direction.pop();
		// Get the current last node of the path
		front = path[path.size()-1];
		if (dir == 1)
		{
			Edge nxt_edge = getIntEdge(front);
			Node n = getOtherNode(front,nxt_edge);
			if (already_seen[n.id()] == 0 && forbiden->value(n.id()) == 0)
			{
				// We continue the path
				new_path = path;
				new_path.push_back(n);
				paths.push(new_path);
				new_already_seen = already_seen;
				new_already_seen[n.id()] = 1;
				alreadySeen.push(new_already_seen);
				direction.push(-1);
			}
		}
		else
		{
			leaving_edges = getExtEdges(front);
			if (leaving_edges.empty())
			{
				// It means we have reached a root
				shortestPath = path;
				Continue = false;
			}
			else
			{
				for (auto e:leaving_edges)
				{
					Node n = getOtherNode(front,e);
					if (already_seen[n.id()] == 0 && forbiden->value(n.id()) == 0)
					{
						// Then the path continues
						new_path = path;
						new_path.push_back(n);
						paths.push(new_path);
						new_already_seen = already_seen;
						new_already_seen[n.id()] = 1;
						alreadySeen.push(new_already_seen);
						direction.push(1);
					}
				}
			}
		}
	}
	if (ADirection == -1)
		std::reverse(shortestPath.begin(),shortestPath.end());
	return shortestPath;
}

/*----------------------------------------------------------------------------*/
std::vector<Node> QuantizationGraph::shortestCycle(TCellID AID)
{
	auto forbiden = m_mesh_representation->getVariable<int,GMDS_NODE>("forbiden");
	Node n0 = m_mesh_representation->get<Node>(AID);
	std::queue<std::vector<Node>> paths;
	std::queue<std::vector<int>> alreadySeen;
	std::queue<int> direction;
	std::vector<Node> path;
	std::vector<int> already_seen(m_mesh_representation->getNbNodes());
	for (int i = 0; i < m_mesh_representation->getNbNodes(); i++)
		already_seen[i] = 0;
	path.push_back(n0);
	already_seen[n0.id()] = 1;
	paths.push(path);
	alreadySeen.push(already_seen);
	direction.push(1);
	int dir;
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
		dir = direction.front();
		direction.pop();
		front = path[path.size()-1];
		if (dir == 1)
		{
			Edge nxt_edge = getIntEdge(front);
			Node n = getOtherNode(front,nxt_edge);
			if (n.id() == n0.id())
			{
				Continue = false;
				shortestCycle = path;
				break;
			}
			else 
			{
				if (already_seen[n.id()] == 0 && forbiden->value(n.id()) == 0)
				{
					// Then the path continues
					new_path = path;
					new_path.push_back(n);
					paths.push(new_path);
					new_already_seen = already_seen;
					new_already_seen[n.id()] = 1;
					alreadySeen.push(new_already_seen);
					direction.push(-1);
				}
			}
		}
		else
		{
			leaving_edges = getExtEdges(front);
			for (auto e:leaving_edges)
			{
				Node n = getOtherNode(front,e);
				if (n.id() == n0.id())
				{
					Continue = false;
					shortestCycle = path;
					break;
				}
				else
				{
					if (already_seen[n.id()] == 0 && forbiden->value(n.id()) == 0)
					{
						// Then the path continues
						new_path = path;
						new_path.push_back(n);
						paths.push(new_path);
						new_already_seen = already_seen;
						new_already_seen[n.id()] = 1;
						alreadySeen.push(new_already_seen);
						direction.push(1);
					}
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
	std::vector<Node> fromNode = shortestHalfPath(AID,1);
	std::vector<Node> toNode = shortestHalfPath(AID,-1);
	std::vector<Node> cycle = shortestCycle(AID);
	if (fromNode.empty() || toNode.empty())
		shortestPath = cycle;
	else
	{
		if ((fromNode.size()+toNode.size()-1 <= cycle.size()) || cycle.empty())
		{
			shortestPath = toNode;
			for (int i = 1; i < fromNode.size(); i++)
				shortestPath.push_back(fromNode[i]);
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
	std::queue<int> direction;
	std::vector<Node> path;
	std::vector<int> already_seen(m_mesh_representation->getNbNodes());
	for (int i = 0; i < m_mesh_representation->getNbNodes(); i++)
		already_seen[i] = 0;
	path.push_back(n0);
	already_seen[n0.id()] = 1;
	paths.push(path);
	paths.push(path);
	alreadySeen.push(already_seen);
	alreadySeen.push(already_seen);
	direction.push(1);
	direction.push(-1);
	bool Continue = true;
	std::vector<Node> problematicPath;
	std::vector<Edge> leaving_edges;
	Node front,lastNode;
	std::vector<Node> new_path;
	std::vector<int> new_already_seen(m_mesh_representation->getNbNodes());
	int dir;
	while (Continue && !paths.empty())
	{
		path = paths.front();
		paths.pop();
		already_seen = alreadySeen.front();
		alreadySeen.pop();
		dir = direction.front();
		direction.pop();
		front = path[path.size()-1];
		leaving_edges.clear();
		if (dir == 1)
		{
			leaving_edges.push_back(getIntEdge(front));
		}
		else
		{
			leaving_edges = getExtEdges(front);
		}
		for (auto e:leaving_edges)
		{
			Node n = getOtherNode(front,e);
			if (already_seen[n.id()] == 0)
			{
				// Then the path continues
				new_path = path;
				new_path.push_back(n);
				paths.push(new_path);
				new_already_seen = already_seen;
				new_already_seen[n.id()] = 1;
				alreadySeen.push(new_already_seen);
				direction.push(-dir);
			}
			else
			{
				// Then we found a cycle, the path stops
				Continue = false;
				problematicPath = path;
				lastNode = n;
				break;
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
	for (auto e1:AN1.get<Edge>())
	{
		if (e1.get<Node>()[1].id() == AN2.id() || e1.get<Node>()[0].id() == AN2.id())
		{
			e = e1;
			break;
		}
	}
	return e;
}

/*----------------------------------------------------------------------------*/
std::vector<Node> QuantizationGraph::middleQuadEdge(std::vector<Node> AV)
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
		return {n1,n2};
	else
	{
		n2 = n1;
		if (i > 0)
			n1 = AV[i-1];
		else
			n1 = AV[AV.size()-1];
		return {n1,n2};
	}
}

/*----------------------------------------------------------------------------*/
bool QuantizationGraph::increaseSolution(TCellID AID)
{
	auto sol = m_mesh_representation->getVariable<int,GMDS_NODE>("solution");
	auto flux = m_mesh_representation->getVariable<int,GMDS_EDGE>("flux");
	std::vector<Node> path = shortestElementaryPath(AID);
	if (path.empty())
		return false;
	for (int i = 0; i < path.size(); i++)
	{
		Node n1 = path[i];
		sol->set(n1.id(),sol->value(n1.id())+1);
		if (i < path.size()-1)
		{
			Node n2 = path[i+1];
			Edge e = getCorrespondingEdge(n1,n2);
			flux->set(e.id(),flux->value(e.id())+1);
		}
		else
		{
			// Check if the path is a cycle
			Node n2 = path[0];
			bool areLinked = false;
			for (auto e:n1.get<Edge>())
			{
				if (e.get<Node>()[0].id() == n2.id() || e.get<Node>()[1].id() == n2.id())
				{
					areLinked = true;
					break;
				}
			}
			if (areLinked)
			{
				Edge e = getCorrespondingEdge(n1,n2);
				flux->set(e.id(),flux->value(e.id())+1);
			}
		}
	}
	// for (auto n:path)
	// 	sol->set(n.id(),sol->value(n.id())+1);
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
				std::vector<Node> probNodes = middleQuadEdge(probCycle);
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
				std::vector<Node> probNodes = middleQuadEdge(probCycle);
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
			std::vector<Node> probNodes = middleQuadEdge(probCycle);
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
			std::vector<Node> probNodes = middleQuadEdge(probCycle);
			problematicCouplesOfNodes.push_back(probNodes);
		}
	}
	return problematicCouplesOfNodes;
}

/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/