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
	// Mark visited nodes as 1
	m_mesh_representation->newVariable<int,GMDS_NODE>("alreadyVisited");
	// Mark pushed nodes as 1
	m_mesh_representation->newVariable<int,GMDS_NODE>("alreadyPushed");
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
void QuantizationGraph::buildQuantizationSolution()
{
	std::cout<<" "<<std::endl;
	std::cout<<"====== Building quantization solution ====="<<std::endl;
	propagateFromRoots();
	addOnCycles();
	std::cout<<"====== Building quantization solution ====="<<std::endl;
	std::cout<<" "<<std::endl;
	// Test
//	auto sol = m_mesh_representation->getVariable<int,GMDS_NODE>("solution");
//	for (auto n_id:m_mesh_representation->nodes())
//		std::cout<<n_id<<" : "<<sol->value(n_id)<<std::endl;
}
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/