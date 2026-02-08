//
// Created by chenyt on 26/09/24.
//

#include "gmds/medialaxis/QuantizationSolver.h"
#include <queue>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
QuantizationSolver::QuantizationSolver(gmds::Mesh &AMesh, double AMeshSize)
{
	// Non conformal quad mesh
	m_mesh = &AMesh;
	// Correspondence edge/half edge
	m_mesh->newVariable<int,GMDS_EDGE>("edge2halfEdge1");
	m_mesh->newVariable<int,GMDS_EDGE>("edge2halfEdge2");
	// Correspondance edge/quantization graph edge
	m_mesh->newVariable<int,GMDS_EDGE>("edge2qgrapheEdge");
	// Length of the edge after quantization
	m_mesh->newVariable<int,GMDS_EDGE>("length");
	// Mark with 1 the nodes forming a T-junction
	m_mesh->newVariable<int,GMDS_NODE>("T-junction");
	// Mark with 1 the nodes forming a T-junction (the previous variable is no longer necessary)
	//m_mesh->newVariable<int,GMDS_NODE>("is_a_T-junction");

	// Target mesh size
	m_mesh_size = AMeshSize;

	// Quantization graph
	m_quantization_graph = new QuantizationGraph();

	// Mesh with a singularity dipole added
	m_fixed_mesh = new Mesh(MeshModel(DIM3 | E | N | F |
	                                           E2N | N2E | F2N | N2F | F2E | E2F));
	// Mark with 1 the faces affected by the dipole
	m_mesh->newVariable<int,GMDS_FACE>("is_affected_by_dipole");

}

/*----------------------------------------------------------------------------*/
void QuantizationSolver::findTJunctions()
{
	auto tj = m_mesh->getVariable<int,GMDS_NODE>("is_a_T-junction");
	// Initialize at -1
	for (auto n_id:m_mesh->nodes())
		tj->set(n_id,-1);
	for (auto f_id:m_mesh->faces())
	{
		Face f = m_mesh->get<Face>(f_id);
		std::vector<Node> nodes = f.get<Node>();
		if (nodes.size() != 4)
		{
			// Then this face contains T-junctions
			int NbTJ = nodes.size() - 4;
			std::vector<Node> t_junctions;
			// Find the NbTJ flatest angles
			while (NbTJ > 0)
			{
				int min_pos;
				double min_value = 100;
				for (int i = 0; i < nodes.size(); i++)
				{
					Node n = nodes[i];
					// Find the two edges of f containing n
					Edge e1,e2;
					for (auto e:f.get<Edge>())
					{
						if (n.id() == e.get<Node>()[0].id() || n.id() == e.get<Node>()[1].id())
						{
							e1 = e;
							break;
						}
					}
					for (auto e:f.get<Edge>())
					{
						if (e.id() != e1.id() && (n.id() == e.get<Node>()[0].id() || n.id() == e.get<Node>()[1].id()))
						{
							e2 = e;
							break;
						}
					}
					double angle = oriented_angle(edge2vec(e1,n),edge2vec(e2,n));
					double dist;
					if (fabs(M_PI-angle) < fabs(M_PI+angle))
						dist = fabs(M_PI-angle);
					else
						dist = fabs(M_PI+angle);
					if (dist < min_value)
					{
						min_pos = i;
						min_value = dist;
					}
				}
				t_junctions.push_back(nodes[min_pos]);
				nodes.erase(nodes.begin() + min_pos);
				NbTJ -= 1;
			}
			for (auto n:t_junctions)
				tj->set(n.id(),f_id);
		}
	}
}

/*----------------------------------------------------------------------------*/
std::vector<std::vector<Edge>> QuantizationSolver::alignedEdgesGroups(Face &AF)
{
	auto tj = m_mesh->getVariable<int,GMDS_NODE>("is_a_T-junction");
	std::vector<Edge> ordered_edges = orientateEdges(AF);
	std::vector<std::vector<Edge>> edges_groups;
	std::vector<Edge> group;
	// group of the first element
	bool finished = false;
	Edge e0 = ordered_edges[0];
	Edge current_edge = ordered_edges[0];
	group.push_back(current_edge);
	ordered_edges.erase(ordered_edges.begin());
	Edge next_edge;
	while (!finished)
	{
		next_edge = ordered_edges[ordered_edges.size()-1];
		Node n = getCommonNode(current_edge,next_edge);
		if (tj->value(n.id()) == AF.id())
		{
			group.push_back(next_edge);
			ordered_edges.erase(ordered_edges.begin()+ordered_edges.size()-1);
			current_edge = next_edge;
		}
		else
			finished = true;
	}
	std::reverse(group.begin(),group.end());
	if (!ordered_edges.empty())
	{
		finished = false;
		current_edge = e0;
		while (!finished)
		{
			next_edge = ordered_edges[0];
			Node n = getCommonNode(current_edge,next_edge);
			if (tj->value(n.id()) == AF.id())
			{
				group.push_back(next_edge);
				ordered_edges.erase(ordered_edges.begin());
				current_edge = next_edge;
			}
			else
				finished = true;
		}
	}
	edges_groups.push_back(group);
	// Group of the other elements
	while (!ordered_edges.empty())
	{
		group.clear();
		current_edge = ordered_edges[0];
		group.push_back(current_edge);
		ordered_edges.erase(ordered_edges.begin());
		finished = false;
		while (!finished)
		{
			next_edge = ordered_edges[0];
			Node n = getCommonNode(current_edge,next_edge);
			if (tj->value(n.id()) == AF.id())
			{
				group.push_back(next_edge);
				ordered_edges.erase(ordered_edges.begin());
				current_edge = next_edge;
				if (ordered_edges.size() == 0)
					finished = true;
			}
			else
				finished = true;
		}
		edges_groups.push_back(group);
	}


	return edges_groups;
}

/*----------------------------------------------------------------------------*/
void
QuantizationSolver::buildHalfEdges()
{
	std::cout<<"> Building half edges"<<std::endl;
	auto e2he1 = m_mesh->getVariable<int,GMDS_EDGE>("edge2halfEdge1");
	auto e2he2 = m_mesh->getVariable<int,GMDS_EDGE>("edge2halfEdge2");
	for (auto e_id:m_mesh->edges())
	{
		e2he1->set(e_id,-1);
		e2he2->set(e_id,-1);
	}

	// Build the half edges
	int Id = 0;
	for (auto f_id:m_mesh->faces())
	{
		Face f = m_mesh->get<Face>(f_id);
		//std::vector<std::vector<Edge>> edges_groups = groupsOfAlignedEdges(f);
		std::vector<std::vector<Edge>> edges_groups = alignedEdgesGroups(f);
		if (edges_groups.size() != 4)
		{
			std::cout<<"The face "<<f_id<<" has "<<edges_groups.size()<<" groups of aligned edges"<<std::endl;
			// // test
			// for (auto edge:f.get<Edge>())
			// 	std::cout<<"test "<<edge.id()<<std::endl;


			// for (auto group:edges_groups)
			// {
			// 	std::cout<<"Group of edges :"<<std::endl;
			// 	for (auto edge:group)
			// 		std::cout<<"> "<<edge.id()<<std::endl;
			// }
			throw GMDSException("buildHalfEdges() : non conformal quads must have 4 half edges");
		}
		for (auto i = 0; i < 4; i++)
		{
			NonConformalHalfEdge new1 = NonConformalHalfEdge(Id+i,f,edges_groups[i]);
			for (auto e:edges_groups[i])
			{
				if (e2he1->value(e.id()) == -1)
					e2he1->set(e.id(),Id+i);
				else
					e2he2->set(e.id(),Id+i);
			}
			// Build next elements
			if (i < 3)
				new1.next(Id+i+1);
			else
				new1.next(Id);
			// Build first and last nodes
			if (i > 0)
				new1.firstNode(getCommonNode(edges_groups[i][0],edges_groups[i-1][edges_groups[i-1].size()-1]));
			else
				new1.firstNode(getCommonNode(edges_groups[i][0],edges_groups[i+3][edges_groups[i+3].size()-1]));
			if (i < 3)
				new1.endNode(getCommonNode(edges_groups[i][edges_groups[i].size()-1],edges_groups[i+1][0]));
			else
				new1.endNode(getCommonNode(edges_groups[i][edges_groups[i].size()-1],edges_groups[0][0]));
			m_half_edges.push_back(new1);
		}
		Id += 4;
	}

	// Set the opposite half edges
	for (auto half_edge:m_half_edges)
	{
		NonConformalHalfEdge he2 = half_edge;
		std::vector<int> opp;
		for (auto e:half_edge.edges())
		{
			int i1 = e2he1->value(e.id());
			int i2 = e2he2->value(e.id());
			if (i1 == half_edge.id())
			{
				bool alreadySeen = false;
				for (auto i:opp)
				{
					if (i == i2)
						alreadySeen = true;
				}
				if (i2 >= 0 && !alreadySeen)
					opp.push_back(i2);
			}
			else
			{
				bool alreadySeen = false;
				for (auto i:opp)
				{
					if (i == i1)
						alreadySeen = true;
				}
				if (i1 >= 0 && !alreadySeen)
					opp.push_back(i1);
			}
		}
		he2.opposite(opp);
		// Replace the old edge by the updated edge
		m_half_edges[half_edge.id()] = he2;
	}

	// test
//	for (auto half_edge:m_half_edges)
//	{
//		std::cout<<"Test ";
//		for (auto i:half_edge.opposite())
//			std::cout<<i<<" ";
//		std::cout<<std::endl;
//	}
}

/*----------------------------------------------------------------------------*/
Edge QuantizationSolver::getCommonEdge(int AHalfEdgeID1, int AHalfEdgeID2)
{
	NonConformalHalfEdge he1 = m_half_edges[AHalfEdgeID1];
	NonConformalHalfEdge he2 = m_half_edges[AHalfEdgeID2];
	Edge e;
	for (auto e1:he1.edges())
	{
		for (auto e2:he2.edges())
		{
			if (e1.id() == e2.id())
			{
				e = e1;
				return e;
			}
		}
	}
	throw GMDSException("getCommonEdge() : non conformal the two given half-edges are not opposite"); 
	return e;
}

/*----------------------------------------------------------------------------*/
void
QuantizationSolver::buildQuantizationGraphNodes()
{
	std::cout<<"> Building quantization graph nodes"<<std::endl;
	for (int i = 0; i < m_half_edges.size(); i++)
	{
		NonConformalHalfEdge he = m_half_edges[i];
		double len = he.firstNode().point().distance(he.endNode().point());
		Node n = m_quantization_graph->newNode();
		m_quantization_graph->setGeometricalLength(n.id(),len);
	}
}

/*----------------------------------------------------------------------------*/
void
QuantizationSolver::buildConnectedComponent(int AHalfEdgeID)
{
	std::cout<<"> Building connected component of half edge "<<AHalfEdgeID<<std::endl;
	auto e2qge = m_mesh->getVariable<int,GMDS_EDGE>("edge2qgrapheEdge");
	std::queue<int> front;
	front.push(AHalfEdgeID);
	int current;
	while (!front.empty())
	{
		current = front.front();
		front.pop();
		if (m_quantization_graph->alreadyVisited(current) == 0)
		{
			m_quantization_graph->markAsVisited(current);
			NonConformalHalfEdge e = m_half_edges[current];
			// Build edges linking to opp
			for (auto opp:e.opposite())
			{
				if (m_quantization_graph->alreadyVisited(opp) == 0)
				{
					Edge newGraphEdge = m_quantization_graph->newEdge(current,opp);
					Edge commonEdge = getCommonEdge(current,opp);
					e2qge->set(commonEdge.id(),newGraphEdge.id());
					if (m_quantization_graph->alreadyPushed(oppositeInQuad(opp)) == 0)
					{
						front.push(oppositeInQuad(opp));
						m_quantization_graph->markAsPushed(oppositeInQuad(opp));
					}
				}
			}
			
			int nxt = oppositeInQuad(current);
			if (m_quantization_graph->alreadyVisited(nxt) == 0)
			{
				// Build edge linking to next(next)
				NonConformalHalfEdge e2 = m_half_edges[nxt];
				m_quantization_graph->markAsVisited(nxt);
				Edge newGraphEdge = m_quantization_graph->newEdge(nxt,current);
				m_quantization_graph->markAsInQuad(newGraphEdge.id());
				// Build edges linking next(next) to opp(next(next))
				for (auto opp:e2.opposite())
				{
					if (m_quantization_graph->alreadyVisited(opp) == 0)
					{
						newGraphEdge = m_quantization_graph->newEdge(opp,nxt);
						Edge commonEdge = getCommonEdge(opp,nxt);
						e2qge->set(commonEdge.id(),newGraphEdge.id());
						if (m_quantization_graph->alreadyPushed(opp) == 0)
						{
							front.push(opp);
							m_quantization_graph->markAsPushed(opp);
						}
					}
				}
			}
		}
	}

}

/*----------------------------------------------------------------------------*/
int
QuantizationSolver::oppositeInQuad(int AHalfEdgeID)
{
	int opp = m_half_edges[m_half_edges[AHalfEdgeID].next()].next();
	return opp;
}

/*----------------------------------------------------------------------------*/
void
QuantizationSolver::buildQuantizationGraph()
{
	std::cout<<" "<<std::endl;
	std::cout<<"========== Building the quantization graph =========="<<std::endl;
	//findTJunctions();
	buildHalfEdges();
	buildQuantizationGraphNodes();
	for (int i = 0; i < m_half_edges.size(); i++)
	{
		if (m_quantization_graph->alreadyVisited(i) == 0)
			buildConnectedComponent(i);
	}
	std::cout<<"====================================================="<<std::endl;
	std::cout<<" "<<std::endl;
	// test
	//m_quantization_graph->display();
}

/*----------------------------------------------------------------------------*/
QuantizationGraph*
QuantizationSolver::getQuantizationGraph()
{
	return m_quantization_graph;
}

/*----------------------------------------------------------------------------*/
void QuantizationSolver::setHalfEdgesLength()
{
	std::vector<int> v(m_half_edges.size());
	for (int i = 0; i < m_half_edges.size(); i++)
		v[i] = m_quantization_graph->quantizationSolutionValue(i);
	m_half_edges_lengths = v;
}

/*----------------------------------------------------------------------------*/
void QuantizationSolver::setEdgesLength()
{
	auto l = m_mesh->getVariable<int,GMDS_EDGE>("length");
	auto e2qge = m_mesh->getVariable<int,GMDS_EDGE>("edge2qgrapheEdge");
	for (auto e_id:m_mesh->edges())
		l->set(e_id,m_quantization_graph->fluxValue(e2qge->value(e_id)));
	for (auto he:m_half_edges)
	{
		if (he.opposite().empty())
		{
			Edge e = he.edges()[0];
			l->set(e.id(),m_half_edges_lengths[he.id()]);
		}
	}
}

/*----------------------------------------------------------------------------*/
std::vector<std::vector<TCellID>> QuantizationSolver::sides()
{
	auto e2he1 = m_mesh->getVariable<int,GMDS_EDGE>("edge2halfEdge1");
	std::vector<std::vector<TCellID>> sides;
	auto alreadySeen = m_mesh->newMark<Edge>();
	for (auto e_id:m_mesh->edges())
	{
		if (!m_mesh->isMarked<Edge>(e_id,alreadySeen))
		{
			Edge e = m_mesh->get<Edge>(e_id);
			if (e.get<Face>().size() == 1)
			{
				// Then it is a boundary edge, we build a new side
				std::vector<TCellID> newSide;
				std::queue<Edge> front;
				front.push(e);
				m_mesh->mark<Edge>(e_id,alreadySeen);
				while (!front.empty())
				{
					Edge e1 = front.front();
					front.pop();
					newSide.push_back(e2he1->value(e1.id()));
					for (auto n:e1.get<Node>())
					{
						for (auto e2:n.get<Edge>())
						{
							if (e2.get<Face>().size() == 1)
							{
								if (!m_mesh->isMarked<Edge>(e2.id(),alreadySeen))
								{
									double alpha = oriented_angle(edge2vec(e1,n),edge2vec(e2,n));
									if (fabs(alpha-M_PI) < 0.1 || fabs(alpha+M_PI) < 0.1)
									{
										front.push(e2);
										m_mesh->mark<Edge>(e2.id(),alreadySeen);
									}
								}
							}
						}
					}
				}
				sides.push_back(newSide);
			}
		}
	}
	return sides;
}

/*----------------------------------------------------------------------------*/
void QuantizationSolver::markTJunctions()
{
	auto TJ = m_mesh->getVariable<int,GMDS_NODE>("T-junction");
	auto e2he1 = m_mesh->getVariable<int,GMDS_EDGE>("edge2halfEdge1");
	auto e2he2 = m_mesh->getVariable<int,GMDS_EDGE>("edge2halfEdge2");
	for (auto n_id:m_mesh->nodes())
	{
		Node n = m_mesh->get<Node>(n_id);
		bool isATJ = false;
		for (auto e:n.get<Edge>())
		{
			NonConformalHalfEdge he1 = m_half_edges[e2he1->value(e.id())];
			if (n.id() != he1.firstNode().id() && n.id() != he1.endNode().id())
			{
				isATJ = true;
				break;
			}
			if (e2he2->value(e.id()) >= 0)
			{
				he1 = m_half_edges[e2he2->value(e.id())];
				if (n.id() != he1.firstNode().id() && n.id() != he1.endNode().id())
				{
					isATJ = true;
					break;
				}
			}
		}
		if (isATJ)
			TJ->set(n_id,1);
	}
}

/*----------------------------------------------------------------------------*/
bool QuantizationSolver::isOnBoundary(Node AN)
{
	bool boo = false;
	for (auto e:AN.get<Edge>())
	{
		if (e.get<Face>().size() == 1)
		{
			boo = true;
			break;
		}
	}
	return boo;
}

/*----------------------------------------------------------------------------*/
std::vector<TCellID> QuantizationSolver::nonZeroHalfEdges()
{
	auto e2he1 = m_mesh->getVariable<int,GMDS_EDGE>("edge2halfEdge1");
	auto e2he2 = m_mesh->getVariable<int,GMDS_EDGE>("edge2halfEdge2");
	auto TJ = m_mesh->getVariable<int,GMDS_NODE>("T-junction");
	std::vector<int> alreadyAdded(m_half_edges.size());
	for (int i = 0; i < m_half_edges.size(); i++)
		alreadyAdded[i] = 0;
	std::vector<TCellID> nzhe;
	for (auto n_id:m_mesh->nodes())
	{
		Node n = m_mesh->get<Node>(n_id);
		if (!isOnBoundary(n))
		{
			if (TJ->value(n_id) == 0)
			{
				std::vector<Edge> adj_edges = n.get<Edge>();
				if (adj_edges.size() != 4)
				{
					// Then it is an interior singularity
					for (auto e:adj_edges)
					{
						int he = e2he1->value(e.id());
						if (alreadyAdded[he] == 0)
						{
							nzhe.push_back(he);
							alreadyAdded[he] = 1;
						}
						he = e2he2->value(e.id());
						if (he >= 0)
						{
							if (alreadyAdded[he] == 0)
							{
								nzhe.push_back(he);
								alreadyAdded[he] = 1;
							}
						}
						
					}
				}
			}
		}
		else
		{
			std::vector<Edge> adj_edges = n.get<Edge>();
			if (adj_edges.size() != 3)
			{
				// Then it is a boundary singularity
				for (auto e:adj_edges)
				{
					if (e.get<Face>().size() == 2)
					{
						int he = e2he1->value(e.id());
						if (alreadyAdded[he] == 0)
						{
							nzhe.push_back(he);
							alreadyAdded[he] = 1;
						}
						he = e2he2->value(e.id());
						if (he >= 0)
						{
							if (alreadyAdded[he] == 0)
							{
								nzhe.push_back(he);
								alreadyAdded[he] = 1;
							}
						}
					}
				}
			}
		}
	}
	return nzhe;
}

/*----------------------------------------------------------------------------*/
std::vector<std::vector<Node>> QuantizationSolver::groupsOfConfundedNodes()
{
	std::cout<<"> Building groups of confunded nodes"<<std::endl;
	auto groupId = m_mesh->getVariable<int,GMDS_NODE>("groupId");
	auto e2he1 = m_mesh->getVariable<int,GMDS_EDGE>("edge2halfEdge1");
	auto e2he2 = m_mesh->getVariable<int,GMDS_EDGE>("edge2halfEdge2");
	std::vector<std::vector<Node>> groups;
	auto visited = m_mesh->newMark<Node>();
	int Id = 0;
	for (auto n_id:m_mesh->nodes())
	{
		if (!m_mesh->isMarked<Node>(n_id,visited))
		{
			// Then we create a new group
			std::vector<Node> newGroup;
			std::queue<Node> toVisit;
			Node n0 = m_mesh->get<Node>(n_id);
			toVisit.push(n0);
			m_mesh->mark<Node>(n_id,visited);
			while (!toVisit.empty())
			{
				Node n = toVisit.front();
				toVisit.pop();
				newGroup.push_back(n);
				groupId->set(n.id(),Id);
				for (auto e:n.get<Edge>())
				{
					int he_id = e2he1->value(e.id());
					if (m_half_edges_lengths[he_id] == 0)
					{
						NonConformalHalfEdge he = m_half_edges[he_id];
						for (auto n1:he.getOrderedNodes())
						{
							if (!m_mesh->isMarked<Node>(n1.id(),visited))
							{
								toVisit.push(n1);
								m_mesh->mark<Node>(n1.id(),visited);
							}
						}
					}
					he_id = e2he2->value(e.id());
					if (he_id >= 0)
					{
						if (m_half_edges_lengths[he_id] == 0)
						{
							NonConformalHalfEdge he = m_half_edges[he_id];
							for (auto n1:he.getOrderedNodes())
							{
								if (!m_mesh->isMarked<Node>(n1.id(),visited))
								{
									toVisit.push(n1);
									m_mesh->mark<Node>(n1.id(),visited);
								}
							}
						}
					}
				}
			}
			groups.push_back(newGroup);
			Id += 1;
		}
	}
	return groups;
}

/*----------------------------------------------------------------------------*/
bool QuantizationSolver::isABoundaryCorner(Node AN)
{
	if (!isOnBoundary(AN))
		return false;
	else
	{
		// Find the two boundary edges adjacent at AN
		Edge e1,e2;
		for (auto e:AN.get<Edge>())
		{
			if (e.get<Face>().size() == 1)
				{
					e1 = e;
					break;
				}
		}
		for (auto e:AN.get<Edge>())
		{
			if (e.get<Face>().size() == 1 && e.id() != e1.id())
				{
					e2 = e;
					break;
				}
		}
		double alpha = oriented_angle(edge2vec(e1,AN),edge2vec(e2,AN));
		if (fabs(alpha-M_PI) > 0.1 && fabs(alpha+M_PI) > 0.1)
			return true;
		else
			return false;
	}
}

/*----------------------------------------------------------------------------*/
void QuantizationSolver::addOldNodesOnFixedMesh()
{
	std::cout<<"> Adding existing nodes on the fixed mesh"<<std::endl;
	for (auto n_id:m_mesh->nodes())
	{
		Node n = m_mesh->get<Node>(n_id);
		m_fixed_mesh->newNode(n.point());
	}
}

/*----------------------------------------------------------------------------*/
void QuantizationSolver::writeFixedMesh(std::basic_string<char> AFileName)
{
	std::cout<<"> Writing the fixed mesh"<<std::endl;
	IGMeshIOService ioService(m_fixed_mesh);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(AFileName);
}

/*----------------------------------------------------------------------------*/
void QuantizationSolver::markAffectedFaces(int AId1, int AId2)
{
	auto var = m_mesh->getVariable<int,GMDS_FACE>("is_affected_by_dipole");
	NonConformalHalfEdge he1 = m_half_edges[AId1];
	NonConformalHalfEdge he2 = m_half_edges[he1.next()];
	NonConformalHalfEdge he3 = m_half_edges[AId2];
	NonConformalHalfEdge he4 = m_half_edges[he3.next()];
	Face f0 = he1.face();
	var->set(f0.id(),1);
	for (auto e:he2.edges())
	{
		for (auto f:e.get<Face>())
		{
			if (f.id() != f0.id())
				var->set(f.id(),1);
		}
	}
	for (auto e:he4.edges())
	{
		for (auto f:e.get<Face>())
		{
			if (f.id() != f0.id())
				var->set(f.id(),1);
		}
	}
}

/*----------------------------------------------------------------------------*/
void QuantizationSolver::addOldFacesOnFixedMesh()
{
	std::cout<<"> Adding existing faces on the fixed mesh"<<std::endl;
	auto var = m_mesh->getVariable<int,GMDS_FACE>("is_affected_by_dipole");
	for (auto f_id:m_mesh->faces())
	{
		if (var->value(f_id) == 0)
		{
			Face f = m_mesh->get<Face>(f_id);
			std::vector<TCellID> nodesIds;
			for (auto n:f.get<Node>())
				nodesIds.push_back(n.id());
			m_fixed_mesh->newFace(nodesIds);
		}
	}
}

/*----------------------------------------------------------------------------*/
void QuantizationSolver::buildOppFacesInFixedMesh(int AHalfEdgeId, Node AN1, Node AN2)
{
	NonConformalHalfEdge he = m_half_edges[AHalfEdgeId];
	Face f = he.face();

	// Build the faces containing AN1 and AN2
	int id1 = -1;
	// Find the face where AN1 is
	for (auto edge:he.edges())
	{
		if (AN1.id() == edge.get<Node>()[0].id() || AN1.id() == edge.get<Node>()[1].id())
			break;
		if (isOnSegment(AN1.point(),edge.get<Node>()[0].point(),edge.get<Node>()[1].point()))
		{
			for (auto f2:edge.get<Face>())
			{
				if (f2.id() != f.id())
					{
						id1 = f2.id();
						break;
					}
			}
		}
	}

	int id2 = -1;
	// Find the face where AN2 is
	for (auto edge:he.edges())
	{
		if (AN2.id() == edge.get<Node>()[0].id() || AN2.id() == edge.get<Node>()[1].id())
			break;
		if (isOnSegment(AN2.point(),edge.get<Node>()[0].point(),edge.get<Node>()[1].point()))
		{
			for (auto f2:edge.get<Face>())
			{
				if (f2.id() != f.id())
					{
						id2 = f2.id();
						break;
					}
			}
		}
	}

	// Now we build the eventual faces containing the new points
	if (id1 >= 0 && id1 == id2)
	{
		Face f2 = m_mesh->get<Face>(id1);
		std::vector<TCellID> new_face = insertPoint(AN1,f2.get<Node>());
		std::vector<Node> nodes;
		for (auto id:new_face)
			nodes.push_back(m_fixed_mesh->get<Node>(id));
		new_face = insertPoint(AN2,nodes);
		m_fixed_mesh->newFace(new_face);
	}
	else
	{
		if (id1 >= 0)
		{
			Face f2 = m_mesh->get<Face>(id1);
			std::vector<TCellID> new_face = insertPoint(AN1,f2.get<Node>());
			m_fixed_mesh->newFace(new_face);
		}
		if (id2 >= 0)
		{
			Face f2 = m_mesh->get<Face>(id2);
			std::vector<TCellID> new_face = insertPoint(AN2,f2.get<Node>());
			m_fixed_mesh->newFace(new_face);
		}
	}

	// Build the other faces
	for (auto edge:he.edges())
	{
		for (auto f2:edge.get<Face>())
		{
			if (f2.id() != f.id() && f2.id() != id1 &&f2.id() != id2)
			{
				std::vector<TCellID> nodes;
				for (auto n:f2.get<Face>())
					nodes.push_back(n.id());
				m_fixed_mesh->newFace(nodes);
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
void QuantizationSolver::buildFixedMesh(int AId1, int AId2)
{
	std::cout<<" "<<std::endl;
	std::cout<<"=== Building the fixed mesh ==="<<std::endl;

	addOldNodesOnFixedMesh();
	markAffectedFaces(AId1, AId2);
	addOldFacesOnFixedMesh();

	NonConformalHalfEdge he1 = m_half_edges[AId1];
	NonConformalHalfEdge he2 = m_half_edges[he1.next()];
	NonConformalHalfEdge he3 = m_half_edges[AId2];
	NonConformalHalfEdge he4 = m_half_edges[he3.next()];

	// Drawing the singularity dipole

	// Corner nodes
	Node c1 = m_fixed_mesh->get<Node>(he1.endNode().id());
	Node c2 = m_fixed_mesh->get<Node>(he2.endNode().id());
	Node c3 = m_fixed_mesh->get<Node>(he3.endNode().id());
	Node c4 = m_fixed_mesh->get<Node>(he4.endNode().id());

	// Directions
	math::Point d1 = c1.point()+(-1.)*c4.point();
	math::Point d2 = c2.point()+(-1.)*c1.point();
	math::Point d3 = c3.point()+(-1.)*c2.point();
	math::Point d4 = c4.point()+(-1.)*c3.point();

	// Nodes on edges
	// Find their distance to corner nodes
	double h1,h2,h3,h4;
	std::vector<Node> nodes_he2 = he2.getOrderedNodes();
	std::vector<Node> nodes_he4 = he4.getOrderedNodes();
	if (nodes_he2.size() == 2)
	{
		h1 = 1./5.;
		h2 = 4./5.;
	}
	if (nodes_he2.size() >= 3)
	{
		if (1./5.+1./10. < nodes_he2[0].point().distance(nodes_he2[1].point())/vec(d2).norm())
			h1 = 1./5.;
		else
			h1 = 0.5*nodes_he2[0].point().distance(nodes_he2[1].point())/vec(d2).norm();
		if (4./5.-1./10. > nodes_he2[0].point().distance(nodes_he2[nodes_he2.size()-2].point())/vec(d2).norm())
			h2 = 4./5.;
		else
			h2 = 1.-0.5*nodes_he2[nodes_he2.size()-1].point().distance(nodes_he2[nodes_he2.size()-2].point())/vec(d2).norm();
	}
	if (nodes_he4.size() == 2)
	{
		h3 = 1./5.;
		h4 = 4./5.;
	}
	if (nodes_he4.size() >= 3)
	{
		if (1./5.+1./10. < nodes_he4[0].point().distance(nodes_he4[1].point())/vec(d4).norm())
			h3 = 1./5.;
		else
			h3 = 0.5*nodes_he4[0].point().distance(nodes_he4[1].point())/vec(d4).norm();
		if (4./5.-1./10. > nodes_he4[0].point().distance(nodes_he4[nodes_he4.size()-2].point())/vec(d4).norm())
			h4 = 4./5.;
		else
			h4 = 1.-0.5*nodes_he4[nodes_he4.size()-1].point().distance(nodes_he4[nodes_he4.size()-2].point())/vec(d4).norm();
	}
	Node b1 = m_fixed_mesh->newNode(c1.point()+h1*d2);
	Node b2 = m_fixed_mesh->newNode(c1.point()+h2*d2);
	Node b3 = m_fixed_mesh->newNode(c3.point()+h3*d4);
	Node b4 = m_fixed_mesh->newNode(c3.point()+h4*d4);
	
	// Build the faces opposite to he2 and he4
	buildOppFacesInFixedMesh(he2.id(),b1,b2);
	buildOppFacesInFixedMesh(he4.id(),b3,b4);
	
	// Directions
	math::Point dg = b4.point()+(-1.)*b1.point();
	math::Point dd = b3.point()+(-1.)*b2.point();

	// Nodes on inner edges
	Node s1 = m_fixed_mesh->newNode(b1.point()+(1./6.)*dg);
	Node s4 = m_fixed_mesh->newNode(b2.point()+(1./6.)*dd);
	Node s5 = m_fixed_mesh->newNode(b2.point()+(5./6.)*dd);
	Node s8 = m_fixed_mesh->newNode(b1.point()+(5./6.)*dg);
	Node s2 = m_fixed_mesh->newNode(s1.point()+(1./3.)*s4.point()+(-1./3.)*s1.point());
	Node s3 = m_fixed_mesh->newNode(s1.point()+(2./3.)*s4.point()+(-2./3.)*s1.point());
	Node s6 = m_fixed_mesh->newNode(s8.point()+(2./3.)*s5.point()+(-2./3.)*s8.point());
	Node s7 = m_fixed_mesh->newNode(s8.point()+(1./3.)*s5.point()+(-1./3.)*s8.point());

	// Directions
	math::Point d5 = s8.point()+(-1.)*s1.point();
	math::Point d7 = s7.point()+(-1.)*s2.point();
	math::Point d6 = (1./2.)*(d5+d7);
	math::Point d8 = s6.point()+(-1.)*s3.point();
	math::Point d9 = s5.point()+(-1.)*s4.point();

	// Inner nodes
	Node i1 = m_fixed_mesh->newNode(s1.point()+(1./4.)*d5);
	Node i2 = m_fixed_mesh->newNode(s2.point()+(1./4.)*d7);
	Node i3 = m_fixed_mesh->newNode(s1.point()+(2./4.)*d5);
	Node i4 = m_fixed_mesh->newNode((1./2.)*s1.point()+(1./2.)*s2.point()+(1./2.)*d6);
	Node i5 = m_fixed_mesh->newNode(s3.point()+(1./2.)*d8);
	Node i6 = m_fixed_mesh->newNode(s4.point()+(1./2.)*d9);
	Node i7 = m_fixed_mesh->newNode(s1.point()+(3./4.)*d5);
	Node i8 = m_fixed_mesh->newNode(s2.point()+(3./4.)*d7);

	// Create new faces
	std::vector<TCellID> new_face;
	for (auto node:he1.getOrderedNodes())
		new_face.push_back(node.id());
	new_face.push_back(b1.id());
	new_face.push_back(s1.id());
	new_face.push_back(i1.id());
	new_face.push_back(i3.id());
	new_face.push_back(i7.id());
	new_face.push_back(s8.id());
	new_face.push_back(b4.id());
	m_fixed_mesh->newFace(new_face);
	
	new_face.clear();
	for (auto node:he3.getOrderedNodes())
		new_face.push_back(node.id());
	new_face.push_back(b3.id());
	new_face.push_back(s5.id());
	new_face.push_back(i6.id());
	new_face.push_back(s4.id());
	new_face.push_back(b2.id());
	m_fixed_mesh->newFace(new_face);
	
	new_face.clear();
	new_face.push_back(b3.id());
	for (auto node:he4.getOrderedNodes())
		{
			if (node.id() != c3.id() && node.id() != c4.id())
				new_face.push_back(node.id());
		}
	new_face.push_back(b4.id());
	new_face.push_back(s8.id());
	new_face.push_back(s7.id());
	new_face.push_back(s6.id());
	new_face.push_back(s5.id());
	m_fixed_mesh->newFace(new_face);

	new_face.clear();
	new_face.push_back(b1.id());
	for (auto node:he2.getOrderedNodes())
		{
			if (node.id() != c1.id() && node.id() != c2.id())
				new_face.push_back(node.id());
		}
	new_face.push_back(b2.id());
	new_face.push_back(s4.id());
	new_face.push_back(s3.id());
	new_face.push_back(s2.id());
	new_face.push_back(s1.id());
	m_fixed_mesh->newFace(new_face);
	
	m_fixed_mesh->newFace({i7.id(),i8.id(),s7.id(),s8.id()});
	m_fixed_mesh->newFace({i3.id(),i4.id(),i8.id(),i7.id()});
	m_fixed_mesh->newFace({i1.id(),i2.id(),i4.id(),i3.id()});
	m_fixed_mesh->newFace({s1.id(),s2.id(),i2.id(),i1.id()});
	m_fixed_mesh->newFace({i8.id(),i5.id(),s6.id(),s7.id()});
	m_fixed_mesh->newFace({i2.id(),i5.id(),i8.id(),i4.id()});
	m_fixed_mesh->newFace({s2.id(),s3.id(),i5.id(),i2.id()});
	m_fixed_mesh->newFace({i5.id(),i6.id(),s5.id(),s6.id()});
	m_fixed_mesh->newFace({s3.id(),s4.id(),i6.id(),i5.id()});

	std::cout<<"==============================="<<std::endl;
	std::cout<<" "<<std::endl;
}

/*----------------------------------------------------------------------------*/
Mesh QuantizationSolver::getFixedMesh()
{
	return *m_fixed_mesh;
}

/*----------------------------------------------------------------------------*/
void QuantizationSolver::setFixedMeshConnectivity()
{
	MeshDoctor doc(m_fixed_mesh);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
}

/*----------------------------------------------------------------------------*/
void QuantizationSolver::forceZeroEdgesToZeroLength()
{
	std::cout<<"> Imposing degenerated half edges to be of length 0"<<std::endl;
	for (int id = 0; id<m_half_edges.size(); id++)
	{
		NonConformalHalfEdge he = m_half_edges[id];
		double l = he.firstNode().point().distance(he.endNode().point());
		if (l < 1e-6)
		{
			m_quantization_graph->markAsZero(id);
			std::cout<<"Imposing a half edge to 0"<<std::endl;
		}
	}
}

/*----------------------------------------------------------------------------*/
std::vector<std::vector<Node>> QuantizationSolver::buildCompleteSolution()
{
	std::cout<<" "<<std::endl;
    std::cout<<"=== Building complete quantization solution ==="<<std::endl;

	buildQuantizationGraph();
	// Set conditions on sides
	std::vector<std::vector<TCellID>> sidesSet = sides();
	m_quantization_graph->nonZeroGroups(sidesSet);
	// Set condition to separate singularities
	markTJunctions();
	std::vector<TCellID> non_zeros = nonZeroHalfEdges();
	m_quantization_graph->nonZeroVerticies(non_zeros);
	m_quantization_graph->updateConnectivity();
	//forceZeroEdgesToZeroLength();
	m_quantization_graph->buildCompleteSolution();
	//m_quantization_graph->roughlyRepairSolution();
	//m_quantization_graph->display(); // Display the graph
	//m_quantization_graph->displaySolution(); // Display the quantization solution
	std::vector<std::vector<Node>> problematicCouples;// = m_quantization_graph->checkSolutionValidity();
	m_quantization_graph->improveSolution(m_mesh_size);

	std::cout<<"===================================="<<std::endl;
    std::cout<<" "<<std::endl;

	return problematicCouples;
}
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/