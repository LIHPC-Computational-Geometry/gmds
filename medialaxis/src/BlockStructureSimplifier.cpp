#include "gmds/medialaxis/BlockStructureSimplifier.h"
#include <queue>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
BlockStructureSimplifier::BlockStructureSimplifier(gmds::Mesh &AMesh)
{
	// Non conformal quad mesh
	m_mesh = &AMesh;
	// Correspondence edge/half edge
	m_mesh->newVariable<int,GMDS_EDGE>("edge2halfEdge1");
	m_mesh->newVariable<int,GMDS_EDGE>("edge2halfEdge2");
	// Mark with 1 the nodes forming a T-junction
	m_mesh->newVariable<int,GMDS_NODE>("T-junction");
	// Id of the corresponding confunded nodes group
	m_mesh->newVariable<int,GMDS_NODE>("groupId");
	// In the mesh simplifying process, gain of the chosing of the node
	m_mesh->newVariable<double,GMDS_NODE>("chosing_gain");

	// Quantization graph
	m_quantization_graph = new QuantizationGraph();

	// Simplified mesh(zero half edges removed)
	m_simplified_mesh = new Mesh(MeshModel(DIM3 | E | N | F |
	                                           E2N | N2E | F2N | N2F | F2E | E2F));

}

/*----------------------------------------------------------------------------*/
void
BlockStructureSimplifier::buildHalfEdges()
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
		std::vector<std::vector<Edge>> edges_groups = groupsOfAlignedEdges(f);
		if (edges_groups.size() != 4)
			throw GMDSException("buildHalfEdges() : non conformal quads must have 4 non conformal edges");
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
void
BlockStructureSimplifier::buildQuantizationGraphNodes()
{
	std::cout<<"> Building quantization graph nodes"<<std::endl;
	for (int i = 0; i < m_half_edges.size(); i++)
		m_quantization_graph->newNode();
}

/*----------------------------------------------------------------------------*/
void
BlockStructureSimplifier::buildConnectedComponent(int AHalfEdgeID)
{
	std::cout<<"> Building connected component of half edge "<<AHalfEdgeID<<std::endl;
	std::queue<int> front;
	front.push(AHalfEdgeID);
	int current;
	while (!front.empty())
	{
		current = front.front();
		front.pop();
		m_quantization_graph->markAsVisited(current);
		NonConformalHalfEdge e = m_half_edges[current];
		// Build forwards edges
		for (auto opp:e.opposite())
		{
			if (m_quantization_graph->alreadyVisited(opp) == 0)
			{
				m_quantization_graph->newEdge(current,opp);
				if (m_quantization_graph->alreadyPushed(oppositeInQuad(opp)) == 0)
				{
					front.push(oppositeInQuad(opp));
					m_quantization_graph->markAsPushed(oppositeInQuad(opp));
				}
			}
		}
		// Build backward edges
		int nxt = oppositeInQuad(current);
		NonConformalHalfEdge e2 = m_half_edges[nxt];
		m_quantization_graph->markAsVisited(nxt);
		Edge newGraphEdge = m_quantization_graph->newEdge(nxt,current);
		m_quantization_graph->markAsInQuad(newGraphEdge.id());
		for (auto opp:e2.opposite())
		{
			if (m_quantization_graph->alreadyVisited(opp) == 0)
			{
				m_quantization_graph->newEdge(opp,nxt);
				if (m_quantization_graph->alreadyPushed(opp) == 0)
				{
					front.push(opp);
					m_quantization_graph->markAsPushed(opp);
				}
			}
		}
	}

}

/*----------------------------------------------------------------------------*/
int
BlockStructureSimplifier::oppositeInQuad(int AHalfEdgeID)
{
	int opp = m_half_edges[m_half_edges[AHalfEdgeID].next()].next();
	return opp;
}

/*----------------------------------------------------------------------------*/
void
BlockStructureSimplifier::buildQuantizationGraph()
{
	std::cout<<" "<<std::endl;
	std::cout<<"========== Building the quantization graph =========="<<std::endl;
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
BlockStructureSimplifier::getQuantizationGraph()
{
	return m_quantization_graph;
}

/*----------------------------------------------------------------------------*/
void BlockStructureSimplifier::setHalfEdgesLength()
{
	std::vector<int> v(m_half_edges.size());
	for (int i = 0; i < m_half_edges.size(); i++)
		v[i] = m_quantization_graph->quantizationSolutionValue(i);
	m_half_edges_lengths = v;
}

/*----------------------------------------------------------------------------*/
std::vector<std::vector<TCellID>> BlockStructureSimplifier::sides()
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
void BlockStructureSimplifier::markTJunctions()
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
bool BlockStructureSimplifier::isOnBoundary(Node AN)
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
std::vector<TCellID> BlockStructureSimplifier::nonZeroHalfEdges()
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
std::vector<std::vector<Node>> BlockStructureSimplifier::groupsOfConfundedNodes()
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
bool BlockStructureSimplifier::isABoundaryCorner(Node AN)
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
void BlockStructureSimplifier::computeChosingGains()
{
	auto cost = m_mesh->getVariable<double,GMDS_NODE>("chosing_gain");
	for (auto n_id:m_mesh->nodes())
	{
		Node n = m_mesh->get<Node>(n_id);
		if (isABoundaryCorner(n))
			cost->set(n_id,16.);
		else
		{
			if (isOnBoundary(n))
				cost->set(n_id,15.);
			else
			{
				if (n.get<Edge>().size() != 4)
					cost->set(n_id,14.);
				else
				{
					std::vector<Edge> adj_edges = n.get<Edge>();
					std::vector<Edge> sorted_edges = sortEdges(n,adj_edges);
					math::Vector u1 = edge2vec(adj_edges[0],n);
					math::Vector u2 = (-1.)*edge2vec(adj_edges[2],n);
					math::Vector v1 = edge2vec(adj_edges[1],n);
					math::Vector v2 = (-1.)*edge2vec(adj_edges[3],n);
					double c = sqrt(fabs(oriented_angle(u1,u2))*fabs(oriented_angle(u1,u2))+fabs(oriented_angle(v1,v2))*fabs(oriented_angle(v1,v2)));
					cost->set(n_id,c);
				}
			}
		}
	}
}

/*----------------------------------------------------------------------------*/
Node BlockStructureSimplifier::representative(std::vector<Node> AV)
{
	auto gain = m_mesh->getVariable<double,GMDS_NODE>("chosing_gain");
	Node rep;
	double max_gain = -1.;
	for (auto n:AV)
	{
		if (gain->value(n.id()) > max_gain)
		{
			rep = n;
			max_gain = gain->value(n.id());
		}
	}
	// if (AV.size() == 1)
	// 	return AV[0];
	// for (auto n:AV)
	// {
	// 	if (isABoundaryCorner(n))
	// 		return n;
	// }
	// for (auto n:AV)
	// {
	// 	if (isOnBoundary(n))
	// 		return n;
	// }
	// return AV[0];
	return rep;
}

/*----------------------------------------------------------------------------*/
void BlockStructureSimplifier::buildSimplifiedMesh()
{
	std::cout<<"=== Building simplified mesh ==="<<std::endl;
	auto groupId = m_mesh->getVariable<int,GMDS_NODE>("groupId");
	// build the nodes
	std::vector<std::vector<Node>> groups = groupsOfConfundedNodes();
	std::vector<TCellID> group2newNode;
	for (auto group:groups)
	{
		Node newNode = m_simplified_mesh->newNode(representative(group).point());
		group2newNode.push_back(newNode.id());
	}
	// Build the faces
	for (auto f_id:m_mesh->faces())
	{
		// Check if it has non zero area
		if (m_half_edges_lengths[4*f_id] > 0 && m_half_edges_lengths[4*f_id+1] > 0)
		{
			std::vector<TCellID> nodes;
			Face f = m_mesh->get<Face>(f_id);
			std::vector<int> seenGroups;
			for (auto n:f.get<Node>())
			{
				int id = groupId->value(n.id());
				bool seen = false;
				for (auto i:seenGroups)
				{
					if (id == i)
					{
						seen = true;
						break;
					}
				}
				if (!seen)
				{
					nodes.push_back(group2newNode[id]);
					seenGroups.push_back(id);
				}
			}
			m_simplified_mesh->newFace(nodes);
		}
	}
	std::cout<<"============================"<<std::endl;
}

/*----------------------------------------------------------------------------*/
void BlockStructureSimplifier::setSimplifiedMeshConnectivity()
{
	MeshDoctor doc(m_simplified_mesh);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
}

/*----------------------------------------------------------------------------*/
Mesh BlockStructureSimplifier::getSimplifiedMesh()
{
	return *m_simplified_mesh;
}

/*----------------------------------------------------------------------------*/
void BlockStructureSimplifier::writeSimplifiedMesh(std::basic_string<char> AFileName)
{
	std::cout<<"> Writing the simplified mesh"<<std::endl;
	IGMeshIOService ioService(m_simplified_mesh);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(AFileName);
}

/*----------------------------------------------------------------------------*/
void BlockStructureSimplifier::execute()
{
    std::cout<<" "<<std::endl;
    std::cout<<"=== Building the simplified mesh ==="<<std::endl;
    
    buildQuantizationGraph();
	// Set conditions on sides
	std::vector<std::vector<TCellID>> sides2 = sides();
	m_quantization_graph->nonZeroGroups(sides2);
	// Set condition to separate singularities
	markTJunctions();
	std::vector<TCellID> non_zeros2 = nonZeroHalfEdges();
	m_quantization_graph->nonZeroVerticies(non_zeros2);
	m_quantization_graph->updateConnectivity();
	m_quantization_graph->buildMinimalSolution();
	m_quantization_graph->roughlyRepairSolution();
	//m_quantization_graph->displaySolution();
	setHalfEdgesLength();
	computeChosingGains();
	buildSimplifiedMesh();

    std::cout<<"===================================="<<std::endl;
    std::cout<<" "<<std::endl;
}
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/