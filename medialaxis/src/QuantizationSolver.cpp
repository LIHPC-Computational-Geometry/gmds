//
// Created by chenyt on 26/09/24.
//

#include "gmds/medialaxis/QuantizationSolver.h"
#include <queue>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
QuantizationSolver::QuantizationSolver(gmds::Mesh &AMesh)
{
	// Non conformal quad mesh
	m_mesh = &AMesh;
	// Correspondence edge/half edge
	m_mesh->newVariable<int,GMDS_EDGE>("edge2halfEdge1");
	m_mesh->newVariable<int,GMDS_EDGE>("edge2halfEdge2");
	// Nodes of the quantized mesh corresponding to each non-conformal face
	m_mesh->newVariable<std::vector<std::vector<int>>,GMDS_FACE>("face2quantizedNodes");

	// Quantization graph
	m_quantization_graph = new QuantizationGraph();

	// Quantized mesh
	m_quantized_mesh = new Mesh(MeshModel(DIM3 | E | N | F |
	                                           E2N | N2E | F2N | N2F | F2E | E2F));

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
QuantizationSolver::buildQuantizationGraphNodes()
{
	std::cout<<"> Building quantization graph nodes"<<std::endl;
	for (int i = 0; i < m_half_edges.size(); i++)
		m_quantization_graph->newNode();
}

/*----------------------------------------------------------------------------*/
void
QuantizationSolver::buildConnectedComponent(int AHalfEdgeID)
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
		m_quantization_graph->newEdge(nxt,current);
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
std::vector<std::vector<int>> QuantizationSolver::halfEdgesGroup(int AID, std::vector<int> &AV)
{
	std::vector<std::vector<int>> groups;
	std::vector<int> group1;
	std::vector<int> group2;
	NonConformalHalfEdge e = m_half_edges[AID];
	bool over = false;
	while(!over)
	{
		group1.push_back(e.id());
		AV[e.id()] = 1;
		for (auto opp:e.opposite())
		{
			if (AV[opp] == 0)
			{
				group2.push_back(opp);
				AV[opp] = 1;
			}
		}
		if (e.opposite().empty())
			over = true;
		else
		{
			NonConformalHalfEdge e2 = m_half_edges[e.opposite()[e.opposite().size()-1]];
			if (e2.firstNode().id() == e.endNode().id())
				over = true;
			else
				e = m_half_edges[m_half_edges[m_half_edges[e.next()].opposite()[0]].next()];
		}
	}
	groups.push_back(group1);
	groups.push_back(group2);
	return groups;
}

/*----------------------------------------------------------------------------*/
std::vector<std::vector<std::vector<int>>> QuantizationSolver::halfEdgesGroups()
{
	// Vector to mark the already seen half edges, initialized at 0
	std::vector<int> V(m_half_edges.size());
	for (int i = 0; i < m_half_edges.size(); i++)
		V[i] = 0;

	// Build the groups
	std::vector<std::vector<std::vector<int>>> groups;
	std::vector<std::vector<int>> group;
	for (int i = 0; i < m_half_edges.size(); i++)
	{
		if (V[i] == 0)
		{
			NonConformalHalfEdge e1 = m_half_edges[i];
			if (e1.opposite().empty())
			{
				group = halfEdgesGroup(i,V);
				groups.push_back(group);
			}
			else
			{
				NonConformalHalfEdge e2 = m_half_edges[e1.opposite()[0]];
				if (e1.firstNode().id() == e2.endNode().id())
				{
					group = halfEdgesGroup(i,V);
					groups.push_back(group);
				}
			}
		}
	}

//	// Test : display the groups
//	for (auto g:groups)
//	{
//		std::cout<<"Group 1 : ";
//		for (auto i:g[0])
//			std::cout<<i<<" ";
//		std::cout<<"Group 2 : ";
//		for (auto i:g[1])
//			std::cout<<i<<" ";
//		std::cout<<std::endl;
//	}

	return groups;
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
void QuantizationSolver::writeQuantizedMesh(std::basic_string<char> AFileName)
{
	std::cout<<"> Writing the quantized mesh"<<std::endl;
	IGMeshIOService ioService(m_quantized_mesh);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(AFileName);
}

/*----------------------------------------------------------------------------*/
void QuantizationSolver::buildQuantizedMeshNodesOnEdges()
{
	std::cout<<"> Building nodes of the quantized mesh that lay on non-conformal edges"<<std::endl;
	// Initialize half_edges_to_nodes
	std::vector<std::vector<TCellID>> v(m_half_edges.size());
	m_half_edges_to_new_nodes = v;
	// Add the nodes of the non-conformal mesh
	for (auto n_id:m_mesh->nodes())
	{
		Node n = m_mesh->get<Node>(n_id);
		m_quantized_mesh->newNode(n.point());
	}
	// Build the new nodes
	std::vector<std::vector<std::vector<int>>> hegroups = halfEdgesGroups();
	for (auto g:hegroups)
	{
		// Get the groups of half edges on each side
		std::vector<int> g1 = g[0];
		std::vector<int> g2 = g[1];
		// Compute the total length of the group
		int L = 0;
		for (auto i:g1)
			L += m_half_edges_lengths[i];
		// If the half edges is on the boundary, ie g2 is empty
		if (g2.empty())
		{
			NonConformalHalfEdge e = m_half_edges[g1[0]];
			std::vector<TCellID> newNodes;
			newNodes.push_back(e.firstNode().id());
			math::Point E = e.endNode().point()+(-1.)*e.firstNode().point();
			math::Point P;
			for (int l = 1; l < L; l++)
			{
				P = e.firstNode().point()+(double (l)/double(L))*E;
				Node n = m_quantized_mesh->newNode(P);
				newNodes.push_back(n.id());
			}
			newNodes.push_back(e.endNode().id());
			m_half_edges_to_new_nodes[e.id()] = newNodes;
		}
		else // g2 is not empty
		{
			NonConformalHalfEdge e1 = m_half_edges[g1[0]];
			NonConformalHalfEdge e2 = m_half_edges[g2[0]];
			int i1 = 0;
			int i2 = 0;
			int L1 = 0;
			int L2 = 0;
			std::vector<TCellID> newNodes1;
			std::vector<TCellID> newNodes2;
			newNodes1.push_back(e1.firstNode().id());
			newNodes2.push_back(e2.endNode().id());
			for (int l = 1; l < L; l++)
			{
				if (l < L1+m_half_edges_lengths[e1.id()] && l < L2+m_half_edges_lengths[e2.id()])
				{
					// Create the new node
					math::Point dir1 = e1.endNode().point()+(-1.)*m_quantized_mesh->get<Node>(newNodes1[newNodes1.size()-1]).point();
					math::Point dir2 = e2.firstNode().point()+(-1.)*m_quantized_mesh->get<Node>(newNodes2[newNodes2.size()-1]).point();
					math::Point P;
					if (vec(dir1).norm()/double(m_half_edges_lengths[e1.id()]) < vec(dir2).norm()/double(m_half_edges_lengths[e2.id()]))
						P = m_quantized_mesh->get<Node>(newNodes1[newNodes1.size()-1]).point()+(1./double(m_half_edges_lengths[e1.id()]))*dir1;
					else
						P = m_quantized_mesh->get<Node>(newNodes2[newNodes2.size()-1]).point()+(1./double(m_half_edges_lengths[e2.id()]))*dir2;
					Node n = m_quantized_mesh->newNode(P);
					newNodes1.push_back(n.id());
					newNodes2.push_back(n.id());
				}
				if (l == L1+m_half_edges_lengths[e1.id()] && l < L2+m_half_edges_lengths[e2.id()])
				{
					newNodes1.push_back(e1.endNode().id());
					newNodes2.push_back(e1.endNode().id());
					L1 = L1 + m_half_edges_lengths[e1.id()];
					e1 = m_half_edges[g1[i1+1]];
					m_half_edges_to_new_nodes[g1[i1]] = newNodes1;
					i1 += 1;
					newNodes1.clear();
					newNodes1.push_back(e1.firstNode().id());
				}
				if (l < L1+m_half_edges_lengths[e1.id()] && l == L2+m_half_edges_lengths[e2.id()])
				{
					newNodes1.push_back(e2.firstNode().id());
					newNodes2.push_back(e2.firstNode().id());
					L2 = L2 + m_half_edges_lengths[e2.id()];
					e2 = m_half_edges[g2[i2+1]];
					std::reverse(newNodes2.begin(),newNodes2.end());
					m_half_edges_to_new_nodes[g2[i2]] = newNodes2;
					i2 += 1;
					newNodes2.clear();
					newNodes2.push_back(e2.endNode().id());
				}
			}
			newNodes1.push_back(m_half_edges[g1[g1.size()-1]].endNode().id());
			newNodes2.push_back(m_half_edges[g2[g2.size()-1]].firstNode().id());
			std::reverse(newNodes2.begin(),newNodes2.end());
			m_half_edges_to_new_nodes[g1[i1]] = newNodes1;
			m_half_edges_to_new_nodes[g2[i2]] = newNodes2;
		}
	}
}

/*----------------------------------------------------------------------------*/
void QuantizationSolver::buildQuantizedMeshInternalNodes()
{
	std::cout<<"> Building nodes of the quantized mesh that are internal to non-conformal faces"<<std::endl;
	auto f2qn = m_mesh->getVariable<std::vector<std::vector<int>>,GMDS_FACE>("face2quantizedNodes");
	for (auto f_id:m_mesh->faces())
	{
		// Get the half edges of the face
		NonConformalHalfEdge e1 = m_half_edges[4*f_id];
		NonConformalHalfEdge e2 = m_half_edges[4*f_id+1];
		NonConformalHalfEdge e3 = m_half_edges[4*f_id+2];
		NonConformalHalfEdge e4 = m_half_edges[4*f_id+3];
		std::vector<std::vector<int>> nodesIDs(m_half_edges_lengths[e2.id()]+1);
		std::vector<int> row(m_half_edges_lengths[e1.id()]+1);
		// first row
		for (int i = 0; i < m_half_edges_lengths[e1.id()]+1; i++)
			row[i] = m_half_edges_to_new_nodes[e1.id()][i];
		nodesIDs[0] = row;
		// Last row
		for (int i = 0; i < m_half_edges_lengths[e3.id()]+1; i++)
			row[i] = m_half_edges_to_new_nodes[e3.id()][m_half_edges_lengths[e3.id()]-i];
		nodesIDs[m_half_edges_lengths[e2.id()]] = row;
		for (int j = 1; j < m_half_edges_lengths[e2.id()]; j++)
		{
			row[0] = m_half_edges_to_new_nodes[e4.id()][m_half_edges_lengths[e4.id()]-j];
			row[m_half_edges_lengths[e1.id()]] = m_half_edges_to_new_nodes[e2.id()][j];
			for (int i = 1; i < m_half_edges_lengths[e1.id()]; i++)
			{
				math::Point P1 = m_quantized_mesh->get<Node>(m_half_edges_to_new_nodes[e1.id()][i]).point();
				math::Point P2 = m_quantized_mesh->get<Node>(m_half_edges_to_new_nodes[e3.id()][m_half_edges_lengths[e3.id()]-i]).point();
				math::Point P = P1 + (double(j)/double(m_half_edges_lengths[e2.id()]))*(P2+(-1.)*P1);
				Node n = m_quantized_mesh->newNode(P);
				row[i] = n.id();
			}
			nodesIDs[j] = row;
		}
		f2qn->set(f_id,nodesIDs);
	}
}

/*----------------------------------------------------------------------------*/
void QuantizationSolver::buildQuantizedMeshFaces()
{
	std::cout<<"> Building faces of the quantized mesh"<<std::endl;
	auto f2qn = m_mesh->getVariable<std::vector<std::vector<int>>,GMDS_FACE>("face2quantizedNodes");
	for (auto f_id:m_mesh->faces())
	{
		std::vector<std::vector<int>> nodes = f2qn->value(f_id);
		int I = nodes.size();
		int J = nodes[0].size();
		TCellID i1,i2,i3,i4;
		for (int i = 0; i < I-1; i++)
		{
			for (int j = 0; j < J-1; j++)
			{
				i1 = nodes[i][j];
				i2 = nodes[i+1][j];
				i3 = nodes[i+1][j+1];
				i4 = nodes[i][j+1];
				m_quantized_mesh->newFace({i1,i2,i3,i4});
			}
		}
	}
}
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/