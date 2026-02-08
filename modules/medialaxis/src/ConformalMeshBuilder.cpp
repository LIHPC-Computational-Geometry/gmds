#include "gmds/medialaxis/ConformalMeshBuilder.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
ConformalMeshBuilder::ConformalMeshBuilder(gmds::Mesh &AMesh, std::vector<NonConformalHalfEdge> AHalfEdges, std::vector<int> ALengths)
{
	// Non conformal quad mesh
	m_mesh = &AMesh;
	// Nodes of the quantized mesh corresponding to each non-conformal face
	m_mesh->newVariable<std::vector<std::vector<int>>,GMDS_FACE>("face2quantizedNodes");
	
	// Corresponding non-conformal half edges
	m_half_edges = AHalfEdges;
	// Half edges length. WARNING : all lengths must be > 0
	m_half_edges_lengths = ALengths;

	// Quantized mesh
	m_quantized_mesh = new Mesh(MeshModel(DIM3 | E | N | F |
	                                           E2N | N2E | F2N | N2F | F2E | E2F));
}

/*----------------------------------------------------------------------------*/
void ConformalMeshBuilder::writeQuantizedMesh(std::basic_string<char> AFileName)
{
	std::cout<<"> Writing the quantized mesh"<<std::endl;
	IGMeshIOService ioService(m_quantized_mesh);
	VTKWriter vtkWriter(&ioService);
	vtkWriter.setCellOptions(N| E| F);
	vtkWriter.setDataOptions(N| E| F);
	vtkWriter.write(AFileName);
}

/*----------------------------------------------------------------------------*/
Mesh ConformalMeshBuilder::getQuantizedMesh()
{
	return *m_quantized_mesh;
}

/*----------------------------------------------------------------------------*/
void ConformalMeshBuilder::setQuantizedMeshConnectivity()
{
	std::cout<<"> Setting quantized mesh connectivity"<<std::endl;
	MeshDoctor doc(m_quantized_mesh);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();
}

/*----------------------------------------------------------------------------*/
void ConformalMeshBuilder::buildQuantizedMeshNodesOnEdges()
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
void ConformalMeshBuilder::buildQuantizedMeshInternalNodes()
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
void ConformalMeshBuilder::buildQuantizedMeshFaces()
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
std::vector<std::vector<int>> ConformalMeshBuilder::halfEdgesGroup(int AID, std::vector<int> &AV)
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
std::vector<std::vector<std::vector<int>>> ConformalMeshBuilder::halfEdgesGroups()
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
void ConformalMeshBuilder::execute()
{
	std::cout<<" "<<std::endl;
	std::cout<<"=== Building the conformal mesh ==="<<std::endl;

	buildQuantizedMeshNodesOnEdges();
	buildQuantizedMeshInternalNodes();
	buildQuantizedMeshFaces();
	setQuantizedMeshConnectivity();

	std::cout<<"==================================="<<std::endl;
	std::cout<<" "<<std::endl;
}
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/