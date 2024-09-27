//
// Created by chenyt on 26/09/24.
//

#include "gmds/medialaxis/QuantizationGraphBuilder.h"
#include <queue>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
QuantizationGraphBuilder::QuantizationGraphBuilder(gmds::Mesh &AMesh)
{
	// Non conformal quad mesh
	m_mesh = &AMesh;
	// Correspondence edge/half edge
	m_mesh->newVariable<int,GMDS_EDGE>("edge2halfEdge1");
	m_mesh->newVariable<int,GMDS_EDGE>("edge2halfEdge2");

	// Quantization graph
	m_quantization_graph = new QuantizationGraph();
}

/*----------------------------------------------------------------------------*/
void QuantizationGraphBuilder::buildHalfEdges()
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
			if (i < 3)
				new1.next(Id+i+1);
			else
				new1.next(Id);
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
void QuantizationGraphBuilder::buildQuantizationGraphNodes()
{
	std::cout<<"> Building quantization graph nodes"<<std::endl;
	for (int i = 0; i < m_half_edges.size(); i++)
		m_quantization_graph->newNode();
}

/*----------------------------------------------------------------------------*/
void QuantizationGraphBuilder::buildConnectedComponent(int AHalfEdgeID)
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
int QuantizationGraphBuilder::oppositeInQuad(int AHalfEdgeID)
{
	int opp = m_half_edges[m_half_edges[AHalfEdgeID].next()].next();
	return opp;
}

/*----------------------------------------------------------------------------*/
void QuantizationGraphBuilder::execute()
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
QuantizationGraph* QuantizationGraphBuilder::getQuantizationGraph()
{
	return m_quantization_graph;
}
/*----------------------------------------------------------------------------*/
}  // end namespace gmds
/*----------------------------------------------------------------------------*/