#ifndef GMDS_ORIENTED_GRAPH_TESTSUITE_H
#define GMDS_ORIENTED_GRAPH_TESTSUITE_H

#include "gtest/gtest.h"
#include <gmds/utils/OrientedGraph.h>
#include <unit_test_config.h>

using namespace gmds;

TEST(OrientedGraphTestSuite, GraphNodeAddOutEdge)
{
	GraphNode node1(1);
	GraphNode node2(2);
	GraphEdge edge(1, &node1, &node2);

	ASSERT_TRUE(node1.addOutEdge(&edge));
	ASSERT_EQ(node1.outEdges().size(), 1);
	ASSERT_EQ(node1.outEdges()[0]->id(), edge.id());
}

TEST(OrientedGraphTestSuite, GraphNodeAddInEdge)
{
	GraphNode node1(1);
	GraphNode node2(2);
	GraphEdge edge(1, &node1, &node2);

	ASSERT_TRUE(node2.addInEdge(&edge));
	ASSERT_EQ(node2.inEdges().size(), 1);
	ASSERT_EQ(node2.inEdges()[0]->id(), edge.id());
}

TEST(OrientedGraphTestSuite, GraphEdgeGetEdgesStartingFrom)
{
	GraphNode node1(1);
	GraphNode node2(2);
	GraphEdge edge1(1, &node1, &node2);
	GraphEdge edge2(2, &node1, &node2);

	ASSERT_EQ(edge1.getEdgesStartingFrom(&node1).size(), 1);
	ASSERT_EQ(edge1.getEdgesStartingFrom(&node2).size(), 0);
	ASSERT_EQ(edge2.getEdgesStartingFrom(&node1).size(), 1);
	ASSERT_EQ(edge2.getEdgesStartingFrom(&node2).size(), 0);
}

TEST(OrientedGraphTestSuite, OrientedGraphAddEdge)
{
	OrientedGraph graph(3);
	ASSERT_TRUE(graph.addEdge(1, 0, 1));
	ASSERT_EQ(graph.nbEdges(), 1);
}

TEST(OrientedGraphTestSuite, OrientedGraphUpdateNodes)
{
	OrientedGraph graph(3);
	graph.addEdge(1, 0, 1);
	graph.updateNodes();

	ASSERT_EQ(graph.node(0)->outEdges().size(), 1);
	ASSERT_EQ(graph.node(1)->inEdges().size(), 1);
}

TEST(OrientedGraphTestSuite, OrientedGraphNode)
{
	OrientedGraph graph(3);
	GraphNode *node = graph.node(1);
	ASSERT_EQ(node->id(), 1);
}

TEST(OrientedGraphTestSuite, OrientedGraphEdge)
{
	OrientedGraph graph(3);
	graph.addEdge(1, 0, 1);
	GraphEdge *edge = graph.edge(0);
	ASSERT_EQ(edge->id(), 1);
}

#endif     // ORIENTED_GRAPH_TEST_H
