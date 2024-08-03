/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/mctsblock/Graph.h>
#include <iostream>
/*----------------------------------------------------------------------------*/
TEST(GraphTestSuite, box)
{
	gmds::mctsblock::Graph g(9);
	g.addEdge(0, 1, 4);
	g.addEdge(0, 7, 8);
	g.addEdge( 1, 2, 8);
	g.addEdge( 1, 7, 11);
	g.addEdge( 2, 3, 7);
	g.addEdge( 2, 8, 2);
	g.addEdge( 2, 5, 4);
	g.addEdge( 3, 4, 9);
	g.addEdge( 3, 5, 14);
	g.addEdge( 4, 5, 10);
	g.addEdge( 5, 6, 2);
	g.addEdge( 6, 7, 1);
	g.addEdge( 6, 8, 6);
	g.addEdge( 7, 8, 7);

	g.computeDijkstra(0);

	auto paths = g.getShortestPath();
	auto weights = g.getShortestPathWeights();

	ASSERT_EQ(paths[1].size(),2);
	ASSERT_EQ(weights[1],4);

	ASSERT_EQ(paths[2].size(),3);
	ASSERT_EQ(weights[2],12);

	ASSERT_EQ(paths[3].size(),4);
	ASSERT_EQ(weights[3],19);

	ASSERT_EQ(paths[4].size(),5);
	ASSERT_EQ(weights[4],21);

	ASSERT_EQ(paths[5].size(),4);
	ASSERT_EQ(weights[5],11);

	ASSERT_EQ(paths[6].size(),3);
	ASSERT_EQ(weights[6],9);

	ASSERT_EQ(paths[7].size(),2);
	ASSERT_EQ(weights[7],8);

	ASSERT_EQ(paths[8].size(),4);
	ASSERT_EQ(weights[8],14);

}