/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/mctsblock/Graph.h>
#include <iostream>
#include <set>
/*----------------------------------------------------------------------------*/
TEST(GraphTestSuite, ordered_graph)
{
	std::set<gmds::TCellID> ids;
	ids.insert(0);
	ids.insert(1);
	ids.insert(2);
	ids.insert(3);
	ids.insert(4);
	ids.insert(5);
	ids.insert(6);
	ids.insert(7);
	ids.insert(8);
	gmds::mctsblock::Graph g(ids);
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
/*--------------------------------------------------*/
   TEST(GraphTestSuite, unordered_graph)
{
	std::set<gmds::TCellID> ids;
	ids.insert(0);
	ids.insert(10);
	ids.insert(20);
	ids.insert(30);
	ids.insert(40);
	ids.insert(50);
	ids.insert(60);
	ids.insert(70);
	ids.insert(80);
	gmds::mctsblock::Graph g(ids);
	g.addEdge(0, 10, 4);
	g.addEdge(0, 70, 8);
	g.addEdge( 10, 20, 8);
	g.addEdge( 10, 70, 11);
	g.addEdge( 20, 30, 7);
	g.addEdge( 20, 80, 2);
	g.addEdge( 20, 50, 4);
	g.addEdge( 30, 40, 9);
	g.addEdge( 30, 50, 14);
	g.addEdge( 40, 50, 10);
	g.addEdge( 50, 60, 2);
	g.addEdge( 60, 70, 1);
	g.addEdge( 60, 80, 6);
	g.addEdge( 70, 80, 7);

	g.computeDijkstra(0);

	auto paths = g.getShortestPath();
	auto weights = g.getShortestPathWeights();

	ASSERT_EQ(paths[10].size(),2);
	ASSERT_EQ(weights[10],4);

	ASSERT_EQ(paths[20].size(),3);
	ASSERT_EQ(weights[20],12);

	ASSERT_EQ(paths[30].size(),4);
	ASSERT_EQ(weights[30],19);

	ASSERT_EQ(paths[40].size(),5);
	ASSERT_EQ(weights[40],21);

	ASSERT_EQ(paths[50].size(),4);
	ASSERT_EQ(weights[50],11);

	ASSERT_EQ(paths[60].size(),3);
	ASSERT_EQ(weights[60],9);

	ASSERT_EQ(paths[70].size(),2);
	ASSERT_EQ(weights[70],8);

	ASSERT_EQ(paths[80].size(),4);
	ASSERT_EQ(weights[80],14);

}