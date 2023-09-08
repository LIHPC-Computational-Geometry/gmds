/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/math/Point.h>
/*----------------------------------------------------------------------------*/
#include <gmds/ant_colony/Env.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;

/*----------------------------------------------------------------------------*/
TEST(EnvTestClass, test_constructor)
{
	// Given a matrix 
	int exist_tens[3][3][3] = {{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}},
                                {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}},
                                {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
	// Flatten 
	int *exist_flatten = &exist_tens[0][0][0];
	gmds::Env env(3, 3, 3, exist_flatten);
	env.writeVTK("env.vtk");
	std::vector<TCellID> front = {45};
}
/*----------------------------------------------------------------------------*/
TEST(EnvTestClass, test_getCommonEdge)
{
	int exist_tens[3][3][3] = {{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}},
                                {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}},
                                {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
	// Flatten 
	int *exist_flatten = &exist_tens[0][0][0];
	gmds::Env env(3, 3, 3, exist_flatten);
	std::vector<std::pair<TCellID,TCellID>> front = {std::pair<TCellID,TCellID>(42, 40)};
	TCellID edg = env.getCommonEdge(front[0].first,front[0].second);
	ASSERT_EQ( edg, 63 );
	std::cout << edg << std::endl;
}


/*----------------------------------------------------------------------------*/
TEST(EnvTestClass, test_execute_bis)
{
	int exist_tens[3][3][3] = {{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}},
                                {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}},
                                {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
	// Flatten 
	int *exist_flatten = &exist_tens[0][0][0];
	gmds::Env env(3, 3, 3, exist_flatten);
	std::vector<std::pair<TCellID,TCellID>> front = {std::pair<TCellID,TCellID>(42, 40)};
	env.execute_bis(front);
}

/*----------------------------------------------------------------------------*/

TEST(SimplexMeshTestClass, test_point_in_two_tetra)
{
}
