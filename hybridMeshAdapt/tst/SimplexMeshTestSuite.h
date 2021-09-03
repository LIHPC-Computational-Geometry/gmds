/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::hybrid;
using namespace gmds::math;
using namespace hybrid;
using namespace hybrid::simplicesCell;
using namespace hybrid::simplicesNode;
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test_of_one_tetra)
{
  //exemple 1
  gmds::BitVector node_ids(4);
  node_ids.fillAll();
  std::vector<gmds::math::Point> coords;
  coords.push_back(Point(0.0,0.0,0.0));
  coords.push_back(Point(1.0,0.0,0.0));
  coords.push_back(Point(0.0,0.0,1.0));
  coords.push_back(Point(0.0,1.0,0.0));

  std::vector<TSimplexID> base({0,0,0,0});

  gmds::BitVector tet_ids(1);
  tet_ids.fillAll();
  std::vector<TSimplexID* > tet_nodes;
  TSimplexID* T1 = new TSimplexID[4]{0,1,2,3};

  tet_nodes.push_back(T1);

  std::vector<TSimplexID* > tet_adj;
  TSimplexID erorId = std::numeric_limits<int>::min();
  TSimplexID* Ta1 = new TSimplexID[4]{erorId,erorId,erorId,erorId};
  tet_adj.push_back(Ta1);

  gmds::BitVector tri_ids(0);
  std::vector<TSimplexID* > tri_nodes;
  std::vector<TSimplexID* > tri_adj  ;

  SimplexMesh mesh(node_ids,coords,base,
                   tet_ids,tet_nodes,tet_adj,
                   tri_ids,tri_nodes,tri_adj);

  SimplicesNode simplicesNode0(&mesh, 0);
  /*operateur de d'affectation semantique + operateur de deplacement semantic*/
  std::vector<TInt>&& ballOf = simplicesNode0.ballOf();

  ASSERT_EQ(ballOf.size(),1);
  ASSERT_EQ(ballOf[0],0);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[0], 0);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[0], 0);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3))[0], 0);

}

/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test_of_two_tetra)
{
  //exemple 2
  gmds::BitVector node_ids(5);
  node_ids.fillAll();
  std::vector<gmds::math::Point> coords;
  coords.push_back(Point(0.0,0.0,0.0));
  coords.push_back(Point(1.0,0.0,0.0));
  coords.push_back(Point(0.0,0.0,1.0));
  coords.push_back(Point(0.0,1.0,0.0));
  coords.push_back(Point(0.0,-1.0,0.0));

  std::vector<TSimplexID> base({0,0,0,0,1});

  gmds::BitVector tet_ids(2);
  tet_ids.fillAll();
  std::vector<TSimplexID* > tet_nodes;
  TSimplexID* T1 = new TSimplexID[4]{0,1,2,3};
  TSimplexID* T2 = new TSimplexID[4]{0,1,2,4};

  tet_nodes.push_back(T1);
  tet_nodes.push_back(T2);

  std::vector<TSimplexID* > tet_adj;
  TSimplexID erorId = std::numeric_limits<int>::min();
  TSimplexID* Ta1 = new TSimplexID[4]{erorId, erorId, erorId, 1};
  TSimplexID* Ta2 = new TSimplexID[4]{erorId, erorId, erorId, 0};

  tet_adj.push_back(Ta1);
  tet_adj.push_back(Ta2);

  gmds::BitVector tri_ids;
  std::vector<TSimplexID* > tri_nodes;
  std::vector<TSimplexID* > tri_adj  ;

  SimplexMesh mesh(node_ids,coords,base,
                   tet_ids,tet_nodes,tet_adj,
                   tri_ids,tri_nodes,tri_adj);

  SimplicesNode simplicesNode0(&mesh, 0);
  std::vector<TInt> ballOf = simplicesNode0.ballOf();

  //ball
  ASSERT_EQ(ballOf.size(),2);
  ASSERT_EQ(ballOf[0],0);
  ASSERT_EQ(ballOf[1],1);

  //shell
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1)).size(), 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[0], 0);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[1], 1);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2)).size(), 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[0], 0);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[1], 1);


  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3))[0],     0);

}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test_of_three_tetra)
{

  //exemple 1
  gmds::BitVector node_ids(6);
  node_ids.fillAll();
  std::vector<gmds::math::Point> coords;
  coords.push_back(Point(0.0,0.0,0.0));
  coords.push_back(Point(1.0,0.0,0.0));
  coords.push_back(Point(0.0,0.0,1.0));
  coords.push_back(Point(0.0,1.0,0.0));
  coords.push_back(Point(0.0,-1.0,0.0));
  coords.push_back(Point(0.0,0.0,-1.0));

  std::vector<TSimplexID> base({0,0,0,0,1,2});

  gmds::BitVector tet_ids(3);
  tet_ids.fillAll();
  std::vector<TSimplexID* > tet_nodes;
  TSimplexID* T1 = new TSimplexID[4]{0,1,2,3};
  TSimplexID* T2 = new TSimplexID[4]{0,1,2,4};
  TSimplexID* T3 = new TSimplexID[4]{0,1,4,5};

  tet_nodes.push_back(T1);
  tet_nodes.push_back(T2);
  tet_nodes.push_back(T3);

  std::vector<TSimplexID* > tet_adj;
  TSimplexID erorId = std::numeric_limits<int>::min();
  TSimplexID* Ta1 = new TSimplexID[4]{ erorId, erorId, erorId,1};
  TSimplexID* Ta2 = new TSimplexID[4]{ erorId, erorId,      2,0};
  TSimplexID* Ta3 = new TSimplexID[4]{ erorId, erorId, erorId,1};

  tet_adj.push_back(Ta1);
  tet_adj.push_back(Ta2);
  tet_adj.push_back(Ta3);

  gmds::BitVector tri_ids;
  std::vector<TSimplexID* > tri_nodes;
  std::vector<TSimplexID* > tri_adj  ;

  SimplexMesh mesh(node_ids,coords,base,
                   tet_ids,tet_nodes,tet_adj,
                   tri_ids,tri_nodes,tri_adj);

  SimplicesNode simplicesNode0(&mesh, 0);
  std::vector<TInt> ballOf = simplicesNode0.ballOf();

  ASSERT_EQ(ballOf.size(),3);
  ASSERT_EQ(ballOf[0],0);
  ASSERT_EQ(ballOf[1],1);
  ASSERT_EQ(ballOf[2],2);

  //ball
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1)).size(), 3);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[0], 0);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[1], 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[2], 2);

  //shell
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4)).size(), 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4))[0], 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4))[1], 2);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3))[0], 0);


  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2)).size(), 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[0], 0);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[1], 1);


}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test_ball_of_four_tetra)
{

  //exemple 1
  gmds::BitVector node_ids(6);
  node_ids.fillAll();
  std::vector<gmds::math::Point> coords;
  coords.push_back(Point(0.0,0.0,0.0));
  coords.push_back(Point(1.0,0.0,0.0));
  coords.push_back(Point(0.0,0.0,1.0));
  coords.push_back(Point(0.0,1.0,0.0));
  coords.push_back(Point(0.0,-1.0,0.0));
  coords.push_back(Point(0.0,0.0,-1.0));

  std::vector<TSimplexID> base({0,0,0,3,1,2});

  gmds::BitVector tet_ids(4);
  tet_ids.fillAll();
  std::vector<TSimplexID* > tet_nodes;
  TSimplexID* T1 = new TSimplexID[4]{0,1,2,3};
  TSimplexID* T2 = new TSimplexID[4]{0,1,2,4};
  TSimplexID* T3 = new TSimplexID[4]{0,1,4,5};
  TSimplexID* T4 = new TSimplexID[4]{0,1,3,5};

  tet_nodes.push_back(T1);
  tet_nodes.push_back(T2);
  tet_nodes.push_back(T3);
  tet_nodes.push_back(T4);

  std::vector<TSimplexID* > tet_adj;
  TSimplexID erorId = std::numeric_limits<int>::min();
  TSimplexID* Ta1 = new TSimplexID[4]{erorId,erorId,     3,1};
  TSimplexID* Ta2 = new TSimplexID[4]{erorId,erorId,     2,0};
  TSimplexID* Ta3 = new TSimplexID[4]{erorId,erorId,     3,1};
  TSimplexID* Ta4 = new TSimplexID[4]{erorId,erorId,     2,0};

  tet_adj.push_back(Ta1);
  tet_adj.push_back(Ta2);
  tet_adj.push_back(Ta3);
  tet_adj.push_back(Ta4);

  gmds::BitVector tri_ids;
  std::vector<TSimplexID* > tri_nodes;
  std::vector<TSimplexID* > tri_adj  ;

  SimplexMesh mesh(node_ids,coords,base,
                   tet_ids,tet_nodes,tet_adj,
                   tri_ids,tri_nodes,tri_adj);



  SimplicesNode simplicesNode0(&mesh, 0);
  std::vector<TInt> ballOf = simplicesNode0.ballOf();

  ASSERT_EQ(ballOf.size(),4);
  ASSERT_EQ(ballOf[0],0);
  ASSERT_EQ(ballOf[1],1);
  ASSERT_EQ(ballOf[2],2);
  ASSERT_EQ(ballOf[3],3);

  //ball
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1)).size(), 4);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[0], 0);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[1], 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[2], 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[3], 3);

  //shell
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3)).size(), 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3))[0], 0);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3))[1], 3);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4)).size(), 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4))[0], 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4))[1], 2);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2)).size(), 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[0], 0);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[1], 1);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 5)).size(), 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 5))[0], 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 5))[1], 3);

}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test_ball_of_eight_tetra)
{

  //exemple 1
  gmds::BitVector node_ids(7);
  node_ids.fillAll();
  std::vector<gmds::math::Point> coords;
  coords.push_back(Point(0.0,0.0,0.0));
  coords.push_back(Point(1.0,0.0,0.0));
  coords.push_back(Point(0.0,0.0,1.0));
  coords.push_back(Point(0.0,1.0,0.0));
  coords.push_back(Point(0.0,-1.0,0.0));
  coords.push_back(Point(0.0,0.0,-1.0));
  coords.push_back(Point(-1.0,0.0,0.0));

  std::vector<TSimplexID> base({7,0,1,4,2,3,6});

  gmds::BitVector tet_ids(8);
  tet_ids.fillAll();
  std::vector<TSimplexID* > tet_nodes;
  TSimplexID* T1 = new TSimplexID[4]{0,1,2,3};
  TSimplexID* T2 = new TSimplexID[4]{0,1,2,4};
  TSimplexID* T3 = new TSimplexID[4]{0,1,4,5};
  TSimplexID* T4 = new TSimplexID[4]{0,1,3,5};

  TSimplexID* T5 = new TSimplexID[4]{0,5,6,3};
  TSimplexID* T6 = new TSimplexID[4]{0,5,6,4};
  TSimplexID* T7 = new TSimplexID[4]{0,4,2,6};
  TSimplexID* T8 = new TSimplexID[4]{0,3,2,6};

  tet_nodes.push_back(T1);
  tet_nodes.push_back(T2);
  tet_nodes.push_back(T3);
  tet_nodes.push_back(T4);
  tet_nodes.push_back(T5);
  tet_nodes.push_back(T6);
  tet_nodes.push_back(T7);
  tet_nodes.push_back(T8);

  std::vector<TSimplexID* > tet_adj;
  TSimplexID erorId = std::numeric_limits<int>::min();
  TSimplexID* Ta1 = new TSimplexID[4]{erorId,     7,     3,1};
  TSimplexID* Ta2 = new TSimplexID[4]{erorId,     6,     2,0};
  TSimplexID* Ta3 = new TSimplexID[4]{erorId,     5,     3,1};
  TSimplexID* Ta4 = new TSimplexID[4]{erorId,     4,     2,0};

  TSimplexID* Ta5 = new TSimplexID[4]{erorId,     7,     3, 5};
  TSimplexID* Ta6 = new TSimplexID[4]{erorId,     6,     2, 4};
  TSimplexID* Ta7 = new TSimplexID[4]{erorId,     7,     5, 1};
  TSimplexID* Ta8 = new TSimplexID[4]{erorId,     6,     4, 0};

  tet_adj.push_back(Ta1);
  tet_adj.push_back(Ta2);
  tet_adj.push_back(Ta3);
  tet_adj.push_back(Ta4);
  tet_adj.push_back(Ta5);
  tet_adj.push_back(Ta6);
  tet_adj.push_back(Ta7);
  tet_adj.push_back(Ta8);

  gmds::BitVector tri_ids;
  std::vector<TSimplexID* > tri_nodes;
  std::vector<TSimplexID* > tri_adj  ;

  SimplexMesh mesh(node_ids,coords,base,
                   tet_ids,tet_nodes,tet_adj,
                   tri_ids,tri_nodes,tri_adj);

  //BALL OF POINT
  SimplicesNode simplicesNode0(&mesh, 0);
  std::vector<TInt> ballOf = simplicesNode0.ballOf();

  //ball
  ASSERT_EQ(ballOf.size(),8);
  ASSERT_EQ(ballOf[0],7);
  ASSERT_EQ(ballOf[1],0);
  ASSERT_EQ(ballOf[2],1);
  ASSERT_EQ(ballOf[3],2);
  ASSERT_EQ(ballOf[4],3);
  ASSERT_EQ(ballOf[5],4);
  ASSERT_EQ(ballOf[6],5);
  ASSERT_EQ(ballOf[7],6);


  //shell with 4 tetra
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1)).size(), 4);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[0], 0);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[1], 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[2], 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[3], 3);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3)).size(), 4);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3))[0], 7);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3))[1], 0);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3))[2], 3);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3))[3], 4);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4)).size(), 4);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4))[0], 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4))[1], 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4))[2], 5);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4))[3], 6);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2)).size(), 4);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[0], 7);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[1], 0);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[2], 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[3], 6);


  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 5)).size(), 4);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 5))[0], 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 5))[1], 3);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 5))[2], 4);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 5))[3], 5);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 6)).size(), 4);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 6))[0], 7);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 6))[1], 4);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 6))[2], 5);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 6))[3], 6);



  //shell with 2 tetra
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 3)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 3))[0], 0);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 3))[1], 3);

  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 4)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 4))[0], 1);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 4))[1], 2);

  ASSERT_EQ(SimplicesNode(&mesh, 5).shell(SimplicesNode(&mesh, 3)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 5).shell(SimplicesNode(&mesh, 3))[0], 3);
  ASSERT_EQ(SimplicesNode(&mesh, 5).shell(SimplicesNode(&mesh, 3))[1], 4);

  ASSERT_EQ(SimplicesNode(&mesh, 5).shell(SimplicesNode(&mesh, 4)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 5).shell(SimplicesNode(&mesh, 4))[0], 2);
  ASSERT_EQ(SimplicesNode(&mesh, 5).shell(SimplicesNode(&mesh, 4))[1], 5);



  ASSERT_EQ(SimplicesNode(&mesh, 6).shell(SimplicesNode(&mesh, 3)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 6).shell(SimplicesNode(&mesh, 3))[0], 4);
  ASSERT_EQ(SimplicesNode(&mesh, 6).shell(SimplicesNode(&mesh, 3))[1], 7);

  ASSERT_EQ(SimplicesNode(&mesh, 6).shell(SimplicesNode(&mesh, 4)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 6).shell(SimplicesNode(&mesh, 4))[0], 6);
  ASSERT_EQ(SimplicesNode(&mesh, 6).shell(SimplicesNode(&mesh, 4))[1], 5);

  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 3)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 3))[0], 0);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 3))[1], 7);

  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4))[0], 1);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4))[1], 6);

}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test_ball_of_two_tetra_with_triangle_between)
{
  //exemple 2 tetra with triagle
  gmds::BitVector node_ids(5);
  node_ids.fillAll();
  std::vector<gmds::math::Point> coords;
  coords.push_back(Point(0.0,0.0,0.0));
  coords.push_back(Point(1.0,0.0,0.0));
  coords.push_back(Point(0.0,0.0,1.0));
  coords.push_back(Point(0.0,1.0,0.0));
  coords.push_back(Point(0.0,-1.0,0.0));

  std::vector<TSimplexID> base({0,0,0,0,1});

  gmds::BitVector tet_ids(2);
  tet_ids.fillAll();
  std::vector<TSimplexID* > tet_nodes;
  TSimplexID* T1 = new TSimplexID[4]{0,1,2,3};
  TSimplexID* T2 = new TSimplexID[4]{0,1,2,4};

  tet_nodes.push_back(T1);
  tet_nodes.push_back(T2);

  std::vector<TSimplexID* > tet_adj;
  TSimplexID erorId = std::numeric_limits<int>::min();
  TSimplexID* Ta1 = new TSimplexID[4]{erorId, erorId, erorId, -1};
  TSimplexID* Ta2 = new TSimplexID[4]{erorId, erorId, erorId, -1};

  tet_adj.push_back(Ta1);
  tet_adj.push_back(Ta2);

  gmds::BitVector tri_ids;
  std::vector<TSimplexID* > tri_nodes;
  std::vector<TSimplexID* > tri_adj  ;

  TSimplexID* t0 = new TSimplexID[4]{0,1,2,0};
  tri_nodes.push_back(nullptr);
  tri_nodes.push_back(t0);

  TSimplexID* t0Adj = new TSimplexID[4]{erorId, erorId, erorId, 1};
  tri_adj.push_back(nullptr);
  tri_adj.push_back(t0Adj);

  SimplexMesh mesh(node_ids,coords,base,
                   tet_ids,tet_nodes,tet_adj,
                   tri_ids,tri_nodes,tri_adj);

  SimplicesNode simplicesNode0(&mesh, 0);
  std::vector<TInt> ballOf = simplicesNode0.ballOf();
  std::vector<TSimplexID> shell = simplicesNode0.shell(SimplicesNode(&mesh,2));

  //ball
  ASSERT_EQ(ballOf.size(),2);
  ASSERT_EQ(ballOf[0],0);
  ASSERT_EQ(ballOf[1],1);

  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf().size(),1);
  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf()[0],0);

  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf().size(),1);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[0],1);

  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[1],1);

  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[1],1);

  //shell
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 1)).size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 1))[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 1))[1],1);

  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 2)).size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 2))[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 2))[1],1);


  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 3)).size(),1);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 3))[0],0);

  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 4)).size(),1);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 4))[0],1);

  ASSERT_EQ(SimplicesNode(&mesh, 3).shell(SimplicesNode(&mesh, 1)).size(),1);
  ASSERT_EQ(SimplicesNode(&mesh, 3).shell(SimplicesNode(&mesh, 1))[0],0);

  ASSERT_EQ(SimplicesNode(&mesh, 3).shell(SimplicesNode(&mesh, 2)).size(),1);
  ASSERT_EQ(SimplicesNode(&mesh, 3).shell(SimplicesNode(&mesh, 2))[0],0);

}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test_ball_of_four_tetra_with_triangles_between_2_of_them)
{
  //exemple  4 tetra with  3 triangle

  gmds::BitVector node_ids(6);
  node_ids.fillAll();
  std::vector<gmds::math::Point> coords;
  coords.push_back(Point(0.0,0.0,0.0));
  coords.push_back(Point(1.0,0.0,0.0));
  coords.push_back(Point(0.0,0.0,1.0));
  coords.push_back(Point(0.0,1.0,0.0));
  coords.push_back(Point(0.0,-1.0,0.0));
  coords.push_back(Point(0.0,0.0,-1.0));

  std::vector<TSimplexID> base({3,0,0,0,1,2});

  gmds::BitVector tet_ids(4);
  tet_ids.fillAll();
  std::vector<TSimplexID* > tet_nodes;
  TSimplexID* T1 = new TSimplexID[4]{0,1,2,3};
  TSimplexID* T2 = new TSimplexID[4]{0,1,2,4};
  TSimplexID* T3 = new TSimplexID[4]{0,1,4,5};
  TSimplexID* T4 = new TSimplexID[4]{0,1,3,5};

  tet_nodes.push_back(T1);
  tet_nodes.push_back(T2);
  tet_nodes.push_back(T3);
  tet_nodes.push_back(T4);

  std::vector<TSimplexID* > tet_adj;
  TSimplexID errorId = std::numeric_limits<int>::min();
  TSimplexID* Ta1 = new TSimplexID[4]{errorId, errorId,     -3, -1};
  TSimplexID* Ta2 = new TSimplexID[4]{errorId, errorId,      2, -1};
  TSimplexID* Ta3 = new TSimplexID[4]{errorId, errorId,     -2,  1};
  TSimplexID* Ta4 = new TSimplexID[4]{errorId, errorId,     -2, -3};

  tet_adj.push_back(Ta1);
  tet_adj.push_back(Ta2);
  tet_adj.push_back(Ta3);
  tet_adj.push_back(Ta4);

  gmds::BitVector tri_ids;
  std::vector<TSimplexID* > tri_nodes;
  std::vector<TSimplexID* > tri_adj  ;

  TSimplexID* t0 = new TSimplexID[4]{0,1,2,0};
  TSimplexID* t1 = new TSimplexID[4]{0,1,5,2};
  TSimplexID* t2 = new TSimplexID[4]{0,1,3,0};

  tri_nodes.push_back(nullptr);
  tri_nodes.push_back(t0);
  tri_nodes.push_back(t1);
  tri_nodes.push_back(t2);

  TSimplexID* t0Adj = new TSimplexID[4]{errorId, errorId,       2, 1};
  TSimplexID* t1Adj = new TSimplexID[4]{errorId, errorId,       1, 3};
  TSimplexID* t2Adj = new TSimplexID[4]{errorId, errorId, errorId, 3};

  tri_adj.push_back(nullptr);
  tri_adj.push_back(t0Adj);
  tri_adj.push_back(t1Adj);
  tri_adj.push_back(t2Adj);

  SimplexMesh mesh(node_ids,coords,base,
                   tet_ids,tet_nodes,tet_adj,
                   tri_ids,tri_nodes,tri_adj);



  //ball
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf().size(),4);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[0],3);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[1],0);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[2],1);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[3],2);

  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf()[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf()[1],3);

  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[0],1);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[1],2);

  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf().size(),4);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[1],1);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[2],2);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[3],3);

  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[1],1);

  ASSERT_EQ(SimplicesNode(&mesh, 5).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 5).ballOf()[0],2);
  ASSERT_EQ(SimplicesNode(&mesh, 5).ballOf()[1],3);

  //shell
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 1)).size(),4);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 1))[0],3);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 1))[1],0);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 1))[2],1);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 1))[3],2);

  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 3)).size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 3))[0],3);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 3))[1],0);


  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 4)).size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 4))[0],1);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 4))[1],2);

  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 5)).size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 5))[0],3);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 5))[1],2);

  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 2)).size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 2))[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 2))[1],1);

  ASSERT_EQ(SimplicesNode(&mesh, 3).shell(SimplicesNode(&mesh, 1)).size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 3).shell(SimplicesNode(&mesh, 1))[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 3).shell(SimplicesNode(&mesh, 1))[1],3);

  ASSERT_EQ(SimplicesNode(&mesh, 3).shell(SimplicesNode(&mesh, 5)).size(),1);
  ASSERT_EQ(SimplicesNode(&mesh, 3).shell(SimplicesNode(&mesh, 5))[0],3);


  //Opposite face of triangles
  ASSERT_EQ(mesh.getOppositeFace(2,1),2); //opposite triangle of triangle  1 by the node 2
  ASSERT_EQ(mesh.getOppositeFace(1,1),errorId);
  ASSERT_EQ(mesh.getOppositeFace(0,1),errorId);

  ASSERT_EQ(mesh.getOppositeFace(5,2),1); //opposite triangle of triangle  2 by the node 5
  ASSERT_EQ(mesh.getOppositeFace(0,2),errorId);
  ASSERT_EQ(mesh.getOppositeFace(1,2),errorId);


  ASSERT_EQ(mesh.getOppositeFace(0,3), errorId); //opposite triangle of triangle 3 by the node 2
  ASSERT_EQ(mesh.getOppositeFace(1,3), errorId);
  ASSERT_EQ(mesh.getOppositeFace(3,3), errorId);



}
