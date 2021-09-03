/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/math/Point.h>
/*----------------------------------------------------------------------------*/
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
#include <gmds/hybridMeshAdapt/CommonInfo.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::hybrid;
using namespace gmds::math;
using namespace hybrid;
using namespace hybrid::simplicesCell;
using namespace hybrid::simplicesNode;

/*----------------------------------------------------------------------------*/
/*TEST(SimplexMeshTestClass, test_add_Node_to_tetra)
{
  std::vector<gmds::math::Point> coords;
  coords.push_back(Point(0.0,0.0,0.0));
  coords.push_back(Point(1.0,0.0,0.0));
  coords.push_back(Point(0.0,0.0,1.0));
  coords.push_back(Point(0.0,1.0,0.0));

  std::vector<TSimplexID* > tet_nodes;
  TSimplexID* T1 = new TSimplexID[4]{0,1,2,3};

  tet_nodes.push_back(T1);

  std::vector<TSimplexID* > tet_adj;
  TSimplexID erorId = std::numeric_limits<int>::min();
  TSimplexID* Ta1 = new TSimplexID[4]{erorId,erorId,erorId,erorId};
  tet_adj.push_back(Ta1);

  std::vector<TSimplexID* > tri_nodes;
  std::vector<TSimplexID* > tri_adj  ;

  SimplexMesh mesh(coords,
                   tet_nodes,tet_adj,
                   tri_nodes,tri_adj);


  SimplicesNode simplicesNode0(&mesh, 0);
  std::vector<TInt> ballOf = simplicesNode0.ballOf();

  ASSERT_EQ(ballOf.size(),1);
  ASSERT_EQ(ballOf[0],0);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[0], 0);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[0], 0);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3))[0], 0);

  //Adding a new Point :  Node 4
  mesh.addNode(Point(0,-1,0));
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4)).size(), 0);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 4)).size(), 0);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4)).size(), 0);


}*/
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test_add_Tet_to_tetra)
{
  SimplexMesh mesh;
  mesh.addNode(Point(0.0,0.0,0.0));
  mesh.addNode(Point(1.0,0.0,0.0));
  mesh.addNode(Point(0.0,0.0,1.0));
  mesh.addNode(Point(0.0,1.0,0.0));
  mesh.addTetraedre(0,1,2,3);

  SimplicesNode simplicesNode0(&mesh, 0);
  std::vector<TInt> ballOf = simplicesNode0.ballOf();

  ASSERT_EQ(ballOf.size(),1);
  ASSERT_EQ(ballOf[0],0);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[0], 0);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[0], 0);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3))[0], 0);

  /*Adding a new Point :  Node 4*/
  mesh.addNode(Point(0,-1,0));

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4)).size(), 0);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 4)).size(), 0);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4)).size(), 0);

  /*Adding a new Tetra :  Tetra 1*/
  mesh.addTetraedre(0,1,2,4);

  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[1],1);

  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[1],1);

  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[1],1);

  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf().size(),1);
  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf()[0],0);

  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf().size(),1);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[0],1);


  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1)).size(), 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[0], 0);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[1], 1);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2)).size(), 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[0], 0);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[1], 1);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4))[0], 1);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3))[0], 0);

  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 3)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 3))[0], 0);

  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4))[0], 1);

  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 3)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 3))[0], 0);

  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 2)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 2))[0], 0);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 2))[1], 1);

  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 4)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 4))[0], 1);


  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4))[0], 1);

}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test_add_two_Tet_to_tetra)
{
  SimplexMesh mesh;
  mesh.addNode(Point(0.0,0.0,0.0));
  mesh.addNode(Point(1.0,0.0,0.0));
  mesh.addNode(Point(0.0,0.0,1.0));
  mesh.addNode(Point(0.0,1.0,0.0));

  mesh.addTetraedre(0,1,2,3);

  SimplicesNode simplicesNode0(&mesh, 0);
  std::vector<TInt> ballOf = simplicesNode0.ballOf();

  ASSERT_EQ(ballOf.size(),1);
  ASSERT_EQ(ballOf[0],0);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[0], 0);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[0], 0);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3))[0], 0);

  /*Adding a new Point :  Node 4 & Node 5*/
  mesh.addNode(Point(0,-1,  0));
  mesh.addNode(Point(0, 0, -1));

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4)).size(), 0);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 4)).size(), 0);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4)).size(), 0);

  /*Adding a new Tetra :  Tetra 1 & Tetra 2*/
  mesh.addTetraedre(0,1,2,4);
  mesh.addTetraedre(0,1,3,5);

  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf().size(),3);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[1],1);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[2],2);

  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf().size(),3);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[1],1);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[2],2);

  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[1],1);

  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf()[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf()[1],2);

  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf().size(),1);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[0],1);

  ASSERT_EQ(SimplicesNode(&mesh, 5).ballOf().size(),1);
  ASSERT_EQ(SimplicesNode(&mesh, 5).ballOf()[0],2);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1)).size(), 3);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[0], 0);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[1], 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[2], 2);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2)).size(), 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[0], 0);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[1], 1);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4))[0], 1);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3)).size(), 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3))[0], 0);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3))[1], 2);

  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 3)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 3))[0], 0);

  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4))[0], 1);

  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 2)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 2))[0], 0);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 2))[1], 1);

  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 4)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 4))[0], 1);

  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4))[0], 1);

  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 3)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 3))[0], 0);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 3))[1], 2);

  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 5)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 5))[0], 2);

  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 5)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 5))[0], 2);

  ASSERT_EQ(SimplicesNode(&mesh, 3).shell(SimplicesNode(&mesh, 5)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 3).shell(SimplicesNode(&mesh, 5))[0], 2);
}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test_add_seven_Tet_to_tetra)
{
  SimplexMesh mesh;
  mesh.addNode(Point(0.0,0.0,0.0));
  mesh.addNode(Point(1.0,0.0,0.0));
  mesh.addNode(Point(0.0,0.0,1.0));
  mesh.addNode(Point(0.0,1.0,0.0));

  mesh.addTetraedre(0,1,2,3);

  SimplicesNode simplicesNode0(&mesh, 0);
  std::vector<TInt> ballOf = simplicesNode0.ballOf();

  ASSERT_EQ(ballOf.size(),1);
  ASSERT_EQ(ballOf[0],0);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[0], 0);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[0], 0);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3))[0], 0);

  /*Adding a new Point :  Node 4 & Node 5*/
  mesh.addNode(Point(0,-1,  0));
  mesh.addNode(Point(0, 0, -1));

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4)).size(), 0);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 4)).size(), 0);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4)).size(), 0);

  mesh.addTetraedre(0,1,2,4); //T1
  mesh.addTetraedre(0,1,4,5); //T2
  mesh.addTetraedre(0,1,3,5); //T3

  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf().size(),4);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[1],1);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[2],2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[3],3);

  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf().size(),4);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[1],1);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[2],2);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[3],3);

  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[1],1);

  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf()[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf()[1],3);

  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[0],1);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[1],2);

  ASSERT_EQ(SimplicesNode(&mesh, 5).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 5).ballOf()[0],2);
  ASSERT_EQ(SimplicesNode(&mesh, 5).ballOf()[1],3);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1)).size(), 4);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[0], 0);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[1], 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[2], 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[3], 3);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2)).size(), 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[0], 0);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[1], 1);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4)).size(), 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4))[0], 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4))[1], 2);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3)).size(), 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3))[0], 0);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3))[1], 3);

  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 3)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 3))[0], 0);

  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4))[0], 1);

  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 2)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 2))[0], 0);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 2))[1], 1);

  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 4)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 4))[0], 1);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 4))[1], 2);

  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4))[0], 1);

  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 3)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 3))[0], 0);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 3))[1], 3);

  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 5)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 5))[0], 2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 5))[1], 3);

  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 5)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 5))[0], 2);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 5))[1], 3);

  ASSERT_EQ(SimplicesNode(&mesh, 3).shell(SimplicesNode(&mesh, 5)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 3).shell(SimplicesNode(&mesh, 5))[0], 3);
}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test_add_Triangle_to_tetra)
{
  /*std::vector<gmds::math::Point> coords;
  coords.push_back(Point(0.0,0.0,0.0));
  coords.push_back(Point(1.0,0.0,0.0));
  coords.push_back(Point(0.0,0.0,1.0));
  coords.push_back(Point(0.0,1.0,0.0));

  std::vector<TSimplexID* > tet_nodes;
  TSimplexID* T1 = new TSimplexID[4]{0,1,2,3};

  tet_nodes.push_back(T1);
  std::vector<TSimplexID* > tri_nodes;

  SimplexMesh mesh(coords,
                   tet_nodes,
                   tri_nodes);


  SimplicesNode simplicesNode0(&mesh, 0);
  std::vector<TInt> ballOf = simplicesNode0.ballOf();

  ASSERT_EQ(ballOf.size(),1);
  ASSERT_EQ(ballOf[0],0);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[0], 0);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[0], 0);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3))[0], 0);
*/
  /*Adding a new Point :  Node 4 */
  /*mesh.addNode(Point(0,-1,0));

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4)).size(), 0);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 4)).size(), 0);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4)).size(), 0);
*/
  /*Adding a new Tetra :  Tetra 1 */
  /*mesh.addTetraedre(0,1,2,4);

  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[0],1);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[1],0);

  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[0],1);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[1],0);

  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[0],1);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[1],0);

  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf().size(),1);
  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf()[0],0);

  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf().size(),1);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[0],1);


  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1)).size(), 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[0], 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 1))[1], 0);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2)).size(), 2);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[0], 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 2))[1], 0);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 4))[0], 1);

  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3)).size(), 1);
  ASSERT_EQ(simplicesNode0.shell(SimplicesNode(&mesh, 3))[0], 0);

  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 3)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 3))[0], 0);

  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4))[0], 1);

  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 3)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 3))[0], 0);

  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 2)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 2))[0], 1);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 2))[1], 0);

  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 4)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 4))[0], 1);


  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 4))[0], 1);
*/
  /*add triangle*/
  /*mesh.addTriangle(0, 1, 2);

  mesh.addNode(Point(0,0,-1));
  mesh.addTriangle(0, 1, 5);

  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[0],1);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[1],0);

  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[0],1);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[1],0);

  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[0],1);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[1],0);

  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf().size(),1);
  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf()[0],0);

  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf().size(),1);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[0],1);

  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 1)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 1))[0], 1);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 1))[1], 0);

  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 2)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 2))[0], 1);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 2))[1], 0);

  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 3)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 3))[0], 0);

  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 4)).size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 4))[0], 1);

  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 2)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 2))[0], 1);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 2))[1], 0);

  mesh.addTetraedre(0, 1, 4, 5);
  mesh.addTetraedre(0, 1, 3, 5);

  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf().size(),4);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[0],3);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[1],0);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[2],1);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[3],2);

  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf().size(),4);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[0],3);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[1],0);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[2],1);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[3],2);

  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[0],1);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[1],0);

  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf()[0],3);
  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf()[1],0);

  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[0],2);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[1],1);

  ASSERT_EQ(SimplicesNode(&mesh, 5).ballOf().size(),2);
  ASSERT_EQ(SimplicesNode(&mesh, 5).ballOf()[0],3);
  ASSERT_EQ(SimplicesNode(&mesh, 5).ballOf()[1],2);

  TInt pt = mesh.addNode(Point(-1, 0, 0));

  mesh.addTetraedre(0, 2, 3, 6);
  mesh.addTetraedre(0, 6, 3, 5);
  mesh.addTetraedre(0, 6, 4, 2);
  mesh.addTetraedre(0, 6, 4, 5);

  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf().size(),8);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[0],4);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[1],0);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[2],1);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[3],2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[4],3);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[5],5);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[6],7);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[7],6);


  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf().size(),4);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[0],3);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[1],0);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[2],1);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[3],2);

  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf().size(),4);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[0],6);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[1],4);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[2],0);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[3],1);

  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf().size(),4);
  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf()[0],5);
  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf()[1],4);
  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf()[2],0);
  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf()[3],3);

  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf().size(),4);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[0],7);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[1],6);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[2],1);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[3],2);

  ASSERT_EQ(SimplicesNode(&mesh, 5).ballOf().size(),4);
  ASSERT_EQ(SimplicesNode(&mesh, 5).ballOf()[0],7);
  ASSERT_EQ(SimplicesNode(&mesh, 5).ballOf()[1],2);
  ASSERT_EQ(SimplicesNode(&mesh, 5).ballOf()[2],3);
  ASSERT_EQ(SimplicesNode(&mesh, 5).ballOf()[3],5);

  ASSERT_EQ(SimplicesNode(&mesh, 6).ballOf().size(),4);
  ASSERT_EQ(SimplicesNode(&mesh, 6).ballOf()[0],7);
  ASSERT_EQ(SimplicesNode(&mesh, 6).ballOf()[1],6);
  ASSERT_EQ(SimplicesNode(&mesh, 6).ballOf()[2],4);
  ASSERT_EQ(SimplicesNode(&mesh, 6).ballOf()[3],5);


  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 1)).size(), 4);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 1))[0], 0);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 1))[1], 1);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 1))[2], 2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 1))[3], 3);

  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 6)).size(), 4);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 6))[0], 4);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 6))[1], 5);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 6))[2], 7);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 6))[3], 6);


  mesh.addTriangle(0,6,2);
  mesh.addTriangle(0,6,5);

  //5 eme triangle
  mesh.addTriangle(0,6,3);
  mesh.addTriangle(0,1,3);
  mesh.addTriangle(0,6,4);
  mesh.addTriangle(0,1,4);

  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf().size(),8);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 1)).size(), 4);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 6)).size(), 4);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 3)).size(), 4);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 4)).size(), 4);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 5)).size(), 4);
  ASSERT_EQ(SimplicesNode(&mesh, 0).shell(SimplicesNode(&mesh, 2)).size(), 4);

  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 2)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 6)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 6).shell(SimplicesNode(&mesh, 5)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 5).shell(SimplicesNode(&mesh, 1)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 3).shell(SimplicesNode(&mesh, 1)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 1).shell(SimplicesNode(&mesh, 4)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 4).shell(SimplicesNode(&mesh, 6)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 6).shell(SimplicesNode(&mesh, 3)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 2).shell(SimplicesNode(&mesh, 3)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 3).shell(SimplicesNode(&mesh, 5)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 5).shell(SimplicesNode(&mesh, 4)).size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 4).shell(SimplicesNode(&mesh, 2)).size(), 2);

  //pt 7
  mesh.addNode(Point(-2, 0, 0));
  mesh.addTetraedre(6, 7, 2, 4);

  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf().size(),5);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[0],8);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[1],6);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[2],7);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[3],2);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[4],1);

  ASSERT_EQ(SimplicesNode(&mesh, 4).shell(SimplicesNode(&mesh, 2)).size(), 3);
  ASSERT_EQ(SimplicesNode(&mesh, 4).shell(SimplicesNode(&mesh, 2))[0], 8);
  ASSERT_EQ(SimplicesNode(&mesh, 4).shell(SimplicesNode(&mesh, 2))[1], 6);
  ASSERT_EQ(SimplicesNode(&mesh, 4).shell(SimplicesNode(&mesh, 2))[2], 1);

  mesh.addTriangle(4, 7, 2);

  ASSERT_EQ(SimplicesNode(&mesh, 7).ballOf().size(),1);

  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf().size(),5);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[0],8);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[1],6);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[2],7);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[3],2);
  ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[4],1);

  ASSERT_EQ(SimplicesNode(&mesh, 4).shell(SimplicesNode(&mesh, 2)).size(), 3);
  ASSERT_EQ(SimplicesNode(&mesh, 4).shell(SimplicesNode(&mesh, 2))[0], 8);
  ASSERT_EQ(SimplicesNode(&mesh, 4).shell(SimplicesNode(&mesh, 2))[1], 6);
  ASSERT_EQ(SimplicesNode(&mesh, 4).shell(SimplicesNode(&mesh, 2))[2], 1);
*/
}
/*----------------------------------------------------------------------------*/

TEST(SimplexMeshTestClass, test_point_in_two_tetra)
{/*
  TSimplexID errorId = std::numeric_limits<int>::min();
  SimplexMesh mesh = SimplexMesh();
  mesh.addNode(math::Point(0.0, 0.0, 0.0)); // Node 0
  mesh.addNode(math::Point(1.0, 0.0, 0.0)); // Node 1
  mesh.addNode(math::Point(0.0, 0.0, 1.0)); // Node 2
  mesh.addNode(math::Point(0.0, 1.0, 0.0)); // Node 3
  mesh.addNode(math::Point(0.0, -1.0, 0.0)); // Node 4
  mesh.addNode(math::Point(0.0, 0.0, -1.0)); // Node 5
  mesh.addNode(math::Point(-1.0, 0.0, 0.0)); // Node 6

  mesh.addTetraedre(0, 1 , 2, 3); // Tetra 0
  mesh.addTetraedre(0, 1 , 2, 4); // Tetra 1
  mesh.addTetraedre(0, 1 , 4, 5); // Tetra 2
  mesh.addTetraedre(0, 1 , 3, 5); // Tetra 3
  mesh.addTetraedre(0, 3 , 5, 6); // Tetra 4
  mesh.addTetraedre(0, 2 , 3, 6); // Tetra 5
  mesh.addTetraedre(0, 2 , 4, 6); // Tetra 6
  mesh.addTetraedre(0, 4 , 5, 6); // Tetra 7




  //simplexContaining has changed, take it into account later...
  //Pt à l'interieur du tetra 0
  TSimplexID simplex = 0;
  mesh.simplexContaining(Point(0.0, 0.0, 0.0), simplex);
  ASSERT_EQ(simplex, 0);
  mesh.simplexContaining(Point(0.001, 0.001, 0.001), simplex);
  ASSERT_EQ(simplex, 0);
  mesh.simplexContaining(Point(0.001, -0.001, 0.001), simplex);
  ASSERT_EQ(simplex, 1);
  mesh.simplexContaining(Point(0.001, -0.001, -0.001), simplex);
  ASSERT_EQ(simplex, 2);
  mesh.simplexContaining(Point(0.001, 0.001, -0.001), simplex);
  ASSERT_EQ(simplex, 3);
  mesh.simplexContaining(Point(-0.001, 0.001, -0.001), simplex);
  ASSERT_EQ(simplex, 4);
  mesh.simplexContaining(Point(-0.001, 0.001, 0.001), simplex);
  ASSERT_EQ(simplex, 5);
  mesh.simplexContaining(Point(-0.001, -0.001, 0.001), simplex);
  ASSERT_EQ(simplex, 6);
  mesh.simplexContaining(Point(-0.001, -0.001, -0.001), simplex);
  ASSERT_EQ(simplex, 7);

  mesh.simplexContaining(Point(0.0, 1.0, 0.0), simplex);
  ASSERT_EQ(simplex, 0);
  mesh.simplexContaining(Point(0.0, 0.0, 1.0), simplex);
  ASSERT_EQ(simplex, 0);


  //Pt à l'interieur du tetra 1
  mesh.simplexContaining(Point(0.0, -1.0, 0.0), simplex);
  ASSERT_EQ(simplex, 1);
  mesh.simplexContaining(Point(0.2, -0.5, 0.2), simplex);
  ASSERT_EQ(simplex, 1);
  mesh.simplexContaining(Point(0.25, -0.25, 0.25), simplex);
  ASSERT_EQ(simplex, 1);



  //Pt à l'exterieur du tetra  0 et 1
  mesh.simplexContaining(Point(-1.1, 0.0, 0.0), simplex);
  ASSERT_EQ(simplex, errorId);
  mesh.simplexContaining(Point(0.0, -1.1, 0.0), simplex);
  ASSERT_EQ(simplex, errorId);
  mesh.simplexContaining(Point(0.0, 0.0, -1.1), simplex);
  ASSERT_EQ(simplex, errorId);
  mesh.simplexContaining(Point(-1.1, -1.1, 0.0), simplex);
  ASSERT_EQ(simplex, errorId);
  mesh.simplexContaining(Point(0.0, -1.1, -1.1), simplex);
  ASSERT_EQ(simplex, errorId);
  mesh.simplexContaining(Point(-1.0, 0.0, -1.1), simplex);
  ASSERT_EQ(simplex, errorId);

  //On ajoute unn triangle entre les tetra du mesh pour voir si le code fonctionne toujours
  mesh.addTriangle(0, 1, 2);
  //On recommence les meme test que plus haut pour voir si ajouter un triangle ne change rien ...

  mesh.simplexContaining(Point(0.0, 0.0, 0.0), simplex);
  ASSERT_EQ(simplex, 0);
  mesh.simplexContaining(Point(0.0, 1.0, 0.0), simplex);
  ASSERT_EQ(simplex, 0);
  mesh.simplexContaining(Point(0.0, 0.0, 1.0), simplex);
  ASSERT_EQ(simplex, 0);
  mesh.simplexContaining(Point(0.25, 0.25, 0.25), simplex);


  //Pt à l'interieur du tetra 1
  mesh.simplexContaining(Point(0.0, -1.0, 0.0), simplex);
  ASSERT_EQ(simplex, 1);
  mesh.simplexContaining(Point(0.2, -0.5, 0.2), simplex);
  ASSERT_EQ(simplex, 1);
  mesh.simplexContaining(Point(0.25, -0.25, 0.25), simplex);
  ASSERT_EQ(simplex, 1);

  mesh.simplexContaining(Point(-1.1, 0.0, 0.0), simplex);
  ASSERT_EQ(simplex, errorId);
  mesh.simplexContaining(Point(0.0, -1.1, 0.0), simplex);
  ASSERT_EQ(simplex, errorId);
  mesh.simplexContaining(Point(0.0, 0.0, -1.1), simplex);
  ASSERT_EQ(simplex, errorId);
  mesh.simplexContaining(Point(-1.0, -1.1, 0.0), simplex);
  ASSERT_EQ(simplex, errorId);
  mesh.simplexContaining(Point(0.0, -1.1, -1.0), simplex);
  ASSERT_EQ(simplex, errorId);
  mesh.simplexContaining(Point(-1.0, 0.0, -1.1), simplex);
  ASSERT_EQ(simplex, errorId);

  //On ajoute un triangle en bord de mesh pour voir si le code précédent fonctionne
  mesh.addTriangle(1, 2, 3);
  //Pt à l'interieur du tetra 0
  mesh.simplexContaining(Point(0.0, 0.0, 0.0), simplex);
  ASSERT_EQ(simplex, 0);
  mesh.simplexContaining(Point(0.0, 1.0, 0.0), simplex);
  ASSERT_EQ(simplex, 0);
  mesh.simplexContaining(Point(0.0, 0.0, 1.0), simplex);
  ASSERT_EQ(simplex, 0);
  mesh.simplexContaining(Point(0.25, 0.25, 0.25), simplex);


  //Pt à l'interieur du tetra 1
  mesh.simplexContaining(Point(0.0, -1.0, 0.0), simplex);
  ASSERT_EQ(simplex, 1);
  mesh.simplexContaining(Point(0.2, -0.5, 0.2), simplex);
  ASSERT_EQ(simplex, 1);
  mesh.simplexContaining(Point(0.25, -0.25, 0.25), simplex);
  ASSERT_EQ(simplex, 1);

  //Pt à l'exterieur du tetra  0 et 1
  mesh.simplexContaining(Point(-1.1, 0.0, 0.0), simplex);
  ASSERT_EQ(simplex, errorId);
  mesh.simplexContaining(Point(0.0, -1.1, 0.0), simplex);
  ASSERT_EQ(simplex, errorId);
  mesh.simplexContaining(Point(0.0, 0.0, -1.1), simplex);
  ASSERT_EQ(simplex, errorId);
  mesh.simplexContaining(Point(-1.0, -1.1, 0.0), simplex);
  ASSERT_EQ(simplex, errorId);
  mesh.simplexContaining(Point(0.0, -1.1, -1.0), simplex);
  ASSERT_EQ(simplex, errorId);
  mesh.simplexContaining(Point(-1.0, 0.0, -1.0), simplex);
  ASSERT_EQ(simplex, errorId);
*/
}
