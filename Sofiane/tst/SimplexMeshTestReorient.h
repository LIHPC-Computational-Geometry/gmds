/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/math/Point.h>
/*----------------------------------------------------------------------------*/
#include <gmds/sofiane/SimplexMesh.h>
#include <gmds/sofiane/CommonInfo.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::hybrid;
using namespace gmds::math;
using namespace hybrid;
using namespace hybrid::simplicesCell;
using namespace hybrid::simplicesNode;
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test_Reorient_Tet)
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



  ASSERT_EQ(SimplicesCell(&mesh, 0).getVolumeOfCell() /
  std::fabs(SimplicesCell(&mesh, 0 ).getVolumeOfCell()), -1);

  /*Tet 0 is already well oriented (normal out of the Tet) so the next function won't change anything*/
  SimplicesCell(&mesh, 0).reorientTet();

  ASSERT_EQ(SimplicesCell(&mesh, 0).getVolumeOfCell() /
  std::fabs(SimplicesCell(&mesh, 0 ).getVolumeOfCell()), 1);

  /*create a Tetra not well oriented*/
  TInt indexNode = mesh.addNode(Point(0,-1,0));
  mesh.addTetraedre(0,1,2,4);


  ASSERT_EQ(SimplicesCell(&mesh, 1).getVolumeOfCell() /
  std::fabs(SimplicesCell(&mesh, 1 ).getVolumeOfCell()), 1);

  /*change the orientation node of tetra 1 in order to get the normal out of this Tetra*/
  //SimplicesCell(&mesh, 1).reorientTet();

  /*reorient all the tetra of the mesh (normal out of the Tetra)*/
  mesh.reorientAllTetra();

  ASSERT_EQ(SimplicesCell(&mesh, 0).getVolumeOfCell() /
  std::fabs(SimplicesCell(&mesh, 0 ).getVolumeOfCell()), 1);

  ASSERT_EQ(SimplicesCell(&mesh, 1).getVolumeOfCell() /
  std::fabs(SimplicesCell(&mesh, 1 ).getVolumeOfCell()), 1);


  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf().size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[1],1);

  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf().size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[0],0);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[1],1);

}
