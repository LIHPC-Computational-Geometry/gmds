/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/math/Point.h>
/*----------------------------------------------------------------------------*/
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
#include <gmds/hybridMeshAdapt/CommonInfo.h>
#include <gmds/hybridMeshAdapt/CavityOperator.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace gmds::hybrid;
using namespace gmds::math;
using namespace hybrid;
using namespace hybrid::simplicesCell;
using namespace hybrid::simplicesNode;
using namespace operators;
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test_deleting_tetra)
{
  TSimplexID errorId = std::numeric_limits<int>::min();
  SimplexMesh mesh;
  mesh.addNode(Point(0.0,0.0,0.0));
  mesh.addNode(Point(1.0,0.0,0.0));
  mesh.addNode(Point(0.0,0.0,1.0));
  mesh.addNode(Point(0.0,1.0,0.0));

  TSimplexID tetra0 = mesh.addTetraedre(0,1,2,3);

  //Adding a new Point :  Node 4
  mesh.addNode(Point(0,-1,0));
  //Adding a new Tetra :  Tetra 1
  TSimplexID tetra1 = mesh.addTetraedre(0,1,2,4);

  //Adding a new Point :  Node 5
  mesh.addNode(Point(0,0,-1));
  //Adding a new Tetra :  Tetra 3
  TSimplexID tetra2 = mesh.addTetraedre(0,1,4,5);
  TSimplexID tetra3 = mesh.addTetraedre(0,1,3,5);
  TSimplexID mesh0 = mesh.addTriangle(0, 1, 2);
  TSimplexID mesh1 = mesh.addTriangle(1, 2, 3);
  TSimplexID mesh2 = mesh.addTriangle(1, 2, 4);
  TSimplexID mesh3 = mesh.addTriangle(0, 1, 5);
  TSimplexID mesh4 = mesh.addTriangle(0, 1, 3);
  TSimplexID mesh5 = mesh.addTriangle(0, 1, 4);

  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf().size(), 4);

  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[0], 0);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[1], 1);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[2], 2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[3], 3);

  ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri().size(), 6);
  ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[0], -2);
  ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[1], -5);
  ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[2], -1);
  ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[3], -3);
  ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[4], -6);
  ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[5], -4);

  //TEST SUPPRESSION TET
/*  std::cout << "OKOK" << std::endl;
  mesh.deleteTetra(tetra0);
std::cout << "OKOK" << std::endl;
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf().size(), 3);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[0], 1);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[1], 2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[2], 3);
std::cout << "OKOK" << std::endl;
  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf().size(), 1);


  ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri().size(), 5);
  ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[0], -3);
  ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[1], -6);
  ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[2], -1);
  ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[3], -4);
  ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[4], -5);
std::cout << "OKOK" << std::endl;
  mesh.deleteTetra(tetra1);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf().size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri().size(), 3);

  ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[0], -4);
  ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[1], -6);
  ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[2], -5);

  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf().size(), 0);

  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf().size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf()[0], 3);


  //AJOUT DE TETRA
  TSimplexID idMesh = mesh.addTetraedre(0,1,2,3);


  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf().size(), 3);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[0], 2);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[1], 3);
  ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[2], idMesh);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf().size(), 3);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[0], 2);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[1], 3);
  ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[2], 1);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf().size(), 1);
  ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[0], idMesh);
  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf().size(), 2);
  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf()[0], 3);
  ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf()[1], idMesh);


  ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri().size(), 3);
  ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[0], -4);
  ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[1], -6);
  ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[2], -5);
    //AJOUT DE TETRA

    idMesh = mesh.addTetraedre(0,1,2,4);
    ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf().size(), 4);
    ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[0], 2);
    ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[1], 0);
    ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[2], 1);
    ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[3], 3);
    ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf().size(), 4);
    ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[0], 2);
    ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[1], 0);
    ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[2], 1);
    ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf()[3], 3);
    ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf().size(), 2);
    ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[0], 1);
    ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf()[1], idMesh);
    ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf().size(), 2);
    ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[0], 2);
    ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf()[1], idMesh);


    ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri().size(), 3);
    ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[0], -4);
    ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[1], -6);
    ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[2], -5);

    ASSERT_EQ(SimplicesNode(&mesh, 2).linksTri().size(), 0);

    TSimplexID triangle = mesh.addTriangle(0,1,2);
    ASSERT_EQ(SimplicesNode(&mesh, 2).linksTri().size(), 1);
    ASSERT_EQ(SimplicesNode(&mesh, 0).linksTri().size(), 4);
    ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri().size(), 4);
    ASSERT_EQ(SimplicesNode(&mesh, 5).linksTri().size(), 1);

    ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf().size(), 4);
    ASSERT_EQ(SimplicesNode(&mesh, 1).ballOf().size(), 4);
    ASSERT_EQ(SimplicesNode(&mesh, 2).ballOf().size(), 2);
    ASSERT_EQ(SimplicesNode(&mesh, 5).ballOf().size(), 2);
    ASSERT_EQ(SimplicesNode(&mesh, 3).ballOf().size(), 2);
    ASSERT_EQ(SimplicesNode(&mesh, 4).ballOf().size(), 2);

    //On delete les trianglesm
    mesh.deleteTriangle(triangle);
    ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri().size(), 3);
    ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[0], -4);
    ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[1], -6);
    ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[2], -5);


    triangle = mesh.addTriangle(0,1,2);
    ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri().size(), 4);
    ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[0], -4);
    ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[1], -6);
    ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[2], -triangle);
    ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[3], -5);
*/
/*
    mesh.deleteNode(3);

    ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf().size(), 2);
    ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[1], 0);
    ASSERT_EQ(SimplicesNode(&mesh, 0).ballOf()[0], 2);
*/
    /*ASSERT_EQ(SimplicesNode(&mesh, 0).linksTri().size(), 1);
    ASSERT_EQ(SimplicesNode(&mesh, 0).linksTri()[0], -1);

    ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri().size(), 2);
    ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[0], -3);
    ASSERT_EQ(SimplicesNode(&mesh, 1).linksTri()[1], -1);*/

}
/*----------------------------------------------------------------------------*/
