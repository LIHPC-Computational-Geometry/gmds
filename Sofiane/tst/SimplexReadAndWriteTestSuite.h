/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <gmds/sofiane/PointSmoothing.h>
#include <gmds/sofiane/PointInsertion.h>
#include <gmds/sofiane/EdgeCollapse.h>
#include <gmds/sofiane/ICriterion.h>
#include <gmds/sofiane/SimplexMesh.h>
#include <gmds/sofiane/ISimplexMeshIOService.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace hybrid;
using namespace operators;
using namespace simplicesNode;
using namespace simplicesCell;
/*----------------------------------------------------------------------------*/
TEST(SimplexReadAndWriteTestClass, DISABLED_test_read_CubeVtk_for_simplexMesh)
{
  /*gmds::hybrid::SimplexMesh simplexMesh;
  gmds::ISimplexMeshIOService ioService(&simplexMesh);

    std::string dir(TEST_SAMPLES_DIR);
    std::string vtk_file = dir+"/simpleCube.vtk";

  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.read(vtk_file);

  ASSERT_EQ(simplexMesh.getNbNodes(), 46);
  ASSERT_EQ(simplexMesh.getNbTetra(), 100);*/
}
/*----------------------------------------------------------------------------*/
TEST(MeshClass, testWriterVTK_for_simplexMesh)
{
    SimplexMesh m;
    m.addNode(math::Point(0.0, 0.0, 0.0));
    m.addNode(math::Point(1.0, 0.0, 0.0));
    m.addNode(math::Point(0.0, 0.0, 1.0));
    m.addNode(math::Point(0.0, 1.0, 0.0));
    m.addNode(math::Point(0.0, -1.0, 0.0));
    m.addNode(math::Point(0.0, 0.0, -1.0));

    m.addTetraedre(0,1,2,3);
    m.addTetraedre(0,1,2,4);
    m.addTetraedre(0,1,5,4);

    m.addTriangle(0,5,3);

    //m.addTriangle(0,3,5);
    gmds::ISimplexMeshIOService ioService(&m);
    gmds::VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::R);
    //vtkWriter.setDataOptions(gmds::N|gmds::R); // a voir quand je rajouterai les valeurs aux sommets triangle tetra cell ....
    vtkWriter.write("test_simplex_mesh.vtk");
}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test03_point_smoothing_in_eight_tet_writer)
{
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



    TInt newNode = mesh.addNode(math::Point(0.1, -0.5, -0.10));
    CriterionRAIS criterionRAIS(new VolumeCriterion());
    PointSmoothing ps(&mesh, SimplicesNode(&mesh, 0), SimplicesNode(&mesh, newNode), criterionRAIS);

    gmds::ISimplexMeshIOService ioService(&mesh);
    gmds::VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::R);
    vtkWriter.write("eight_Tet_Point_Smoothing00.vtk");

}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test00_point_insertion_on_two_tetra)
{
    SimplexMesh mesh = SimplexMesh();
    mesh.addNode(math::Point(0.0, 0.0, 0.0)); // Node 0
    mesh.addNode(math::Point(1.0, 0.0, 0.0)); // Node 1
    mesh.addNode(math::Point(0.0, 0.0, 1.0)); // Node 2
    mesh.addNode(math::Point(0.0, 1.0, 0.0)); // Node 3
    mesh.addNode(math::Point(0.0, -1.0, 0.0)); // Node 4

    mesh.addTetraedre(0, 1 , 2, 3); // Tetra 0
    mesh.addTetraedre(0, 1 , 2, 4); // Tetra 1

    TInt newNode0 = mesh.addNode(math::Point(0.2, 0.5, 0.2));
    CriterionRAIS criterionRAIS(new VolumeCriterion());
    PointInsertion(&mesh, SimplicesNode(&mesh, newNode0), criterionRAIS);
    gmds::ISimplexMeshIOService ioService(&mesh);
    gmds::VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::R);
    vtkWriter.write("eight_Tet_point_insertion00.vtk");

}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test01_point_insertion_on_two_tetra)
{
    SimplexMesh mesh = SimplexMesh();
    mesh.addNode(math::Point(0.0, 0.0, 0.0)); // Node 0
    mesh.addNode(math::Point(1.0, 0.0, 0.0)); // Node 1
    mesh.addNode(math::Point(0.0, 0.0, 1.0)); // Node 2
    mesh.addNode(math::Point(0.0, 1.0, 0.0)); // Node 3
    mesh.addNode(math::Point(0.0, -1.0, 0.0)); // Node 4

    mesh.addTetraedre(0, 1 , 2, 3); // Tetra 0
    mesh.addTetraedre(0, 1 , 2, 4); // Tetra 1

    TInt newNode0 = mesh.addNode(math::Point(0.2, 0.5, 0.2));
    TInt newNode2 = mesh.addNode(math::Point(0.1, -0.3, 0.1));
    CriterionRAIS criterionRAIS(new VolumeCriterion());

    PointInsertion pi0(&mesh, SimplicesNode(&mesh, newNode0), criterionRAIS);
    PointInsertion pi2(&mesh, SimplicesNode(&mesh, newNode2), criterionRAIS);


    gmds::ISimplexMeshIOService ioService(&mesh);
    gmds::VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::R);
    vtkWriter.write("eight_Tet_point_insertion01.vtk");
}
TEST(SimplexMeshTestClass, test02_point_insertion_on_two_tetra)
{
    SimplexMesh mesh = SimplexMesh();
    mesh.addNode(math::Point(0.0, 0.0, 0.0)); // Node 0
    mesh.addNode(math::Point(1.0, 0.0, 0.0)); // Node 1
    mesh.addNode(math::Point(0.0, 0.0, 1.0)); // Node 2
    mesh.addNode(math::Point(0.0, 1.0, 0.0)); // Node 3
    mesh.addNode(math::Point(0.0, -1.0, 0.0)); // Node 4

    mesh.addTetraedre(0, 1 , 2, 3); // Tetra 0
    mesh.addTetraedre(0, 1 , 2, 4); // Tetra 1

    TInt newNode0 = mesh.addNode(math::Point(0.2, 0.5, 0.2));
    TInt newNode2 = mesh.addNode(math::Point(0.1, 0.2, 0.1));
    TInt newNode3 = mesh.addNode(math::Point(0.3, 0.15, 0.5));

    CriterionRAIS criterionRAIS(new VolumeCriterion());
    PointInsertion pi0(&mesh, SimplicesNode(&mesh, newNode0), criterionRAIS);
    PointInsertion pi2(&mesh, SimplicesNode(&mesh, newNode2), criterionRAIS);
    PointInsertion pi3(&mesh, SimplicesNode(&mesh, newNode3), criterionRAIS);

    //TInt newNode4 = mesh.addNode(math::Point(0.5, 0.15, 0.5));
    //PointInsertion<VolumeCriterion> pi4(&mesh, SimplicesNode(&mesh, newNode2));

    gmds::ISimplexMeshIOService ioService(&mesh);
    gmds::VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::R);
    vtkWriter.write("eight_Tet_point_insertion02.vtk");
}
/******************************************************************************/
TEST(SimplexMeshTestClass, test00_edge_collapse_on_two_tetra)
{
  SimplexMesh mesh = SimplexMesh();
  mesh.addNode(math::Point(0.0, 0.0, 0.0)); // Node 0
  mesh.addNode(math::Point(1.0, 0.0, 0.0)); // Node 1
  mesh.addNode(math::Point(0.0, 0.0, 1.0)); // Node 2
  mesh.addNode(math::Point(0.0, 1.0, 0.0)); // Node 3
  mesh.addNode(math::Point(0.0, -1.0, 0.0)); // Node 4

  mesh.addTetraedre(0, 1 , 2, 3); // Tetra 0
  mesh.addTetraedre(0, 1 , 2, 4); // Tetra 1

  CriterionRAIS criterionRAIS(new VolumeCriterion());
  TInt newNode0 = mesh.addNode(math::Point(0.2, 0.5, 0.2));
  TInt newNode2 = mesh.addNode(math::Point(0.1, 0.2, 0.1));
  TInt newNode3 = mesh.addNode(math::Point(0.3, 0.15, 0.5));

  PointInsertion pi0(&mesh, SimplicesNode(&mesh, newNode0), criterionRAIS);
  PointInsertion pi2(&mesh, SimplicesNode(&mesh, newNode2), criterionRAIS);
  PointInsertion pi3(&mesh, SimplicesNode(&mesh, newNode3), criterionRAIS);
  EdgeCollapse   ec3(&mesh, SimplicesNode(&mesh, newNode3), SimplicesNode(&mesh, newNode2), criterionRAIS);

  gmds::ISimplexMeshIOService ioService(&mesh);
  gmds::VTKWriter vtkWriter(&ioService);
  vtkWriter.setCellOptions(gmds::N|gmds::R);
  vtkWriter.write("eight_Tet_edge_collapse00.vtk");
}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, DISABLED_test00_simple_cube_minus_sphere)
{
  SimplexMesh simplexMesh = SimplexMesh();
  gmds::ISimplexMeshIOService ioService(&simplexMesh);

  std::string dir(TEST_SAMPLES_DIR);
  //std::string vtk_file = dir+"/cube_minus_sphere_42399.vtk";
  std::string vtk_file = dir+"/cube_minus_sphere_14857.vtk";
  //std::string vtk_file = dir+"/Cube.vtk";
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.read(vtk_file);

  //Modify the structure//
  simplexMesh.buildRobustLayerMesh(2);

  gmds::ISimplexMeshIOService ioServiceWriter(&simplexMesh);
  gmds::VTKWriter vtkWriter(&ioServiceWriter);
  vtkWriter.setCellOptions(gmds::R|gmds::N);
  vtkWriter.write("cube_minus_sphereWriteTest003.vtk");

  //TODO marquer les tetra qui sont en bord de la sphere

}
/*----------------------------------------------------------------------------*/
