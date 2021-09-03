/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <fstream>
#include <bitset>
/*----------------------------------------------------------------------------*/
#include <gmds/hybridMeshAdapt/PointSmoothing.h>
#include <gmds/hybridMeshAdapt/PointInsertion.h>
#include <gmds/hybridMeshAdapt/EdgeCollapse.h>
#include <gmds/hybridMeshAdapt/ICriterion.h>
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
#include <gmds/hybridMeshAdapt/ISimplexMeshIOService.h>
#include <gmds/hybridMeshAdapt/EdgeInsertion.h>
#include <gmds/hybridMeshAdapt/DelaunayPointInsertion.h>
#include <gmds/hybridMeshAdapt/Octree.h>
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
#include <unit_test_config.h>
/*----------------------------------------------------------------------------*/
//#include <chrono>
//#include "sys/types.h"
//#include "sys/sysinfo.h"
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace hybrid;
using namespace operators;
using namespace simplicesNode;
using namespace simplicesTriangle;
using namespace simplicesCell;
/*----------------------------------------------------------------------------*/
TEST(SimplexReadAndWriteTestClass, test_write_simple_triangle)
{
  gmds::hybrid::SimplexMesh simplexMesh;

  simplexMesh.addNode(math::Point(0.0, 0.0, 0.0)); //0
  simplexMesh.addNode(math::Point(1.0, 0.0, 0.0)); //1
  simplexMesh.addNode(math::Point(0.0, 0.0, 1.0)); //2

  simplexMesh.addTriangle(0, 1, 2);

  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKWriter vtkWriter(&ioService);
  vtkWriter.setCellOptions(gmds::N|gmds::F);
  vtkWriter.write("simple_triangle.vtk");
}
/*----------------------------------------------------------------------------*/
TEST(SimplexReadAndWriteTestClass, test)
{
  gmds::hybrid::SimplexMesh simplexMesh;

  simplexMesh.addNode(math::Point(-10.6112, -5, -10.0744)); //0


  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKWriter vtkWriter(&ioService);
  vtkWriter.setCellOptions(gmds::N|gmds::F);
  vtkWriter.write("simple_.vtk");
}
/*----------------------------------------------------------------------------*/
TEST(SimplexReadAndWriteTestClass, test_write_simple_triangles)
{
  gmds::hybrid::SimplexMesh simplexMesh;

  simplexMesh.addNode(math::Point(0.0, 0.0, 0.0)); // Node 0
  simplexMesh.addNode(math::Point(1.0, 0.0, 0.0)); // Node 1
  simplexMesh.addNode(math::Point(0.0, 0.0, 1.0)); // Node 2
  simplexMesh.addNode(math::Point(0.0, 1.0, 0.0)); // Node 3
  simplexMesh.addNode(math::Point(0.0, -1.0, 0.0)); // Node 4
  simplexMesh.addNode(math::Point(0.0, 0.0, -1.0)); // Node 5
  simplexMesh.addNode(math::Point(1.0, 0.0, 1.0)); // Node 6

  simplexMesh.addTetraedre(0, 1 , 2, 3); // Tetra 0
  simplexMesh.addTetraedre(0, 1 , 2, 4); // Tetra 1
  simplexMesh.addTetraedre(0, 1 , 4, 5); // Tetra 2
  simplexMesh.addTetraedre(0, 1 , 3, 5); // Tetra 3

  // the triangle must be oriented with a face of tetraedron
  simplexMesh.addTriangle(0, 1, 2); //-1
  simplexMesh.addTriangle(0, 1, 5); //-2
  simplexMesh.addTriangle(0, 1, 3); //-
  simplexMesh.addTriangle(1, 4, 2); //-4
  simplexMesh.addTriangle(1, 2, 3); //-5

  simplexMesh.addTetraedre(1, 3 , 2, 6); // Tetra 4*/
  Octree oc(&simplexMesh, 1);
  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKWriter vtkWriter(&ioService);
  vtkWriter.setCellOptions(gmds::N|gmds::F|gmds::R);
  vtkWriter.write("simple_triangles.vtk");
}
/*----------------------------------------------------------------------------*/
TEST(MeshClass, test_volume)
{
  gmds::hybrid::SimplexMesh simplexMesh;
  simplexMesh.addNode(math::Point(0.0, 0.0, 0.0)); // Node 0
  simplexMesh.addNode(math::Point(2.0, 0.0, 0.0)); // Node 1
  simplexMesh.addNode(math::Point(0.0, 0.0, 2.0)); // Node 2
  simplexMesh.addNode(math::Point(0.0, 2.0, 0.0)); // Node 3

  simplexMesh.addTetraedre(0, 1 , 2, 3); // Tetra 0

  SimplicesCell cell = SimplicesCell(&simplexMesh, 0);
  double volumeCell = cell.getVolumeOfCell();
}
TEST(MeshClass, test_write_hull_of_simpleCube)
{
    /*gmds::hybrid::SimplexMesh simplexMesh;

    std::string dir(TEST_SAMPLES_DIR);
    gmds::ISimplexMeshIOService ioService(&simplexMesh);
    //std::string vtk_file = dir+"/simpleCube.vtk";
    std::string vtk_file = dir+"/ModeleEnS/ModeleEnS.vtk";
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::R|gmds::N);
    vtkReader.setDataOptions(gmds::N);
    vtkReader.read(vtk_file);
    simplexMesh.buildAdjInfoGlobal();

    simplexMesh.buildSimplexHull();

    gmds::VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::F);
    vtkWriter.setDataOptions(gmds::F);
    vtkWriter.write("ModeleEnS_Hull.vtk");*/
}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test_nodes_labeling)
{/*
  SimplexMesh simplexMesh = SimplexMesh();
  std::string dir(TEST_SAMPLES_DIR);
  //std::string vtk_file = dir+"/halfCylinder/HalfCylinder.vtk";
  //std::string vtk_file = dir+"/Dome/Dome.vtk";
  //std::string vtk_file = dir+"/ModeleEnS/ModeleEnS.vtk";
  //std::string vtk_file = dir+"/CylindreTordu/CylindreTordu.vtk";
  //std::string vtk_file = dir+"/ModeleCAD1/ModeleCAD1.vtk";
  std::string vtk_file = dir+"/ModeleCAD2/ModeleCAD2.vtk";

  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.setDataOptions(gmds::N);
  vtkReader.read(vtk_file);
  simplexMesh.buildAdjInfoGlobal();
  simplexMesh.buildSimplexHull();

  Octree oc(&simplexMesh, 50);
  simplexMesh.setOctree(&oc);

  std::vector<math::Point> points{};
  std::ifstream input;
  //std::string filename = dir+"/halfCylinder/HalfCylinder_points_2147.txt";
  //std::string filename = dir+"/Dome/Dome_generated_2922.txt";
  //std::string filename = dir+"/ModeleEnS/GeneratedPoints4466.txt";
  //std::string filename = dir+"/CylindreTordu/CylindreTordu_generated_5686.txt";
  //std::string filename = dir+"/ModeleCAD1/ModeleCAD1_generated_3033.txt"; <--- probleme avec les indices des surface de base
  std::string filename = dir+"/ModeleCAD2/ModeleCAD2_generated_5199.txt";

  input.open(filename.c_str());
  double x,y,z;
  char ch;
  std::string line;
  std::string delim = " ";
  while(std::getline(input, line))
  {
    std::string strX = line.substr(0, line.find(delim));
    line = line.substr(line.find(delim)+1);
    std::string strY = line.substr(0, line.find(delim));
    line = line.substr(line.find(delim)+1);
    std::string strZ = line.substr(0, line.find(delim));

    x = std::stod(strX);
    y = std::stod(strY);
    z = std::stod(strZ);
    math::Point point(x,y,z);
    points.push_back(point);
  }

  //start points labeling
  SimplexMesh simplexMeshNodes = SimplexMesh();
  std::vector<int> labelPoints{};
  std::vector<int> topoIndex{};
  std::vector<TSimplexID> nearCell{};
  simplexMesh.pointsLabeling(points, labelPoints, topoIndex, nearCell);
  for(auto const & point : points)
  {
    TInt node = simplexMeshNodes.addNode(point);
  }

  gmds::Variable<int>* BND_NODES_TOPO = simplexMeshNodes.newVariable<int, SimplicesNode>("BND_NODES_TOPO");
  gmds::Variable<int>* BND_NODES_INDICES = simplexMeshNodes.newVariable<int, SimplicesNode>("BND_NODES_INDICES");
  if(labelPoints.size() == topoIndex.size())
  {
    const gmds::BitVector& nodeBitVector = simplexMeshNodes.getBitVectorNodes();
    unsigned int cpt = 0;
    for(unsigned int nodeIds = 0 ; nodeIds < nodeBitVector.capacity() ; nodeIds++)
    {
      if(nodeBitVector[nodeIds] == 1)
      {
        simplexMeshNodes.addTetraedre(nodeIds, nodeIds, nodeIds, nodeIds);
        BND_NODES_TOPO->set(nodeIds, labelPoints[cpt]);
        BND_NODES_INDICES->set(nodeIds, topoIndex[cpt]);

        cpt++;
      }
    }
  }
  else
  {
    std::cout << "labelPoints.size() != topoIndex.size()" << std::endl;
  }

  gmds::ISimplexMeshIOService ioServiceNodes(&simplexMeshNodes);
  gmds::VTKWriter vtkWriter(&ioServiceNodes);
  vtkWriter.setCellOptions(gmds::N|gmds::R);
  vtkWriter.setDataOptions(gmds::N|gmds::R);
  vtkWriter.write("ModeleCAD2_5199_NodesLabellingOCTREE.vtk");*/
}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test_DelaunayInsertion_With_Labelized_Nodes)
{
  /*SimplexMesh simplexMesh = SimplexMesh();
  std::string dir(TEST_SAMPLES_DIR);
  std::string vtk_file = dir+"/halfCylinder/HalfCylinder.vtk";
  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.setDataOptions(gmds::N);
  vtkReader.read(vtk_file);
  simplexMesh.buildAdjInfoGlobal();
  simplexMesh.buildSimplexHull();

  std::vector<math::Point> points{};
  std::ifstream input;
  std::string filename = dir+"/halfCylinder/HalfCylinder_points_6683.txt";
  input.open(filename.c_str());
  double x,y,z;
  char ch;
  std::string line;
  std::string delim = " ";
  while(std::getline(input, line))
  {
    std::string strX = line.substr(0, line.find(delim));
    line = line.substr(line.find(delim)+1);
    std::string strY = line.substr(0, line.find(delim));
    line = line.substr(line.find(delim)+1);
    std::string strZ = line.substr(0, line.find(delim));

    x = std::stod(strX);
    y = std::stod(strY);
    z = std::stod(strZ);
    math::Point point(x,y,z);
    points.push_back(point);
  }

  //start points labeling
  std::vector<int> labelPoints{};
  std::vector<int> topoIndex{};
  simplexMesh.pointsLabeling(points, labelPoints, topoIndex);
  gmds::Variable<int>* BND_NODES_TOPO = simplexMesh.newVariable<int, SimplicesNode>("BND_NODES_TOPO");
  gmds::Variable<int>* BND_NODES_INDICES = simplexMesh.newVariable<int, SimplicesNode>("BND_NODES_INDICES");
  if(labelPoints.size() != topoIndex.size())
  {
    std::cout << "labelPoints.size() != topoIndex.size()" << std::endl;
    return;
  }


  //Modification of the structure
  //////////////////////////////////////////////////////////////////////////////
  Variable<Eigen::Matrix3d>* var = simplexMesh.newVariable<Eigen::Matrix3d, SimplicesNode>("metric");
  Eigen::Matrix3d m =  Eigen::MatrixXd::Identity(3, 3);
  m <<  1.0, 0.0, .0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
  var->setValuesTo(m);
  CriterionRAIS criterionRAIS(new VolumeCriterion());
  for(unsigned int idxPt = 0 ; idxPt < points.size() ;idxPt++)
  {
    math::Point point = points[idxPt];
    bool alreadyAdd = false;
    std::vector<TSimplexID> tetraContenaingPt{};
    SimplicesNode::nodeNeighborInfo nodeInfo;
    TInt node = simplexMesh.addNodeAndcheck(point, tetraContenaingPt, alreadyAdd);
    BND_NODES_TOPO->set(node, labelPoints[idxPt]);
    BND_NODES_INDICES->set(node, topoIndex[idxPt]);

    if(!alreadyAdd)
    {
      simplexMesh.getVariable<Eigen::Matrix3d, SimplicesNode>("metric")->value(node) = m;
      DelaunayPointInsertion DI(&simplexMesh, SimplicesNode(&simplexMesh, node), criterionRAIS, tetraContenaingPt);
    }
  }
  //////////////////////////////////////////////////////////////////////////////

  gmds::VTKWriter vtkWriter(&ioService);
  vtkWriter.setCellOptions(gmds::N|gmds::F|gmds::R);
  vtkWriter.setDataOptions(gmds::N|gmds::F|gmds::R);
  vtkWriter.write("DelaunayPointInsertion_LabelizedNode_E_S_V_6683Nodes.vtk");*/
}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test_DelaunayInsertion_With_Labelized_Nodes_With_EdgeCollapse)
{
  SimplexMesh simplexMesh = SimplexMesh();
  std::string dir(TEST_SAMPLES_DIR);
  //std::string vtk_file = dir+"/halfCylinder/HalfCylinder.vtk";
  //std::string vtk_file = dir+"/Dome/Dome.vtk";
  //std::string vtk_file = dir+"/ModeleEnS/ModeleEnS.vtk";
  //std::string vtk_file = dir+"/CylindreTordu/CylindreTordu.vtk";
  //std::string vtk_file = dir+"/ModeleCAD1/ModeleCAD1.vtk";
  std::string vtk_file = dir+"/ModeleCAD2/ModeleCAD2.vtk";


  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.setDataOptions(gmds::N);
  vtkReader.read(vtk_file);
  simplexMesh.buildAdjInfoGlobal();
  simplexMesh.buildSimplexHull();

  Octree oc(&simplexMesh, 50);
  simplexMesh.setOctree(&oc);

  std::vector<math::Point> points{};
  std::ifstream input;
  //std::string filename = dir+"/halfCylinder/HalfCylinder_points_6683.txt";
  //std::string filename = dir+"/Dome/Dome_generated_2922.txt";
  //std::string filename = dir+"/ModeleEnS/GeneratedPoints4466.txt";
  //std::string filename = dir+"/CylindreTordu/CylindreTordu_generated_5686.txt";
  //std::string filename = dir+"/ModeleCAD1/ModeleCAD1_generated_3033.txt";
  std::string filename = dir+"/ModeleCAD2/ModeleCAD2_generated_5199.txt";

  input.open(filename.c_str());
  double x,y,z;
  char ch;
  std::string line;
  std::string delim = " ";
  while(std::getline(input, line))
  {
    std::string strX = line.substr(0, line.find(delim));
    line = line.substr(line.find(delim)+1);
    std::string strY = line.substr(0, line.find(delim));
    line = line.substr(line.find(delim)+1);
    std::string strZ = line.substr(0, line.find(delim));

    x = std::stod(strX);
    y = std::stod(strY);
    z = std::stod(strZ);
    math::Point point(x,y,z);
    points.push_back(point);
  }

  //start points labeling
  std::vector<int> labelPoints{};
  std::vector<int> topoIndex{};
  std::vector<TSimplexID> nearCells{};
  std::map<TSimplexID, std::vector<unsigned int>> indexMap{};
  simplexMesh.pointsLabeling(points, labelPoints, topoIndex, nearCells);
  for(unsigned int index = 0 ; index < nearCells.size() ; index++)
  {
    TSimplexID simplex = nearCells[index];
    indexMap[simplex].push_back(index);
  }


  gmds::Variable<int>* BND_NODES_TOPO = simplexMesh.newVariable<int, SimplicesNode>("BND_NODES_TOPO");
  gmds::Variable<int>* BND_NODES_INDICES = simplexMesh.newVariable<int, SimplicesNode>("BND_NODES_INDICES");
  if(labelPoints.size() != topoIndex.size())
  {
    std::cout << "labelPoints.size() != topoIndex.size()" << std::endl;
    return;
  }

  Variable<int>* BND_VERTEX_COLOR  = simplexMesh.getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
  Variable<int>* BND_CURVE_COLOR   = simplexMesh.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR = simplexMesh.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");

  //Modification of the structure
  //////////////////////////////////////////////////////////////////////////////
  Variable<Eigen::Matrix3d>* var = simplexMesh.newVariable<Eigen::Matrix3d, SimplicesNode>("metric");
  Eigen::Matrix3d m =  Eigen::MatrixXd::Identity(3, 3);
  m <<  1.0, 0.0, .0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
  var->setValuesTo(m);
  CriterionRAIS criterionRAIS(new VolumeCriterion());
  gmds::BitVector nodesAdded(simplexMesh.nodesCapacity());
  for(unsigned int idxPt = 0 ; idxPt < points.size() ; idxPt++)
  {
    const gmds::BitVector & triIds   = simplexMesh.getBitVectorTri();
    const gmds::BitVector & nodesIds = simplexMesh.getBitVectorNodes();

    math::Point point = points[idxPt];
    TSimplexID nearCell = nearCells[idxPt];
    std::vector<unsigned int> indices = indexMap[nearCell];


    bool alreadyAdd = false;
    std::vector<TSimplexID> tetraContenaingPt{};
    TInt node = simplexMesh.addNodeAndcheck(point, tetraContenaingPt, alreadyAdd, nearCell);

    if(labelPoints[idxPt]      == SimplexMesh::topo::CORNER)  {BND_VERTEX_COLOR->set(node, topoIndex[idxPt]);}
    else if(labelPoints[idxPt] == SimplexMesh::topo::RIDGE)   {BND_CURVE_COLOR->set(node, topoIndex[idxPt]);}
    else if(labelPoints[idxPt] == SimplexMesh::topo::SURFACE) {BND_SURFACE_COLOR->set(node, topoIndex[idxPt]);}

    if(!alreadyAdd)
    {
      simplexMesh.getVariable<Eigen::Matrix3d, SimplicesNode>("metric")->value(node) = m;
      if(labelPoints[idxPt] != SimplexMesh::topo::CORNER)
      {
        bool status = false;
        std::vector<TSimplexID> deletedSimplex{};
        DelaunayPointInsertion DI(&simplexMesh, SimplicesNode(&simplexMesh, node), criterionRAIS, tetraContenaingPt, status, nodesAdded, deletedSimplex);
        if(status)
        {
          TSimplexID nearCell = simplexMesh.getSimplexFromBase(node);
          for(auto const index : indices)
          {
            nearCells[index] = nearCell;
          }

          if(nodesAdded.capacity() != nodesIds.capacity())
          {
            nodesAdded.resize(nodesIds.capacity());
          }
          nodesAdded.assign(node);
        }
      }
    }
    else
    {
      if(nodesAdded.capacity() != nodesIds.capacity())
      {
        nodesAdded.resize(nodesIds.capacity());
      }
      nodesAdded.assign(node);
    }
  }

  gmds::VTKWriter vtkWriterDI(&ioService);
  vtkWriterDI.setCellOptions(gmds::N|gmds::F|gmds::R);
  vtkWriterDI.setDataOptions(gmds::N|gmds::F|gmds::R);
  vtkWriterDI.write("ModeleCAD2_DI_5199.vtk");
  //////////////////////////////////////////////////////////////////////////////
  //simplexMesh.fillBNDVariable();
  //////////////////////////////////////////////////////////////////////////////
  std::cout << "edgesRemove start" << std::endl;
  simplexMesh.edgesRemove(nodesAdded);
  std::cout << "edgesRemove end" << std::endl;

  gmds::VTKWriter vtkWriter(&ioService);
  vtkWriter.setCellOptions(gmds::N|gmds::F|gmds::R);
  vtkWriter.setDataOptions(gmds::N|gmds::F|gmds::R);
  //vtkWriter.write("Dome_ER_2922.vtk");
  //vtkWriter.write("HalfCylinder_ER_6683.vtk");
  //vtkWriter.write("ModeleEnS_ER_4466.vtk");
  //vtkWriter.write("CylindreTordu_ER_5686.vtk");
  //vtkWriter.write("ModeleCAD_ER_5686.vtk");
  vtkWriter.write("ModeleCAD2_ER_5199.vtk");
  std::cout << "write" << std::endl;
}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test02_point_insertion_on_two_tetra)
{
    /*SimplexMesh mesh = SimplexMesh();
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
    vtkWriter.write("eight_Tet_point_insertion02.vtk");*/
}
/******************************************************************************/
TEST(SimplexMeshTestClass, test00_edge_collapse_on_two_tetra)
{
  /*SimplexMesh mesh = SimplexMesh();
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
  vtkWriter.write("eight_Tet_edge_collapse00.vtk");*/
}
/*----------------------------------------------------------------------------*/
/*TEST(SimplexMeshTestClass, test00_simple_cube_minus_sphere)
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
  simplexMesh.buildRobustLayerMeshOrderedNode01(1);

  gmds::ISimplexMeshIOService ioServiceWriter(&simplexMesh);
  gmds::VTKWriter vtkWriter(&ioServiceWriter);
  vtkWriter.setCellOptions(gmds::R|gmds::N);
  vtkWriter.write("cube_minus_sphereWriteTest000.vtk");

  //TODO marquer les tetra qui sont en bord de la sphere

}*/
/*----------------------------------------------------------------------------*/
/*TEST(SimplexMeshTestClass, insertion_of_quadTet_in_TestMesh)
{
  SimplexMesh simplexMesh = SimplexMesh();
  gmds::ISimplexMeshIOService ioService(&simplexMesh);

  std::string dir(TEST_SAMPLES_DIR);
  std::string vtk_file = dir+"/simpleCube.vtk";
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.read(vtk_file);

  //CriterionRAIS criterionRAIS(new VolumeCriterion());
  //TInt newNode0 = simplexMesh.addNode(math::Point(0.8, 0.2, 0.8));
  //PointInsertion pi0(&simplexMesh, SimplicesNode(&simplexMesh, newNode0), criterionRAIS);

  //Modify the structure//
  std::vector<math::Point> pointCoords{};
  pointCoords.push_back(0.2 *math::Point(0.2, 0.2, 0.2));
  pointCoords.push_back(0.2 *math::Point(0.2, 0.2, 0.8));
  pointCoords.push_back(0.2 *math::Point(0.8, 0.2, 0.8));
  pointCoords.push_back(0.2 *math::Point(0.8, 0.2, 0.2));
  pointCoords.push_back(0.2 *math::Point(0.2, 0.8, 0.2));
  pointCoords.push_back(0.2 *math::Point(0.2, 0.8, 0.8));
  pointCoords.push_back(0.2 *math::Point(0.8, 0.8, 0.8));
  pointCoords.push_back(0.2 *math::Point(0.8, 0.8, 0.2));

  std::vector<TInt> nodesQuad0{};
  for(auto const & point : pointCoords)
  {
    TInt node = simplexMesh.addNode(point);
    nodesQuad0.push_back(node);
    CriterionRAIS criterionRAIS(new VolumeCriterion());
    PointInsertion pi(&simplexMesh, SimplicesNode(&simplexMesh, node), criterionRAIS);
  }

  std::vector<math::Point> pointCoords1{};
  pointCoords1.push_back(0.2 *math::Point(0.8, 0.2, 0.2));
  pointCoords1.pVariable<Eigen::Matrix3d>* var = simplexMesh.newVariable<Eigen::Matrix3d, SimplicesNode>("metric");
  Eigen::Matrix3d m =  Eigen::MatrixXd::Identity(3, 3);
  m <<  1.0, 0.0, .0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0;
  var->setValuesTo(m);ush_back(0.2 *math::Point(0.8, 0.2, 0.8));
  pointCoords1.push_back(0.2 *math::Point(1.4, 0.2, 0.8));
  pointCoords1.push_back(0.2 *math::Point(1.4, 0.2, 0.2));
  pointCoords1.push_back(0.2 *math::Point(0.8, 0.8, 0.2));
  pointCoords1.push_back(0.2 *math::Point(0.8, 0.8, 0.8));
  pointCoords1.push_back(0.2 *math::Point(1.4, 0.8, 0.8));
  pointCoords1.push_back(0.2 *math::Point(1.4, 0.8, 0.2));

  std::vector<TInt> nodesQuad1{};
  for(auto const & point : pointCoords1)
  {
    TInt node = simplexMesh.addNode(point);
    nodesQuad1.push_back(node);
    CriterionRAIS criterionRAIS(new VolumeCriterion());
    PointInsertion pi(&simplexMesh, SimplicesNode(&simplexMesh, node), criterionRAIS);
  }


  std::vector<math::Point> pointCoords2{};
  pointCoords2.push_back(0.2 *math::Point(0.2, 0.2, 0.8));
  pointCoords2.push_back(0.2 *math::Point(0.2, 0.2, 1.4));
  pointCoords2.push_back(0.2 *math::Point(0.8, 0.2, 1.4));
  pointCoords2.push_back(0.2 *math::Point(0.8, 0.2, 0.8));
  pointCoords2.push_back(0.2 *math::Point(0.2, 0.8, 0.8));
  pointCoords2.push_back(0.2 *math::Point(0.2, 0.8, 1.4));
  pointCoords2.push_back(0.2 *math::Point(0.8, 0.8, 1.4));
  pointCoords2.push_back(0.2 *math::Point(0.8, 0.8, 0.8));

  std::vector<TInt> nodesQuad2{};
  for(auto const & point : pointCoords2)
  {
    TInt node = simplexMesh.addNode(point);
    nodesQuad2.push_back(node);
    CriterionRAIS criterionRAIS(new VolumeCriterion());
    PointInsertion pi(&simplexMesh, SimplicesNode(&simplexMesh, node), criterionRAIS);
  }

  std::vector<math::Point> pointCoords3{};
  pointCoords3.push_back(0.2 *math::Point(0.8, 0.2, 0.8));
  pointCoords3.push_back(0.2 *math::Point(0.8, 0.2, 1.4));
  pointCoords3.push_back(0.2 *math::Point(1.4, 0.2, 1.4));
  pointCoords3.push_back(0.2 *math::Point(1.4, 0.2, 0.8));
  pointCoords3.push_back(0.2 *math::Point(0.8, 0.8, 0.8));
  pointCoords3.push_back(0.2 *math::Point(0.8, 0.8, 1.4));
  pointCoords3.push_back(0.2 *math::Point(1.4, 0.8, 1.4));
  pointCoords3.push_back(0.2 *math::Point(1.4, 0.8, 0.8));

  std::vector<TInt> nodesQuad3{};
  for(auto const & point : pointCoords3)
  {
    TInt node = simplexMesh.addNode(point);
    nodesQuad3.push_back(node);
    CriterionRAIS criterionRAIS(new VolumeCriterion());
    PointInsertion pi(&simplexMesh, SimplicesNode(&simplexMesh, node), criterionRAIS);
  }


  std::vector<math::Point> pointCoordsRotate{};
  pointCoordsRotate.push_back(math::Point(0.2, 0.2, 0.6));
  pointCoordsRotate.push_back(math::Point(0.5, 0.2, 0.8));
  pointCoordsRotate.push_back(math::Point(0.8, 0.2, 0.5));
  pointCoordsRotate.push_back(math::Point(0.5, 0.2, 0.3));
  pointCoordsRotate.push_back(math::Point(0.2, 0.8, 0.52));
  pointCoordsRotate.push_back(math::Point(0.5, 0.8, 0.8));
  pointCoordsRotate.push_back(math::Point(0.8, 0.8, 0.5));
  pointCoordsRotate.push_back(math::Point(0.5, 0.8, 0.2));

  std::vector<TInt> nodesQuadRotate{};
  for(auto const & point : pointCoordsRotate)
  {
    TInt node = simplexMesh.addNode(point);
    nodesQuadRotate.push_back(node);
    CriterionRAIS criterionRAIS(new VolumeCriterion());
    PointInsertion pi(&simplexMesh, SimplicesNode(&simplexMesh, node), criterionRAIS);
  }

  CriterionRAIS criterionRAIS(new VolumeCriterion());
  simplexMesh.buildHexaedre(nodesQuad0, criterionRAIS);
  simplexMesh.buildHexaedre(nodesQuad1, criterionRAIS);
  simplexMesh.buildHexaedre(nodesQuad2, criterionRAIS);
  simplexMesh.buildHexaedre(nodesQuad3, criterionRAIS);
  //simplexMesh.buildHexaedre(nodesQuadRotate, criterionRAIS);

  gmds::ISimplexMeshIOService ioServiceWriter(&simplexMesh);
  gmds::VTKWriter vtkWriter(&ioServiceWriter);
  vtkWriter.setCellOptions(gmds::R|gmds::N);
  vtkWriter.write("simpleCubeToRebuildPInsertion.vtk");

}*/
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test_insertion_edge_new_alfgo)
{
    /*SimplexMesh mesh = SimplexMesh();
    mesh.addNode(math::Point(0.0, 0.0, 0.0)); // Node 0
    mesh.addNode(math::Point(1.0, 0.0, 0.0)); // Node 1
    mesh.addNode(math::Point(0.0, 0.0, 1.0)); // Node 2
    mesh.addNode(math::Point(0.1, 0.2, 0.2)); // Node 3
    mesh.addNode(math::Point(0.5, -0.3, 0.5)); // Node 4

    mesh.addNode(math::Point(0.1, -0.2, 0.2)); // Node 5
    mesh.addNode(math::Point(0.2, -0.5, 0.6)); // Node 5

    mesh.addTetraedre(0, 1 , 2, 3); // Tetra 0
    mesh.addTetraedre(0, 1 , 2, 4); // Tetra 1
    mesh.addTetraedre(0, 2 , 4, 5); // Tetra 0
    mesh.addTetraedre(2, 4 , 5, 6); // Tetra 0


    CriterionRAIS criterionRAIS(new VolumeCriterion());
    SimplicesNode P    = SimplicesNode(&mesh, 3);
    SimplicesNode Pbis = SimplicesNode(&mesh, 6);
    EdgeInsertion(&mesh, P, Pbis, criterionRAIS);
    mesh.checkMesh();

    gmds::ISimplexMeshIOService ioService(&mesh);
    gmds::VTKWriter vtkWriter(&ioService);
    vtkWriter.setCellOptions(gmds::N|gmds::R);
    vtkWriter.write("insertion_Edge_test.vtk");*/
}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, insertion_of_quadTet_in_TestMesh)
{
  /*SimplexMesh simplexMesh = SimplexMesh();
  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  std::string dir(TEST_SAMPLES_DIR);
  std::string vtk_file = dir+"/simpleCube.vtk";
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.read(vtk_file);
  simplexMesh.buildAdjInfoGlobal();
  if(!simplexMesh.checkMesh())
  {
    std::cout << "simplexMesh.checkMesh()" << std::endl;
    return;
  }

  //Modify the structure//
  std::vector<std::vector<math::Point>> pointCoordsHexa{};
  std::vector<math::Point> pointCoords{};
  unsigned int layerX = 10;
  unsigned int layerY = 10;
  unsigned int layerZ = 10;
  double border       = 0.1;

  double sizeShapeX = 1.0;
  double sizeShapeY = 1.0;
  double sizeShapeZ = 1.0;

  std::vector<std::vector<TInt>> nodesQuads{};
  auto start = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_secondsAddnode;
  std::chrono::duration<double> elapsed_secondsDelaunay;
  TInt node;

  Variable<Eigen::Matrix3d>* var = simplexMesh.newVariable<Eigen::Matrix3d, SimplicesNode>("metric");
  Eigen::Matrix3d m =  Eigen::MatrixXd::Identity(3, 3);
  m <<  1.0, 0.0, .0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0;
  var->setValuesTo(m);
  for(unsigned int u = 0; u < layerX; u++)
  {
    for(unsigned int v = 0; v < layerY; v++)
    {
      for(unsigned int w = 0; w < layerZ; w++)
      {
        pointCoords.resize(8);

        double sizeX = (sizeShapeX - 2.0 *border) / (double)layerX;
        double sizeY = (sizeShapeY - 2.0 *border) / (double)layerY;
        double sizeZ = (sizeShapeZ - 2.0 *border) / (double)layerZ;

        pointCoords[0] = (math::Point(border + u * sizeX , border + v * sizeY, border + w * sizeZ));
        pointCoords[1] = (math::Point(border + (u + 1) * sizeX , border + v * sizeY, border + w * sizeZ));
        pointCoords[2] = (math::Point(border + (u + 1) * sizeX , border + v * sizeY, border + (w + 1) * sizeZ));
        pointCoords[3] = (math::Point(border + u * sizeX , border + v * sizeY, border + (w + 1) * sizeZ));

        pointCoords[4] = (math::Point(border + u * sizeX , border + (v + 1) * sizeY, border + w * sizeZ));
        pointCoords[5] = (math::Point(border + (u + 1) * sizeX , border + (v + 1) * sizeY, border + w * sizeZ));
        pointCoords[6] = (math::Point(border + (u + 1) * sizeX , border + (v + 1) * sizeY, border + (w + 1) * sizeZ));
        pointCoords[7] = (math::Point(border + u * sizeX , border + (v + 1) * sizeY, border + (w + 1) * sizeZ));

        std::vector<TInt> nodesQuad{};
        CriterionRAIS criterionRAIS(new VolumeCriterion());
        for(auto const & point : pointCoords)
        {
          bool alreadyAdd = false;
          std::vector<TSimplexID> tetraContenaingPt{};
          auto startAdd = std::chrono::steady_clock::now();
          node = simplexMesh.addNodeAndcheck(point, tetraContenaingPt, alreadyAdd);
          auto endAdd = std::chrono::steady_clock::now();
          elapsed_secondsAddnode += endAdd - startAdd;
          nodesQuad.push_back(node);

          if(!alreadyAdd)
          {
            simplexMesh.getVariable<Eigen::Matrix3d, SimplicesNode>("metric")->value(node) = m;
            auto startDelaunay = std::chrono::steady_clock::now();
            DelaunayPointInsertion DI(&simplexMesh, SimplicesNode(&simplexMesh, node), criterionRAIS, tetraContenaingPt);
            auto endDelaunay = std::chrono::steady_clock::now();
            elapsed_secondsDelaunay += endDelaunay - startDelaunay;
          }
        }
        nodesQuads.push_back(nodesQuad);
        //simplexMesh.buildHexaedre(nodesQuad, criterionRAIS);
        pointCoords.clear();
      }
    }
  }
  if(!simplexMesh.checkMesh())
  {
    std::cout << "simplexMesh.checkMesh()" << std::endl;
    return;
  }
  auto end = std::chrono::steady_clock::now();

  double edgePerformance = 0.0;
  double hexaPerformance = 0.0;
  simplexMesh.hexaBuildPerformance(nodesQuads, edgePerformance, hexaPerformance);
  std::cout << "edge performance = " << edgePerformance << "%"<< std::endl;
  std::cout << "hexa performance = " << hexaPerformance << "%"<< std::endl;

  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
  std::cout << "elapsed time: " << elapsed_secondsAddnode.count() << "s\n";
  std::cout << "elapsed time: " << elapsed_secondsDelaunay.count() << "s\n";

  gmds::ISimplexMeshIOService ioServiceWriter(&simplexMesh);
  gmds::VTKWriter vtkWriter(&ioServiceWriter);
  vtkWriter.setCellOptions(gmds::R|gmds::N);
  vtkWriter.write("test.vtk");*/
}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, DISABLED_insertion_nodes_in_mesh)
{
  /*test that would verify if the detection of the cavity initial by projection in the faces &
    the insertion by delaunay is valid or not
  */

    SimplexMesh simplexMesh = SimplexMesh();
    gmds::ISimplexMeshIOService ioService(&simplexMesh);
    std::string dir(TEST_SAMPLES_DIR);

    //extract point data of the test's file
    std::vector<math::Point> points{};
    std::ifstream input;
    std::string filename = dir+"/halfCylinder/HalfCylinder_points_522.txt";
    input.open(filename.c_str());

    double x,y,z;
    char ch;
    std::string line;
    std::string delim = " ";
    while(std::getline(input, line))
    {
      std::string strX = line.substr(0, line.find(delim));
      line = line.substr(line.find(delim)+1);
      std::string strY = line.substr(0, line.find(delim));
      line = line.substr(line.find(delim)+1);
      std::string strZ = line.substr(0, line.find(delim));

      x = std::stod(strX);
      y = std::stod(strY);
      z = std::stod(strZ);
      math::Point point(x,y,z);
      points.push_back(point);

    }

    //mesh that will be use for the detection and insertion by delaunay algorithm
    std::string vtk_file = dir+"/halfCylinder/HalfCylinder.vtk";
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::R|gmds::N);
    vtkReader.setDataOptions(gmds::N);
    vtkReader.read(vtk_file);
    simplexMesh.buildAdjInfoGlobal();


    //Modify the structure/
    Variable<Eigen::Matrix3d>* var = simplexMesh.newVariable<Eigen::Matrix3d, SimplicesNode>("metric");
    Eigen::Matrix3d m =  Eigen::MatrixXd::Identity(3, 3);
    m <<  1.0, 0.0, .0,
          0.0, 1.0, 0.0,
          0.0, 0.0, 1.0;
    var->setValuesTo(m);

    CriterionRAIS criterionRAIS(new VolumeCriterion());
    std::vector<TInt> nodes{};
    nodes.reserve(points.size());
    for(auto const & point : points)
    {
      bool alreadyAdd = false;
      std::vector<TSimplexID> tetraContenaingPt{};
      SimplicesNode::nodeNeighborInfo nodeInfo;
      TInt node = simplexMesh.addNodeAndcheck(point, tetraContenaingPt, alreadyAdd);
      nodes.push_back(node);
      //if(!alreadyAdd)
      {
        simplexMesh.getVariable<Eigen::Matrix3d, SimplicesNode>("metric")->value(node) = m;
        bool status = false;
        //DelaunayPointInsertion DI(&simplexMesh, SimplicesNode(&simplexMesh, node), criterionRAIS, tetraContenaingPt, status);
        if(!(simplexMesh.checkMeshLocal(SimplicesNode(&simplexMesh, node))))
        {
          std::cout << "PROBLEME with the node : " << node << std::endl;
          return;
        }
      }
    }



    std::string name = "insertion_nodes_in_mesh_TEST.vtk";
    gmds::ISimplexMeshIOService ioServiceWriter(&simplexMesh);
    gmds::VTKWriter vtkWriter(&ioServiceWriter);
    vtkWriter.setCellOptions(gmds::R|gmds::N);
    vtkWriter.write(name);
}
/*----------------------------------------------------------------------------*/
/*TEST(SimplexMeshTestClass, insertion_grid_Point_in_tet)
{
  SimplexMesh simplexMesh = SimplexMesh();
  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  std::string dir(TEST_SAMPLES_DIR);

  //read vtk file
  std::string vtk_file = dir+"/data_tet_points/tet.vtk";
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.read(vtk_file);
  simplexMesh.buildAdjInfoGlobal();


  //Modify the structure/
  Variable<Eigen::Matrix3d>* var = simplexMesh.newVariable<Eigen::Matrix3d, SimplicesNode>("metric");
  Eigen::Matrix3d m =  Eigen::MatrixXd::Identity(3, 3);
  m <<  1.0, 0.0, .0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0;
  var->setValuesTo(m);

  CriterionRAIS criterionRAIS(new VolumeCriterion());

  double layerX = 20;
  double layerY = 20;
  double layerZ = 20;

  double sizeShapeX = 2.5;
  double sizeShapeY = 5.0;
  double sizeShapeZ = 5.0;

  std::vector<std::vector<TInt>> nodesQuads{};
  std::vector<math::Point> pointCoords{};
  math::Point origineBox = math::Point(-5.0, -2.5, -2.5);

  std::vector<TInt> addedNode{};

  for(unsigned int u = 0; u < layerX; u++)
  {
    for(unsigned int v = 0; v < layerY; v++)
    {
      for(unsigned int w = 0; w < layerZ; w++)
      {
        pointCoords.resize(8, origineBox);

        double sizeX = (sizeShapeX ) / (double)layerX;
        double sizeY = (sizeShapeY ) / (double)layerY;
        double sizeZ = (sizeShapeZ ) / (double)layerZ;

        pointCoords[0] = origineBox + (math::Point(u * sizeX ,v * sizeY, w * sizeZ));
        pointCoords[1] = origineBox + (math::Point((u + 1) * sizeX ,v * sizeY, w * sizeZ));
        pointCoords[2] = origineBox + (math::Point((u + 1) * sizeX ,v * sizeY, (w + 1) * sizeZ));
        pointCoords[3] = origineBox + (math::Point(u * sizeX , v * sizeY, (w + 1) * sizeZ));

        pointCoords[4] = origineBox + (math::Point(u * sizeX , (v + 1) * sizeY, w * sizeZ));
        pointCoords[5] = origineBox + (math::Point((u + 1) * sizeX ,(v + 1) * sizeY, w * sizeZ));
        pointCoords[6] = origineBox + (math::Point((u + 1) * sizeX ,(v + 1) * sizeY, (w + 1) * sizeZ));
        pointCoords[7] = origineBox + (math::Point(u * sizeX , (v + 1) * sizeY, (w + 1) * sizeZ));

        std::vector<TInt> nodesQuad{};
        for(auto const & point : pointCoords)
        {
          bool alreadyAdd = false;
          std::vector<TSimplexID> tetraContenaingPt{};
          TInt node = simplexMesh.addNodeAndcheck(point, tetraContenaingPt, alreadyAdd);
          nodesQuad.push_back(node);
          addedNode.push_back(node);

          if(!alreadyAdd)
          {
            simplexMesh.getVariable<Eigen::Matrix3d, SimplicesNode>("metric")->value(node) = m;
            DelaunayPointInsertion DI(&simplexMesh, SimplicesNode(&simplexMesh, node), criterionRAIS, tetraContenaingPt);
          }
        }
        nodesQuads.push_back(nodesQuad);
        pointCoords.clear();
      }
    }
  }

  gmds::BitVector addedNodeBitVector(simplexMesh.nodesCapacity());
  for(auto const & node : addedNode)
  {
    addedNodeBitVector.assign(node);
  }

  //double edgePerformance = 0.0;
  //double hexaPerformance = 0.0;
  //std::cout << "Avant le forçage des arêtes : " << std::endl;
  //simplexMesh.hexaBuildPerformance(nodesQuads, edgePerformance, hexaPerformance);
  //std::cout << "edge performance = " << edgePerformance << "%"<< std::endl;
  //std::cout << "hexa performance = " << hexaPerformance << "%"<< std::endl;



  std::cout << "EDGE COLLAPSE STARTING  " << std::endl;
  simplexMesh.edgesRemove(addedNodeBitVector);
  std::cout << "EDGE COLLAPSE ENDING  " << std::endl;

  std::cout << "EDGE COLLAPSE STARTING  " << std::endl;
  simplexMesh.edgesRemove(addedNodeBitVector);
  std::cout << "EDGE COLLAPSE ENDING  " << std::endl;

  std::cout << "EDGE COLLAPSE STARTING  " << std::endl;
  simplexMesh.edgesRemove(addedNodeBitVector);
  std::cout << "EDGE COLLAPSE ENDING  " << std::endl;

  std::cout << "EDGE COLLAPSE STARTING  " << std::endl;
  simplexMesh.edgesRemove(addedNodeBitVector);
  std::cout << "EDGE COLLAPSE ENDING  " << std::endl;

  //edgePerformance = 0.0;
  //hexaPerformance = 0.0;
  //std::cout << "Avant le forçage des arêtes : " << std::endl;
  //simplexMesh.hexaBuildPerformance(nodesQuads, edgePerformance, hexaPerformance);
  //std::cout << "edge performance = " << edgePerformance << "%"<< std::endl;
  //std::cout << "hexa performance = " << hexaPerformance << "%"<< std::endl;

  gmds::ISimplexMeshIOService ioServiceWriter(&simplexMesh);
  gmds::VTKWriter vtkWriter(&ioServiceWriter);
  vtkWriter.setCellOptions(gmds::R|gmds::N);
  vtkWriter.write("CHOIX_2_newALGO.vtk");
}*/
/*----------------------------------------------------------------------------*/
/*TEST(SimplexMeshTestClass, insertion_of_quadTet_in_TestMesh2)
{
  SimplexMesh simplexMesh = SimplexMesh();
  gmds::ISimplexMeshIOService ioService(&simplexMesh);

  std::string dir(TEST_SAMPLES_DIR);
  std::string vtk_file = dir+"/simpleCube.vtk";
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.read(vtk_file);


  std::vector<math::Point> pointCoordsRotate{};
  pointCoordsRotate.push_back(math::Point(0.2, 0.2, 0.6));
  pointCoordsRotate.push_back(math::Point(0.5, 0.2, 0.8));
  pointCoordsRotate.push_back(math::Point(0.8, 0.2, 0.5));
  pointCoordsRotate.push_back(math::Point(0.5, 0.2, 0.3));
  pointCoordsRotate.push_back(math::Point(0.2, 0.8, 0.52));
  pointCoordsRotate.push_back(math::Point(0.5, 0.8, 0.8));
  pointCoordsRotate.push_back(math::Point(0.8, 0.8, 0.5));
  pointCoordsRotate.push_back(math::Point(0.5, 0.8, 0.2));

  std::vector<TInt> nodesQuadRotate{};
  for(auto const & point : pointCoordsRotate)
  {
    TInt node = simplexMesh.addNode(point);
    nodesQuadRotate.push_back(node);
    CriterionRAIS criterionRAIS(new VolumeCriterion());
    PointInsertion pi(&simplexMesh, SimplicesNode(&simplexMesh, node), criterionRAIS);
  }

  CriterionRAIS criterionRAIS(new VolumeCriterion());
  simplexMesh.buildHexaedre(nodesQuadRotate, criterionRAIS);

  gmds::ISimplexMeshIOService ioServiceWriter(&simplexMesh);
  gmds::VTKWriter vtkWriter(&ioServiceWriter);
  vtkWriter.setCellOptions(gmds::R|gmds::N);
  vtkWriter.write("simpleCubeToRebuildPInsertionRotationTest.vtk");

}*/
/*----------------------------------------------------------------------------*/
/*TEST(SimplexMeshTestClass, custom_SimplexMesh)
{
  SimplexMesh m;
  math::Point p0 = math::Point(0.0, -0.5, 0.0);
  math::Point p1 = math::Point(1.0, 0.3, 0.0);
  math::Point p2 = math::Point(0.0, 1.0, 0.0);
  math::Point p3 = math::Point(0.0, 0.0, 1.0);
  math::Point p4 = math::Point(0.2, 0.3, -1.0);
  math::Point p5 = math::Point(1.0, -0.5, 0.0);
  math::Point p6 = math::Point(1.0, -0.2, -0.5);
  math::Point p7 = math::Point(1.0, 0.2, -1.3);
  math::Point p8 = math::Point(1.0, 0.5, -0.7);
  math::Point p9 = math::Point(1.0, -0.6, -1.0);




  m.addNode(p0);
  m.addNode(p1);
  m.addNode(p2);
  m.addNode(p3);
  m.addNode(p4);
  m.addNode(p5);
  m.addNode(p6);
  m.addNode(p7);
  m.addNode(p8);
  m.addNode(p9);

  m.addTetraedre(0,1,2,3);
  m.addTetraedre(0,1,2,4);
  m.addTetraedre(0,1,4,5);
  m.addTetraedre(1,4,5,6);
  //m.addTetraedre(0,1,3,5);
  //m.addTetraedre(4,5,6,7);
  //m.addTetraedre(0,4,5,7);
  m.addTetraedre(1,4,6,7);
  //m.addTetraedre(1,4,7,8);
  //m.addTetraedre(0,5,7,9);

  //math::Point p_inter = math::Point(0.2, 0.18, -0.2);
  //TInt pi = m.addNode(p_inter);
  //std::vector<TSimplexID> initCavity{0,1};
  //CriterionRAIS criterionRAIS(new VolumeCriterion());
  //PointInsertion(&m, SimplicesNode(&m, pi), criterionRAIS, initCavity);
  //m.deleteTetra(0);

  CriterionRAIS criterionRAIS(new VolumeCriterion());
  EdgeInsertion::EdgeInsertion(&m, SimplicesNode(&m, 3), SimplicesNode(&m, 7), criterionRAIS);
  //m.addTetraedre(3,4,5,7);


  gmds::ISimplexMeshIOService ioService(&m);
  gmds::VTKWriter vtkWriter(&ioService);
  vtkWriter.setCellOptions(gmds::N|gmds::R);
  //vtkWriter.setDataOptions(gmds::N|gmds::R); // a voir quand je rajouterai les valeurs aux sommets triangle tetra cell ....
  vtkWriter.write("customSimplexMeshInit.vtk");

}*/
/*----------------------------------------------------------------------------*/
/*
if(nodeBitVector[nodeCav] == 0  && (*BND_SURFACE_COLOR)[nodeCav] == surfaceNbr && nodeCav != data.node)
{
  {
    setCustomCavity.insert(nodeCav);
  }
}
*/
