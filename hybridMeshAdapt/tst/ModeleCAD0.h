/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <fstream>
#include <bitset>
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
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
TEST(SimplexMeshTestClass, Add_Surface_Node_DI)
{
  CriterionRAIS criterionRAIS(new VolumeCriterion());
  TInt border = std::numeric_limits<TInt>::min();


  SimplexMesh simplexMeshNode = SimplexMesh();
  math::Point point = math::Point(-67.6934, 8, 5);
  simplexMeshNode.addNode(point);
  simplexMeshNode.addTetraedre(0,0,0,0);
  gmds::ISimplexMeshIOService ioServiceNode(&simplexMeshNode);
  gmds::VTKWriter vtkWriterNode(&ioServiceNode);
  vtkWriterNode.setCellOptions(gmds::N|gmds::F|gmds::R);
  vtkWriterNode.write("Surface_Node.vtk");

  SimplexMesh simplexMesh = SimplexMesh();
  std::string dir(TEST_SAMPLES_DIR);
  std::string vtk_mesh = dir+"/ModeleCAD5/ModeleCAD5.vtk";

  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.setDataOptions(gmds::N);
  vtkReader.read(vtk_mesh);
  simplexMesh.buildAdjInfoGlobal();
  simplexMesh.initializeEdgeStructure();
  simplexMesh.buildSimplexHull();

  //return;
  Octree oc(&simplexMesh, 50);
  simplexMesh.setOctree(&oc);
  Variable<int>* BND_VERTEX_COLOR  = simplexMesh.getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
  Variable<int>* BND_CURVE_COLOR   = simplexMesh.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR = simplexMesh.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");

  //Modification of the structure
  //////////////////////////////////////////////////////////////////////////////
  Variable<Eigen::Matrix3d>* var = simplexMesh.newVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  Eigen::Matrix3d m =  Eigen::MatrixXd::Identity(3, 3);
  m <<  1.0 , 0.0, .0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
  var->setValuesTo(m);

  bool alreadyAdd = false;
  std::vector<TSimplexID> tetraContenaingPt{};
  TInt node = simplexMesh.addNodeAndcheck(point, tetraContenaingPt, alreadyAdd);
  BND_SURFACE_COLOR->set(node, 10);

  std::vector<TSimplexID> cavity{};
  std::vector<TSimplexID> markedSimplex{};
  bool status = false;
  std::vector<TSimplexID> deletedSimplex{};
  const std::multimap<TInt, TInt>& facesAlreadyBuilt{};
  gmds::BitVector nodesAdded(simplexMesh.nodesCapacity());

  DelaunayPointInsertion DI(&simplexMesh, SimplicesNode(&simplexMesh, node), criterionRAIS, cavity, status, nodesAdded, deletedSimplex, facesAlreadyBuilt, markedSimplex);
  std::cout << "Status -> " << status << std::endl;

  gmds::VTKWriter vtkWriterDI(&ioService);
  vtkWriterDI.setCellOptions(gmds::N|gmds::F|gmds::R);
  vtkWriterDI.setDataOptions(gmds::N|gmds::F|gmds::R);
  vtkWriterDI.write("Surface_Node_Added_ModeleCAD5.vtk");
}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, Compute_Octree_On_Model)
{
  CriterionRAIS criterionRAIS(new VolumeCriterion());
  TInt border = std::numeric_limits<TInt>::min();

  /*************************************M5**********************************/
  SimplexMesh simplexMesh0 = SimplexMesh();
  std::string dir(TEST_SAMPLES_DIR);
  std::string vtk_mesh = dir+"/ModeleCAD5/ModeleCAD5.vtk";

  gmds::ISimplexMeshIOService ioService(&simplexMesh0);
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.setDataOptions(gmds::N);
  vtkReader.read(vtk_mesh);
  simplexMesh0.buildAdjInfoGlobal();
  simplexMesh0.initializeEdgeStructure();
  simplexMesh0.buildSimplexHull();

  Octree oc100(&simplexMesh0, 100, "ModeleCAD5_Octree_100");
  Octree oc50(&simplexMesh0, 50, "ModeleCAD5_Octree_50");

  /*************************************B31**********************************/
  SimplexMesh simplexMesh1 = SimplexMesh();
  vtk_mesh = dir+"/ModeleCAD4/ModeleCAD4.vtk";
  ioService = gmds::ISimplexMeshIOService(&simplexMesh1);
  vtkReader = gmds::VTKReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.setDataOptions(gmds::N);
  vtkReader.read(vtk_mesh);
  simplexMesh1.buildAdjInfoGlobal();
  simplexMesh1.initializeEdgeStructure();
  simplexMesh1.buildSimplexHull();

  oc100 = Octree(&simplexMesh1, 100, "B31_Octree_100");
  oc50 = Octree(&simplexMesh1, 50, "B31_Octree_50");

}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, Add_Surface_Node_PI)
{
  CriterionRAIS criterionRAIS(new VolumeCriterion());
  TInt border = std::numeric_limits<TInt>::min();


  SimplexMesh simplexMeshNode = SimplexMesh();
  math::Point point = math::Point(0.68, 0.066,1.0);
  simplexMeshNode.addNode(point);
  simplexMeshNode.addTetraedre(0,0,0,0);
  gmds::ISimplexMeshIOService ioServiceNode(&simplexMeshNode);
  gmds::VTKWriter vtkWriterNode(&ioServiceNode);
  vtkWriterNode.setCellOptions(gmds::N|gmds::F|gmds::R);
  vtkWriterNode.write("Surface_Node.vtk");

  SimplexMesh simplexMesh = SimplexMesh();
  std::string dir(TEST_SAMPLES_DIR);
  std::string vtk_mesh = dir+"/simpleCubeV2/CUBE.vtk";

  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.setDataOptions(gmds::N);
  vtkReader.read(vtk_mesh);
  simplexMesh.buildAdjInfoGlobal();
  simplexMesh.initializeEdgeStructure();
  simplexMesh.buildSimplexHull();

  //return;
  Octree oc(&simplexMesh, 50);
  simplexMesh.setOctree(&oc);
  Variable<int>* BND_VERTEX_COLOR  = simplexMesh.getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
  Variable<int>* BND_CURVE_COLOR   = simplexMesh.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR = simplexMesh.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");

  //Modification of the structure
  //////////////////////////////////////////////////////////////////////////////
  Variable<Eigen::Matrix3d>* var = simplexMesh.newVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  Eigen::Matrix3d m =  Eigen::MatrixXd::Identity(3, 3);
  m <<  1.0, 0.0, .0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
  var->setValuesTo(m);

  bool alreadyAdd = false;
  std::vector<TSimplexID> tetraContenaingPt{};
  TInt node = simplexMesh.addNodeAndcheck(point, tetraContenaingPt, alreadyAdd);
  BND_SURFACE_COLOR->set(node, 1);

  std::vector<TSimplexID> cavity{1322, 11724, 1317, 9007, 3299, 1331};
  std::vector<TSimplexID> markedSimplex{};
  bool status = false;
  std::vector<TSimplexID> deletedSimplex{};
  const std::multimap<TInt, TInt>& facesAlreadyBuilt{};
  gmds::BitVector nodesAdded(simplexMesh.nodesCapacity());

  DelaunayPointInsertion DI(&simplexMesh, SimplicesNode(&simplexMesh, node), criterionRAIS, cavity, status, nodesAdded, deletedSimplex, facesAlreadyBuilt, markedSimplex, false);
  std::cout << "Status -> " << status << std::endl;

  gmds::VTKWriter vtkWriterDI(&ioService);
  vtkWriterDI.setCellOptions(gmds::N|gmds::F|gmds::R);
  vtkWriterDI.setDataOptions(gmds::N|gmds::F|gmds::R);
  vtkWriterDI.write("Surface_Node_Added.vtk");
}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, Add_Curve_Node_PI)
{
  CriterionRAIS criterionRAIS(new VolumeCriterion());
  TInt border = std::numeric_limits<TInt>::min();


  SimplexMesh simplexMeshNode = SimplexMesh();
  math::Point point = math::Point(0.75, 0.0, 1.0);
  simplexMeshNode.addNode(point);
  simplexMeshNode.addTetraedre(0,0,0,0);
  gmds::ISimplexMeshIOService ioServiceNode(&simplexMeshNode);
  gmds::VTKWriter vtkWriterNode(&ioServiceNode);
  vtkWriterNode.setCellOptions(gmds::N|gmds::F|gmds::R);
  vtkWriterNode.write("Curve_Node.vtk");

  SimplexMesh simplexMesh = SimplexMesh();
  std::string dir(TEST_SAMPLES_DIR);
  std::string vtk_mesh = dir+"/simpleCubeV2/CUBE.vtk";

  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.setDataOptions(gmds::N);
  vtkReader.read(vtk_mesh);
  simplexMesh.buildAdjInfoGlobal();
  simplexMesh.initializeEdgeStructure();
  simplexMesh.buildSimplexHull();

  //return;
  Octree oc(&simplexMesh, 50);
  simplexMesh.setOctree(&oc);
  Variable<int>* BND_VERTEX_COLOR  = simplexMesh.getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
  Variable<int>* BND_CURVE_COLOR   = simplexMesh.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR = simplexMesh.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");

  //Modification of the structure
  //////////////////////////////////////////////////////////////////////////////
  Variable<Eigen::Matrix3d>* var = simplexMesh.newVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  Eigen::Matrix3d m =  Eigen::MatrixXd::Identity(3, 3);
  m <<  1.0, 0.0, .0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
  var->setValuesTo(m);

  bool alreadyAdd = false;
  std::vector<TSimplexID> tetraContenaingPt{};
  TInt node = simplexMesh.addNodeAndcheck(point, tetraContenaingPt, alreadyAdd);
  BND_CURVE_COLOR->set(node, 12);

  std::vector<TSimplexID> cavity{1322, 9007, 3299, 1331, 547, 3300};
  //std::vector<TSimplexID> cavity{};
  std::vector<TSimplexID> markedSimplex{};
  bool status = false;
  std::vector<TSimplexID> deletedSimplex{};
  const std::multimap<TInt, TInt>& facesAlreadyBuilt{};
  gmds::BitVector nodesAdded(simplexMesh.nodesCapacity());

  DelaunayPointInsertion DI(&simplexMesh, SimplicesNode(&simplexMesh, node), criterionRAIS, cavity, status, nodesAdded, deletedSimplex, facesAlreadyBuilt, markedSimplex, false);
  std::cout << "Status -> " << status << std::endl;
  gmds::VTKWriter vtkWriterDI(&ioService);
  vtkWriterDI.setCellOptions(gmds::N|gmds::F|gmds::R);
  vtkWriterDI.setDataOptions(gmds::N|gmds::F|gmds::R);
  vtkWriterDI.write("Curve_Node_Added.vtk");
}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, Add_Curve_Node_PI_2)
{
  CriterionRAIS criterionRAIS(new VolumeCriterion());
  TInt border = std::numeric_limits<TInt>::min();


  SimplexMesh simplexMeshNode = SimplexMesh();
  math::Point point = math::Point(0.95, 0.0, 1.0);
  simplexMeshNode.addNode(point);
  simplexMeshNode.addTetraedre(0,0,0,0);
  gmds::ISimplexMeshIOService ioServiceNode(&simplexMeshNode);
  gmds::VTKWriter vtkWriterNode(&ioServiceNode);
  vtkWriterNode.setCellOptions(gmds::N|gmds::F|gmds::R);
  vtkWriterNode.write("Curve_Node_2.vtk");

  SimplexMesh simplexMesh = SimplexMesh();
  std::string dir(TEST_SAMPLES_DIR);
  std::string vtk_mesh = dir+"/simpleCubeV2/CUBE.vtk";

  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.setDataOptions(gmds::N);
  vtkReader.read(vtk_mesh);
  simplexMesh.buildAdjInfoGlobal();
  simplexMesh.initializeEdgeStructure();
  simplexMesh.buildSimplexHull();

  //return;
  Octree oc(&simplexMesh, 50);
  simplexMesh.setOctree(&oc);
  Variable<int>* BND_VERTEX_COLOR  = simplexMesh.getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
  Variable<int>* BND_CURVE_COLOR   = simplexMesh.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR = simplexMesh.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");

  //Modification of the structure
  //////////////////////////////////////////////////////////////////////////////
  Variable<Eigen::Matrix3d>* var = simplexMesh.newVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  Eigen::Matrix3d m =  Eigen::MatrixXd::Identity(3, 3);
  m <<  1.0, 0.0, .0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
  var->setValuesTo(m);

  bool alreadyAdd = false;
  std::vector<TSimplexID> tetraContenaingPt{};
  TInt node = simplexMesh.addNodeAndcheck(point, tetraContenaingPt, alreadyAdd);
  BND_CURVE_COLOR->set(node, 12);

  std::vector<TSimplexID> cavity{842,833,837};
  std::vector<TSimplexID> markedSimplex{};
  bool status = false;
  std::vector<TSimplexID> deletedSimplex{};
  const std::multimap<TInt, TInt>& facesAlreadyBuilt{};
  gmds::BitVector nodesAdded(simplexMesh.nodesCapacity());

  DelaunayPointInsertion DI(&simplexMesh, SimplicesNode(&simplexMesh, node), criterionRAIS, cavity, status, nodesAdded, deletedSimplex, facesAlreadyBuilt, markedSimplex, false);
  std::cout << "Status -> " << status << std::endl;
  gmds::VTKWriter vtkWriterDI(&ioService);
  vtkWriterDI.setCellOptions(gmds::N|gmds::F|gmds::R);
  vtkWriterDI.setDataOptions(gmds::N|gmds::F|gmds::R);
  vtkWriterDI.write("Curve_Node_Added_2.vtk");
}
/*----------------------------------------------------------------------------*/
TEST(SimplexMeshTestClass, test_point_insertion_on_modelCAD0)
{
  /*TInt border = std::numeric_limits<TInt>::min();
  SimplexMesh simplexMesh = SimplexMesh();
  std::string dir(TEST_SAMPLES_DIR);
  std::string vtk_mesh = dir+"/ModeleCAD0/ModeleCAD0.vtk";

  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.setDataOptions(gmds::N);
  vtkReader.read(vtk_mesh);
  simplexMesh.buildAdjInfoGlobal();
  simplexMesh.buildSimplexHull();

  //return;
  Octree oc(&simplexMesh, 50);
  simplexMesh.setOctree(&oc);
  Variable<int>* BND_VERTEX_COLOR  = simplexMesh.getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
  Variable<int>* BND_CURVE_COLOR   = simplexMesh.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR = simplexMesh.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");

  //Modification of the structure
  //////////////////////////////////////////////////////////////////////////////
  Variable<Eigen::Matrix3d>* var = simplexMesh.newVariable<Eigen::Matrix3d, SimplicesNode>("metric");
  Eigen::Matrix3d m =  Eigen::MatrixXd::Identity(3, 3);
  m <<  1.0, 0.0, .0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
  var->setValuesTo(m);
  //////////////////////////////////////////////////////////////////////////////

  std::string vtk_node = dir+"/ModeleCAD0/ModeleCAD0_hexa_generatedPoints_372.vtk";
  SimplexMesh simplexNodes = SimplexMesh();
  gmds::ISimplexMeshIOService ioServiceNodes(&simplexNodes);
  gmds::VTKReader vtkReaderNodes(&ioServiceNodes);
  vtkReaderNodes.setCellOptions(gmds::R|gmds::N|gmds::F);
  vtkReaderNodes.setDataOptions(gmds::N);
  vtkReaderNodes.read(vtk_node);
  Variable<int>* BND_CURVE_COLOR_NODES   = simplexNodes.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR_NODES = simplexNodes.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");

  const gmds::BitVector& nodesToAddIds = simplexNodes.getBitVectorNodes();
  const gmds::BitVector& nodePresentInMesh = simplexMesh.getBitVectorNodes();
  std::cout << "INITIAL NODE SIZE IN MESH --> " << nodePresentInMesh.size() << std::endl;

  CriterionRAIS criterionRAIS(new VolumeCriterion());
  gmds::BitVector nodesAdded(simplexMesh.nodesCapacity());
  std::vector<TInt> nodes(nodesToAddIds.capacity(), -1);

  for(unsigned int idx = 0 ; idx < nodesToAddIds.capacity() ; idx++)
  {
    if(nodesToAddIds[idx] != 0)
    {
      const gmds::BitVector & nodesIds = simplexMesh.getBitVectorNodes();
      math::Point point = SimplicesNode(&simplexNodes, idx).getCoords();

      bool alreadyAdd = false;
      std::vector<TSimplexID> tetraContenaingPt{};
      TInt node = simplexMesh.addNodeAndcheck(point, tetraContenaingPt, alreadyAdd);
      nodes[idx] = node;

      if(!alreadyAdd)
      {
        if((*BND_CURVE_COLOR_NODES)[idx] != 0) {BND_CURVE_COLOR->set(node, (*BND_CURVE_COLOR_NODES)[idx]);}
        else if((*BND_SURFACE_COLOR_NODES)[idx] != 0) {BND_SURFACE_COLOR->set(node, (*BND_SURFACE_COLOR_NODES)[idx]);}

        simplexMesh.getVariable<Eigen::Matrix3d, SimplicesNode>("metric")->value(node) = m;
        bool status = false;
        std::vector<TSimplexID> deletedSimplex{};
        const std::multimap<TInt, std::pair<TInt, TInt>> facesAlreadyBuilt{};

        DelaunayPointInsertion DI(&simplexMesh, SimplicesNode(&simplexMesh, node), criterionRAIS, tetraContenaingPt, status, nodesAdded, deletedSimplex, facesAlreadyBuilt);
        if(status)
        {
          if(nodesAdded.capacity() != nodesIds.capacity())
          {
            nodesAdded.resize(nodesIds.capacity());
          }
          nodesAdded.assign(node);
        }
        else
        {
          simplexMesh.deleteNode(node);
          nodes.pop_back();
          nodes.push_back(border);
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
  }

  //std::cout << "INSERTION FINISH..." << std::endl;
  //std::cout << "NODE SIZE IN MESH AFTER DI --> " << nodePresentInMesh.size() << std::endl;
  //std::cout << "inserted node --> " << nodesAdded.size() << " | ~ " << (double)nodesAdded.size() / 372.0 << std::endl;
  gmds::VTKWriter vtkWriterDI(&ioService);
  vtkWriterDI.setCellOptions(gmds::N|gmds::F|gmds::R);
  vtkWriterDI.setDataOptions(gmds::N|gmds::F|gmds::R);
  vtkWriterDI.write("ModeleCAD0_DI_372.vtk");


  //std::cout << "edgesRemove start" << std::endl;
  std::clock_t start;
  double duration;
  start = std::clock();
  std::vector<TInt> deletedNodes{};
  unsigned int tmp = 0;
  for(;;)
  {
    unsigned int edgesRemoved = simplexMesh.edgesRemove(nodesAdded, deletedNodes);
    if(edgesRemoved == tmp || edgesRemoved == 0)
    {
      break;
    }
    tmp = edgesRemoved;
  }
  gmds::VTKWriter vtkWriterER(&ioService);
  vtkWriterER.setCellOptions(gmds::N|gmds::R);
  vtkWriterER.setDataOptions(gmds::N|gmds::R);
  vtkWriterER.write("ModeleCAD0_ER_372.vtk");

  unsigned int nodeSize = 0;
  unsigned int nodeReinsertedSize = 0;
  gmds::BitVector nBV(nodesAdded.capacity());
  for(auto const deletedNode : deletedNodes)
  {
    if(nodesAdded[deletedNode] == 1 && nBV[deletedNode] == 0)
    {
      nBV.assign(deletedNode);
      nodeSize++;
      bool alreadyAdd = false;
      std::vector<TSimplexID> tetraContenaingPt{};
      std::vector<TInt> deletedNodes{};
      const std::multimap<TInt, std::pair<TInt, TInt>> facesAlreadyBuilt{};
      bool status = false;
      PointInsertion(&simplexMesh, SimplicesNode(&simplexMesh, deletedNode), criterionRAIS, status, tetraContenaingPt, nodesAdded, deletedNodes, facesAlreadyBuilt);
      if(status)
      {
        nodeReinsertedSize++;
      }
    }
  }
  duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
  std::cout << "Edge removed duration -> " << duration << std::endl;

  //std::cout << "nBV.size() --> " << nBV.size() << std::endl;
  //std::cout << "nodeReinsertedSize --> " << nodeReinsertedSize << std::endl;
  const gmds::BitVector& nodesIds = simplexMesh.getBitVectorNodes();
  //std::cout << "nodesIds.capacity() --> " << nodesIds.capacity() << std::endl;
  std::cout << "EDGES REMOVAL FINISH..." << std::endl;

  gmds::VTKWriter vtkWriterDI2(&ioService);
  vtkWriterDI2.setCellOptions(gmds::N|gmds::R);
  vtkWriterDI2.setDataOptions(gmds::N|gmds::R);
  vtkWriterDI2.write("ModeleCAD0_DI_BIS_372.vtk");

  ///////////////////////////EDGE BUILDER START HERE///////////////////////////
  std::multimap<TInt, TInt> edges{};
  const gmds::BitVector triIdx = simplexNodes.getBitVectorTri();
  for(unsigned int tri = 1 ; tri < triIdx.capacity() ; tri++)
  {
    if(triIdx[tri] != 0)
    {
      std::vector<TInt> e  = SimplicesTriangle(&simplexNodes,tri).getNodes();
      if(!(nodes[e[0]] == -1 || nodes[e[1]] == -1))
      {
        std::vector<TInt> oe{std::min(nodes[e[0]], nodes[e[1]]), std::max(nodes[e[0]], nodes[e[1]])};
        std::pair<std::multimap<TInt, TInt>::iterator, std::multimap<TInt, TInt>::iterator> p;
        p = edges.equal_range(oe.front());
        bool edgesAlreadyAdded = false;
        for(std::multimap<TInt, TInt>::iterator it = p.first ; it != p.second ; it++)
        {
          if(it->second == oe.back())
          {
            edgesAlreadyAdded = true;
            break;
          }
        }

        if(!edgesAlreadyAdded)
        {
          edges.insert(std::pair<TInt, TInt>(oe.front(), oe.back()));
        }
      }
      else
      {
        std::cout << "e[0,1] --> " << e[0] << " | " << e[1] << std::endl;
      }
    }
  }


  std::vector<std::vector<TSimplexID>> hexesNodes{};
  std::vector<std::vector<TInt>> nodesHex = simplexNodes.getHexadronData();
  gmds::BitVector markedTet;

  double hexBuiltCpt = 0;
  unsigned int faceBuildCpt = 0;
  unsigned int iter = 0;
  duration;
  start = std::clock();
  for(;;)
  {
    std::cout << "ITERATION -> " << iter++ << std::endl;
    hexBuiltCpt = 0;
    hexesNodes.clear();
    unsigned int tmp = 0;
    for(;;)
    {
      unsigned int edgeBuild = simplexMesh.buildEdges(edges, nodesAdded);
      if(edgeBuild == tmp || edgeBuild == 0)
      {
        break;
      }
      tmp = edgeBuild;
    }
    //std::cout << "BUILD EDGE FINISH..." << std::endl;

    /////////////////////HEXA'S FACES BUILDER START HERE /////////////////////////
    std::multimap<TInt, std::pair<TInt, TInt>> facesAlreadyBuilt{};
    unsigned int faceBuiltTmp = 0;
    for(auto const h : nodesHex)
    {
      TInt n0 = nodes[h[0]]; TInt n1 = nodes[h[1]]; TInt n2 = nodes[h[2]]; TInt n3 = nodes[h[3]];
      TInt n4 = nodes[h[4]]; TInt n5 = nodes[h[5]]; TInt n6 = nodes[h[6]]; TInt n7 = nodes[h[7]];
      std::vector<TInt> hexeNodes{n0, n1, n2, n3, n4, n5, n6, n7};
      if(simplexMesh.buildFace(hexeNodes, nodesAdded, facesAlreadyBuilt))
      {
        faceBuiltTmp++;
      }
    }
    markedTet = simplexMesh.getBitVectorTet();

    for(auto const h : nodesHex)
    {
      TInt n0 = nodes[h[0]]; TInt n1 = nodes[h[1]]; TInt n2 = nodes[h[2]]; TInt n3 = nodes[h[3]];
      TInt n4 = nodes[h[4]]; TInt n5 = nodes[h[5]]; TInt n6 = nodes[h[6]]; TInt n7 = nodes[h[7]];
      std::vector<TInt> hexeNodes{n0, n1, n2, n3, n4, n5, n6, n7};
      std::vector<TSimplexID> simplices = simplexMesh.hex2tet(hexeNodes);
      if(simplices.size() != 0)
      {
        for(auto const simplex : simplices){ markedTet.unselect(simplex);}
        hexBuiltCpt = hexBuiltCpt + 1.0;
        hexesNodes.push_back(hexeNodes);
      }
    }
    if(faceBuiltTmp == faceBuildCpt || faceBuiltTmp == 0)
    {
      break;
    }
    faceBuildCpt = faceBuiltTmp;
  }

  duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
  std::cout << "Edge BUILD + FaceBuild + hex2tet duration -> " << duration << std::endl;

  std::cout << "hex build % -> " << hexBuiltCpt / (double)nodesHex.size() << std::endl;
  std::cout << "hexBuiltCpt -> " << hexBuiltCpt << std::endl;
  std::cout << "hexesNodes.size() --> " << nodesHex.size() << std::endl;
  //////////////////

  //ASSERT_EQ(hexBuiltCpt, 197);

  simplexMesh.setHexadronData(hexesNodes);
  simplexMesh.setMarkedTet(markedTet);
  std::cout << "adding hex ..." << std::endl;

  gmds::VTKWriter vtkWriterHT(&ioService);
  vtkWriterHT.setCellOptions(gmds::N|gmds::R);
  vtkWriterHT.setDataOptions(gmds::N|gmds::R);
  vtkWriterHT.write("ModeleCAD0_HEX.vtk");*/
}
