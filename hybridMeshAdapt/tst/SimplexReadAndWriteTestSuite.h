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

/*TEST(SimplexMeshTestClass, test_point_insertion_on_modelCAD5)
{
  TInt border = std::numeric_limits<TInt>::min();
  SimplexMesh simplexMesh = SimplexMesh();
  std::string dir(TEST_SAMPLES_DIR);
  //std::string vtk_mesh = dir+"/ModeleCAD5/ModeleCAD5.vtk";
  std::string vtk_mesh = dir+"/ModeleCAD5Bis/ModeleCAD5Bis.vtk";

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

  //std::string vtk_node = dir+"/ModeleCAD5/ModeleCAD5_edge_generatedPoints_15701.vtk";
  std::string vtk_node = dir+"/ModeleCAD5Bis/ModeleCAD5_edge_generatedPoints_111016.vtk";
  //std::string vtk_node = dir+"/ModeleCAD5/ModeleCAD5_edge_generatedPoints_109316.vtk";

  SimplexMesh simplexNodes = SimplexMesh();
  gmds::ISimplexMeshIOService ioServiceNodes(&simplexNodes);
  gmds::VTKReader vtkReaderNodes(&ioServiceNodes);
  vtkReaderNodes.setCellOptions(gmds::R|gmds::N|gmds::F);
  vtkReaderNodes.setDataOptions(gmds::N);
  vtkReaderNodes.read(vtk_node);
  Variable<int>* BND_CURVE_COLOR_NODES   = simplexNodes.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR_NODES = simplexNodes.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");

  const gmds::BitVector& nodesToAddIds = simplexNodes.getBitVectorNodes();

  CriterionRAIS criterionRAIS(new VolumeCriterion());
  gmds::BitVector nodesAdded(simplexMesh.nodesCapacity());

  std::vector<std::vector<TInt>> edges{};
  const gmds::BitVector triIdx = simplexNodes.getBitVectorTri();
  for(unsigned int tri = 1 ; tri < triIdx.capacity() ; tri++)
  {
    if(triIdx[tri] != 0)
    {
      std::vector<TInt> nodes = SimplicesTriangle(&simplexNodes,tri).getNodes();
      std::vector<TInt> edge{nodes[0], nodes[1]};
      edges.push_back(edge);
    }
  }



  std::vector<TInt> nodes{};
  for(unsigned int idx = 0 ; idx < nodesToAddIds.capacity() ; idx++)
  {
    std::cout << std::endl;
    std::cout << std::endl;
    if(nodesToAddIds[idx] != 0)
    {
      std::cout << "Idx --> " << idx << std::endl;
      const gmds::BitVector & nodesIds = simplexMesh.getBitVectorNodes();
      math::Point point = SimplicesNode(&simplexNodes, idx).getCoords();

      bool alreadyAdd = false;
      std::vector<TSimplexID> tetraContenaingPt{};
      std::cout << "Before addNode" << std::endl;
      std::cout << "coords : " << point << std::endl;
      TInt node = simplexMesh.addNodeAndcheck(point, tetraContenaingPt, alreadyAdd);
      std::cout << "After addNode" << std::endl;
      nodes.push_back(node);

      if(idx == 49768 || idx == 62862)
      {
        simplexMesh.deleteNode(node);
        nodes.pop_back();
        nodes.push_back(border);
        continue;
      }
      if(!alreadyAdd)
      {
        if((*BND_CURVE_COLOR_NODES)[idx] != 0) {BND_CURVE_COLOR->set(node, (*BND_CURVE_COLOR_NODES)[idx]);}
        else if((*BND_SURFACE_COLOR_NODES)[idx] != 0) {BND_SURFACE_COLOR->set(node, (*BND_SURFACE_COLOR_NODES)[idx]);}

        //forcage de node a cause d'une mauvaise labelisation dans fram3D pour le model M5 et les 8396 points
        //if(node == 10506)
        //{
        //  BND_SURFACE_COLOR->set(node, 9);
        //}
        //else if(node == 8448)
        //{
        //  BND_SURFACE_COLOR->set(node, 5);
        //}

        simplexMesh.getVariable<Eigen::Matrix3d, SimplicesNode>("metric")->value(node) = m;
        bool status = false;
        std::vector<TSimplexID> deletedSimplex{};

        std::cout << "node -> " << node << " | " << (*BND_CURVE_COLOR_NODES)[idx] << " | " << (*BND_SURFACE_COLOR_NODES)[idx] << std::endl;
        DelaunayPointInsertion DI(&simplexMesh, SimplicesNode(&simplexMesh, node), criterionRAIS, tetraContenaingPt, status, nodesAdded, deletedSimplex);
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

  gmds::VTKWriter vtkWriterDI(&ioService);
  vtkWriterDI.setCellOptions(gmds::N|gmds::F|gmds::R);
  vtkWriterDI.setDataOptions(gmds::N|gmds::F|gmds::R);
  //vtkWriterDI.write("ModeleCAD5_DI_15701.vtk");
  vtkWriterDI.write("ModeleCAD5_DI_111016.vtk");


  std::cout << "INSERTION FINISH" << std::endl;
  std::cout << "edgesRemove start" << std::endl;
  simplexMesh.edgesRemove(nodesAdded);
  std::cout << "edgesRemove end" << std::endl;

  gmds::VTKWriter vtkWriterER(&ioService);
  vtkWriterER.setCellOptions(gmds::N|gmds::R);
  vtkWriterER.setDataOptions(gmds::N|gmds::R);
  //vtkWriterER.write("ModeleCAD5_ER_15701_SCV.vtk");
  vtkWriterER.write("ModeleCAD5_ER_111016_SCV.vtk");

  for(auto & edge : edges)
  {
    edge.push_back(0);
    edge.push_back(0);
  }
  for(unsigned int idx = 0 ; idx < nodes.size() ; idx++)
  {
    TInt node = nodes[idx];
    if(node != border)
    {
      for(auto & edge : edges)
      {
        if(edge[2] == 0)
        {
            if(edge[0] == idx){edge[0] = node; edge[2] = 1;}
        }
        if(edge[3] == 0)
        {
            if(edge[1] == idx){edge[1] = node; edge[3] = 1;}
        }
      }
    }
  }

  std::vector<std::vector<TInt>> edgesGlobal(edges.size());
  for(unsigned int idx = 0 ; idx < edges.size() ; idx++)
  {
    std::vector<TInt> edge = edges[idx];
    edgesGlobal[idx].push_back(edge[0]);
    edgesGlobal[idx].push_back(edge[1]);
  }
  std::cout << "buildEdges start !" << std::endl;
  simplexMesh.buildEdges(edgesGlobal, nodesAdded);
  std::cout << "buildEdges done !" << std::endl;


  gmds::VTKWriter vtkWriterEB(&ioService);
  vtkWriterEB.setCellOptions(gmds::N|gmds::R);
  vtkWriterEB.setDataOptions(gmds::N|gmds::R);
  //vtkWriterEB.write("ModeleCAD5_EB_15701_SCV.vtk");
  vtkWriterEB.write("ModeleCAD5_EB_111016_SCV.vtk");



}*/

/*TEST(SimplexMeshTestClass, test_edges_remove_on_modelCAD5)
{
  TInt border = std::numeric_limits<TInt>::min();
  //Read the mesh file
  SimplexMesh simplexMesh = SimplexMesh();
  std::string dir(TEST_SAMPLES_DIR);
  std::string vtk_mesh = dir+"/ModeleCAD5Bis/ModeleCAD5_DI_111016.vtk";

  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.setDataOptions(gmds::N);
  vtkReader.read(vtk_mesh);
  simplexMesh.buildAdjInfoGlobal();
  simplexMesh.buildSimplexHull();

  //Read the node added node file
  std::string vtk_node = dir+"/ModeleCAD5Bis/ModeleCAD5_edge_generatedPoints_111016.vtk";
  SimplexMesh simplexNodes = SimplexMesh();
  gmds::ISimplexMeshIOService ioServiceNodes(&simplexNodes);
  gmds::VTKReader vtkReaderNodes(&ioServiceNodes);
  vtkReaderNodes.setCellOptions(gmds::R|gmds::N|gmds::F);
  vtkReaderNodes.setDataOptions(gmds::N);
  vtkReaderNodes.read(vtk_node);
  Variable<int>* BND_CURVE_COLOR_NODES   = simplexNodes.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR_NODES = simplexNodes.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");

  const gmds::BitVector& nodeMeshIds = simplexMesh.getBitVectorNodes();
  const gmds::BitVector& nodesIds = simplexNodes.getBitVectorNodes();
  gmds::BitVector nodesAdded(simplexMesh.nodesCapacity());
  for(unsigned int idx = 0 ; idx < nodeMeshIds.capacity() ; idx++)
  {
    if(nodeMeshIds[idx] != 0)
    {
      math::Point pt = SimplicesNode(&simplexMesh, idx).getCoords();
      for(unsigned int idxNode = 0 ; idxNode < nodesIds.capacity() ; idxNode++)
      {
        if(nodesIds[idxNode] != 0)
        {
          math::Point ptNode = SimplicesNode(&simplexNodes, idxNode).getCoords();
          if(pt == ptNode)
          {
            nodesAdded.assign(idx);
          }

        }
      }
    }
  }

  std::cout << "Edge Remove have started" << std::endl;
  simplexMesh.edgesRemove(nodesAdded);
  std::cout << "Edge Remove Finish" << std::endl;

  gmds::VTKWriter vtkWriterER(&ioService);
  vtkWriterER.setCellOptions(gmds::N|gmds::R);
  vtkWriterER.setDataOptions(gmds::N|gmds::R);
  //vtkWriterER.write("ModeleCAD5_ER_15701_SCV.vtk");
  vtkWriterER.write("ModeleCAD5_ER_111016.vtk");
}*/

/*TEST(SimplexMeshTestClass, test_edge_builder_on_modelCAD5)
{
  TInt border = std::numeric_limits<TInt>::min();
  //Read the mesh file
  SimplexMesh simplexMesh = SimplexMesh();
  std::string dir(TEST_SAMPLES_DIR);
  std::string vtk_mesh = dir+"/ModeleCAD5Bis/ModeleCAD5_ER_111016.vtk";

  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.setDataOptions(gmds::N);
  vtkReader.read(vtk_mesh);
  simplexMesh.buildAdjInfoGlobal();
  simplexMesh.buildSimplexHull();

  //Read the node added node file
  std::string vtk_node = dir+"/ModeleCAD5Bis/ModeleCAD5_edge_generatedPoints_111016.vtk";
  SimplexMesh simplexNodes = SimplexMesh();
  gmds::ISimplexMeshIOService ioServiceNodes(&simplexNodes);
  gmds::VTKReader vtkReaderNodes(&ioServiceNodes);
  vtkReaderNodes.setCellOptions(gmds::R|gmds::N|gmds::F);
  vtkReaderNodes.setDataOptions(gmds::N);
  vtkReaderNodes.read(vtk_node);
  Variable<int>* BND_CURVE_COLOR_NODES   = simplexNodes.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR_NODES = simplexNodes.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");

  std::vector<std::vector<TInt>> edges{};
  const gmds::BitVector triIdx = simplexNodes.getBitVectorTri();
  const gmds::BitVector& nodeMeshIds = simplexMesh.getBitVectorNodes();
  const gmds::BitVector& nodesIds = simplexNodes.getBitVectorNodes();
  gmds::BitVector nodesAdded(simplexMesh.nodesCapacity());
  std::vector<TInt> nodesToMesh(nodesIds.size());

  for(unsigned int idx = 0 ; idx < nodeMeshIds.capacity() ; idx++)
  {
    if(nodeMeshIds[idx] != 0)
    {
      math::Point pt = SimplicesNode(&simplexMesh, idx).getCoords();
      for(unsigned int idxNode = 0 ; idxNode < nodesIds.capacity() ; idxNode++)
      {
        if(nodesIds[idxNode] != 0)
        {
          math::Point ptNode = SimplicesNode(&simplexNodes, idxNode).getCoords();
          if(pt == ptNode)
          {
            nodesToMesh[idxNode] = idx;
            nodesAdded.assign(idx);
          }
        }
      }
    }
  }

  for(unsigned int tri = 1 ; tri < triIdx.capacity() ; tri++)
  {
    if(triIdx[tri] != 0)
    {
      std::vector<TInt> nodes = SimplicesTriangle(&simplexNodes,tri).getNodes();
      std::vector<TInt> edge{nodesToMesh[nodes[0]], nodesToMesh[nodes[0]]};
      edges.push_back(edge);
    }
  }

  std::cout << "buildEdges start !" << std::endl;
  simplexMesh.buildEdges(edges, nodesAdded);
  std::cout << "buildEdges done !" << std::endl;


  gmds::VTKWriter vtkWriterEB(&ioService);
  vtkWriterEB.setCellOptions(gmds::N|gmds::R);
  vtkWriterEB.setDataOptions(gmds::N|gmds::R);
  //vtkWriterEB.write("ModeleCAD5_EB_15701_SCV.vtk");
  vtkWriterEB.write("ModeleCAD5_EB_111016.vtk");


}*/

TEST(SimplexMeshTestClass, test_isHexaBuild)
{
  TInt border = std::numeric_limits<TInt>::min();
  //Read the mesh file
  SimplexMesh simplexMesh = SimplexMesh();
  std::string dir(TEST_SAMPLES_DIR);
  std::string vtk_mesh = dir+"/ModeleCAD5Bis/ModeleCAD5_ER_111016.vtk";

  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.setDataOptions(gmds::N);
  std::cout << "reading ..." << std::endl;
  vtkReader.read(vtk_mesh);
  std::cout << "buildAdjInfoGlobal ..." << std::endl;
  simplexMesh.buildAdjInfoGlobal();
  //simplexMesh.buildSimplexHull();
  std::cout << "simplexMesh.isHexaEdgeBuild(faces) ..." << std::endl;
  std::vector<TInt> face0{41182, 80289, 80286, 91520}; std::vector<TInt> face3{80289, 80286, 80287, 89154};
  std::vector<TInt> face1{91520, 80286, 80287, 73998}; std::vector<TInt> face4{89154, 80287, 73998, 89142};
  std::vector<TInt> face2{41182, 91520, 73998, 89142}; std::vector<TInt> face5{41182, 80289, 89154, 89142};
  std::vector<std::vector<TInt>> faces{face0, face1, face2, face3, face4, face5};
  std::cout << "Is hexa built : " << simplexMesh.isHexaEdgeBuild(faces) << std::endl;

  std::vector<TInt> hexNodes{41182, 91520, 80286, 80289, 89142, 73998, 80287, 89154};
  std::vector<TSimplexID> simplices = simplexMesh.hex2tet(hexNodes);
  std::cout << "hexahedron composed of simplex " << std::endl;
  for(auto const simplex : simplices)
  {
    std::cout << simplex << std::endl;
  }
  std::cout << std::endl;

  simplexMesh.deleteAllSimplicesBut(simplices);
  gmds::VTKWriter vtkWriterEB(&ioService);
  vtkWriterEB.setCellOptions(gmds::N|gmds::R);
  vtkWriterEB.setDataOptions(gmds::N|gmds::R);
  vtkWriterEB.write("ModeleCAD5_Test_extraction_hex2tet.vtk");
}
