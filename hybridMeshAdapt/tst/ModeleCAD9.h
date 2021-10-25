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
TEST(SimplexMeshTestClass, test_point_insertion_on_modelCAD9)
{/*
  TInt border = std::numeric_limits<TInt>::min();
  SimplexMesh simplexMesh = SimplexMesh();
  std::string dir(TEST_SAMPLES_DIR);
  std::string vtk_mesh = dir+"/ModeleCAD9/ModeleCAD9.vtk";

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

  std::string vtk_node = dir+"/ModeleCAD9/ModeleCAD9_hexa_generatedPoints_140458.vtk";
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

      if(!alreadyAdd)
      {
        if((*BND_CURVE_COLOR_NODES)[idx] != 0) {BND_CURVE_COLOR->set(node, (*BND_CURVE_COLOR_NODES)[idx]);}
        else if((*BND_SURFACE_COLOR_NODES)[idx] != 0) {BND_SURFACE_COLOR->set(node, (*BND_SURFACE_COLOR_NODES)[idx]);}

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
  vtkWriterDI.write("ModeleCAD9_DI_125489.vtk");


  std::cout << "INSERTION FINISH" << std::endl;
  std::cout << "edgesRemove start" << std::endl;
  simplexMesh.edgesRemove(nodesAdded);
  std::cout << "edgesRemove end" << std::endl;

  gmds::VTKWriter vtkWriterER(&ioService);
  vtkWriterER.setCellOptions(gmds::N|gmds::R);
  vtkWriterER.setDataOptions(gmds::N|gmds::R);
  vtkWriterER.write("ModeleCAD9_ER_125489.vtk");


  //for(auto & edge : edges)
  //{
  //  edge.push_back(0);
  //  edge.push_back(0);
  //}
  //for(unsigned int idx = 0 ; idx < nodes.size() ; idx++)
  //{
  //  TInt node = nodes[idx];
  //  if(node != border)
  //  {
  //    for(auto & edge : edges)
  //    {
  //      if(edge[2] == 0)
  //      {
  //          if(edge[0] == idx){edge[0] = node; edge[2] = 1;}
  //      }
  //        if(edge[3] == 0)
  //      {
  //          if(edge[1] == idx){edge[1] = node; edge[3] = 1;}
  //      }
  //    }
  //  }
  //}

  //std::vector<std::vector<TInt>> edgesGlobal(edges.size());
  //for(unsigned int idx = 0 ; idx < edges.size() ; idx++)
  //{
  //  std::vector<TInt> edge = edges[idx];
  //  edgesGlobal[idx].push_back(edge[0]);
  //  edgesGlobal[idx].push_back(edge[1]);
  //}
  //std::cout << "buildEdges start !" << std::endl;
  //simplexMesh.buildEdges(edgesGlobal, nodesAdded);
  //std::cout << "buildEdges done !" << std::endl;


  //gmds::VTKWriter vtkWriterEB(&ioService);
  //vtkWriterEB.setCellOptions(gmds::N|gmds::R);
  //vtkWriterEB.setDataOptions(gmds::N|gmds::R);
  //vtkWriterEB.write("ModeleCAD9_EB_15701_SCV.vtk");
  //vtkWriterEB.write("ModeleCAD9_EB_111016_SCV.vtk");
  */
}

TEST(SimplexMeshTestClass, test_hex_exctration_on_modelCAD9)
{
  TInt border = std::numeric_limits<TInt>::min();
  //Read the mesh file
  SimplexMesh simplexMesh = SimplexMesh();
  std::string dir(TEST_SAMPLES_DIR);
  std::string vtk_mesh = dir+"/ModeleCAD9/ModeleCAD9_ER_125489.vtk";

  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.setDataOptions(gmds::N);
  vtkReader.read(vtk_mesh);
  simplexMesh.buildAdjInfoGlobal();
  simplexMesh.buildSimplexHull();


  //Read the node added node file
  std::string vtk_node = dir+"/ModeleCAD9/ModeleCAD9_hexa_generatedPoints_140458.vtk";
  SimplexMesh simplexNodes = SimplexMesh();
  gmds::ISimplexMeshIOService ioServiceNodes(&simplexNodes);
  gmds::VTKReader vtkReaderNodes(&ioServiceNodes);
  vtkReaderNodes.setCellOptions(gmds::R|gmds::N|gmds::F);
  vtkReaderNodes.setDataOptions(gmds::N);
  vtkReaderNodes.read(vtk_node);

  std::vector<std::vector<TInt>> edges{};
  std::vector<std::vector<TInt>> nodesHex = simplexNodes.getHexadronData();
  const gmds::BitVector& triIdx = simplexNodes.getBitVectorTri();
  const gmds::BitVector& nodeMeshIds = simplexMesh.getBitVectorNodes();
  const gmds::BitVector& nodesIds = simplexNodes.getBitVectorNodes();
  gmds::BitVector nodesAdded(simplexMesh.nodesCapacity());
  const gmds::BitVector& tetMeshIds = simplexMesh.getBitVectorTet();
  std::vector<TInt> nodesToMesh(nodeMeshIds.capacity(), -1);

  /*
  std::cout << "START LOADING nodesToMesh" << std::endl;
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
  std::ofstream myfile;
  myfile.open ("markedNodeModeleCAD9_EB.txt");
  for(auto const idx :nodesToMesh)
  {
    myfile << idx << std::endl;
  }
  myfile.close();
  std::cout << "END LOADING nodesToMesh" << std::endl;
  */
  std::ifstream input;
  std::string filename = "markedNodeModeleCAD9_EB.txt";
  input.open(filename.c_str());
  std::string line;
  while(std::getline(input, line))
  {
    TInt idx = std::stod(line);
    nodesToMesh.push_back(idx);
  }

  for(auto const idx : nodesToMesh)
  {
    if(idx != -1)
    {
      nodeAdded.assign(idx);
    }
  }

  for(unsigned int tri = 1 ; tri < triIdx.capacity() ; tri++)
  {
    if(triIdx[tri] != 0)
    {
      std::vector<TInt> nodes = SimplicesTriangle(&simplexNodes,tri).getNodes();
      if(!(nodesToMesh[nodes[0]] == -1 || nodesToMesh[nodes[1]] == -1))
      std::vector<TInt> edge{nodesToMesh[nodes[0]], nodesToMesh[nodes[1]]};
      edges.push_back(edge);
    }
  }

  std::cout << "buildEdges start !" << std::endl;
  simplexMesh.buildEdges(edges, nodesAdded);
  std::cout << "buildEdges done !" << std::endl;













  std::vector<std::vector<TSimplexID>> hexesNodes{};
  gmds::BitVector markedTet(simplexMesh.getBitVectorTet().capacity());
  double hexBuiltCpt = 0.0;
  double remainTet   = 0.0;

  for(auto const h : nodesHex)
  {
    TInt n0 = nodesToMesh[h[0]]; TInt n1 = nodesToMesh[h[1]]; TInt n2 = nodesToMesh[h[2]]; TInt n3 = nodesToMesh[h[3]];
    TInt n4 = nodesToMesh[h[4]]; TInt n5 = nodesToMesh[h[5]]; TInt n6 = nodesToMesh[h[6]]; TInt n7 = nodesToMesh[h[7]];
    if(n0 == -1 || n1 == -1 || n2 == -1 || n3 == -1 || n4 == -1 || n5 == -1 || n6 == -1 || n7 == -1)
    {
      continue;
    }
    std::vector<TInt> hexeNodes{n0, n1, n2, n3, n4, n5, n6, n7};
    if(simplexMesh.buildFace(hexeNodes, nodesAdded))
    {
      std::vector<TSimplexID> simplices = simplexMesh.hex2tet(hexeNodes);
      for(auto const simplex : simplices){markedTet.assign(simplex);}
      if(simplices.size() != 0)
      {
        hexBuiltCpt = hexBuiltCpt + 1.0;
        hexesNodes.push_back(hexeNodes);
      }
      else
      {
        std::cout << "HEX NOT BUILT" << std::endl;
      }
    }
    else
    {
      std::cout << "simplexMesh.buildFace(hexeNodes, nodesAdded) == false" << std::endl;
    }
  }

  ////////////////// Result
  for(unsigned int tet = 0 ; tet < tetMeshIds.capacity() ; tet++)
  {
    if(tetMeshIds[tet] != 0)
    {
      if(markedTet[tet] == 0)
      {
        remainTet = remainTet + 1;
      }
    }
  }
  std::cout << "tet remaining % -> " << remainTet / (double)tetMeshIds.size() << std::endl;
  std::cout << "hex build % -> " << hexBuiltCpt / (double)hexesNodes.size() << std::endl;
  std::cout << "hexBuiltCpt -> " << hexBuiltCpt << std::endl;
  std::cout << "hexesNodes.size() --> " << hexesNodes.size() << std::endl;
  //////////////////
  simplexMesh.setHexadronData(hexesNodes);
  simplexMesh.setMarkedTet(markedTet);
  gmds::Variable<int>* TET_COLOR = simplexMesh.newVariable<int, SimplicesCell>("TET_COLOR");
  TET_COLOR->setValuesTo(1000);

  std::cout << "adding hex ..." << std::endl;
  gmds::VTKWriter vtkWriterHT(&ioService);
  vtkWriterHT.setCellOptions(gmds::N|gmds::R);
  vtkWriterHT.setDataOptions(gmds::N|gmds::R);
  //vtkWriterHT.write("ModeleCAD5_HEX_TET_FF3.vtk");
  vtkWriterHT.write("ModeleCAD9_HEX_TET_FF.vtk");
}
