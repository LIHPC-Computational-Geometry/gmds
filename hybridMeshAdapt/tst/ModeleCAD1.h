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
TEST(SimplexMeshTestClass, test_point_insertion_on_modelCAD1)
{
  TInt border = std::numeric_limits<TInt>::min();
  SimplexMesh simplexMesh = SimplexMesh();
  std::string dir(TEST_SAMPLES_DIR);
  std::string vtk_mesh = dir+"/ModeleCAD1/ModeleCAD1.vtk";

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

  std::string vtk_node = dir+"/ModeleCAD1/ModeleCAD1_hexa_generatedPoints_646.vtk";
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

  std::cout << "INSERTION FINISH..." << std::endl;
  std::cout << "NODE SIZE IN MESH AFTER DI --> " << nodePresentInMesh.size() << std::endl;
  std::cout << "inserted node --> " << nodesAdded.size() << " | ~ " << (double)nodesAdded.size() / 646.0 << std::endl;
  gmds::VTKWriter vtkWriterDI(&ioService);
  vtkWriterDI.setCellOptions(gmds::N|gmds::F|gmds::R);
  vtkWriterDI.setDataOptions(gmds::N|gmds::F|gmds::R);
  vtkWriterDI.write("ModeleCAD1_DI_646.vtk");


  std::cout << "edgesRemove start" << std::endl;
  std::vector<TInt> deletedNodes{};
  simplexMesh.edgesRemove(nodesAdded, deletedNodes);
  gmds::VTKWriter vtkWriterER(&ioService);
  vtkWriterER.setCellOptions(gmds::N|gmds::R);
  vtkWriterER.setDataOptions(gmds::N|gmds::R);
  vtkWriterER.write("ModeleCAD1_ER_646.vtk");

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

  //std::cout << "nBV.size() --> " << nBV.size() << std::endl;
  //std::cout << "nodeReinsertedSize --> " << nodeReinsertedSize << std::endl;
  const gmds::BitVector& nodesIds = simplexMesh.getBitVectorNodes();
  //std::cout << "nodesIds.capacity() --> " << nodesIds.capacity() << std::endl;
  std::cout << "EDGES REMOVAL FINISH..." << std::endl;

  gmds::VTKWriter vtkWriterDI2(&ioService);
  vtkWriterDI2.setCellOptions(gmds::N|gmds::R);
  vtkWriterDI2.setDataOptions(gmds::N|gmds::R);
  vtkWriterDI2.write("ModeleCAD1_DI_BIS_646.vtk");

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
  gmds::BitVector markedTet(simplexMesh.getBitVectorTet().capacity());
  double hexBuiltCpt = 0;
  unsigned int iterMax = 4;
  for(unsigned int iter = 0 ; iter < iterMax ; iter++)
  {
    markedTet.clear();
    hexBuiltCpt = 0;
    hexesNodes.clear();
    simplexMesh.buildEdges(edges, nodesAdded);
    simplexMesh.buildEdges(edges, nodesAdded);
    simplexMesh.buildEdges(edges, nodesAdded);
    simplexMesh.buildEdges(edges, nodesAdded);
    std::cout << "BUILD EDGE FINISH..." << std::endl;

    /////////////////////HEXA'S FACES BUILDER START HERE /////////////////////////
    std::multimap<TInt, std::pair<TInt, TInt>> facesAlreadyBuilt{};
    for(auto const h : nodesHex)
    {
      TInt n0 = nodes[h[0]]; TInt n1 = nodes[h[1]]; TInt n2 = nodes[h[2]]; TInt n3 = nodes[h[3]];
      TInt n4 = nodes[h[4]]; TInt n5 = nodes[h[5]]; TInt n6 = nodes[h[6]]; TInt n7 = nodes[h[7]];
      std::vector<TInt> hexeNodes{n0, n1, n2, n3, n4, n5, n6, n7};
      //std::cout << "n0,1,2,3,4,5,6,7 -> " << n0 << " | " << n1 << " | " << n2 << " | " << n3 << " | " << n4 << " | " << n5 << " | " << n6 << " | " << n7 << " | " << std::endl;
      simplexMesh.buildFace(hexeNodes, nodesAdded, facesAlreadyBuilt);
    }

    for(auto const h : nodesHex)
    {
      TInt n0 = nodes[h[0]]; TInt n1 = nodes[h[1]]; TInt n2 = nodes[h[2]]; TInt n3 = nodes[h[3]];
      TInt n4 = nodes[h[4]]; TInt n5 = nodes[h[5]]; TInt n6 = nodes[h[6]]; TInt n7 = nodes[h[7]];
      std::vector<TInt> hexeNodes{n0, n1, n2, n3, n4, n5, n6, n7};
      std::vector<TSimplexID> simplices = simplexMesh.hex2tet(hexeNodes);
      if(simplices.size() != 0)
      {
        for(auto const simplex : simplices){markedTet.assign(simplex);}
        hexBuiltCpt = hexBuiltCpt + 1.0;
        hexesNodes.push_back(hexeNodes);
      }
      else
      {
        std::cout << "HEX NOT BUILT" << std::endl;
      }
    }
  }

  std::cout << "hex build % -> " << hexBuiltCpt / (double)nodesHex.size() << std::endl;
  std::cout << "hexBuiltCpt -> " << hexBuiltCpt << std::endl;
  std::cout << "hexesNodes.size() --> " << nodesHex.size() << std::endl;
  //////////////////
  ASSERT_EQ(hexBuiltCpt, 421);
  return;


  simplexMesh.setHexadronData(hexesNodes);
  simplexMesh.setMarkedTet(markedTet);

  std::cout << "adding hex ..." << std::endl;

  gmds::VTKWriter vtkWriterHT(&ioService);
  vtkWriterHT.setCellOptions(gmds::N|gmds::R);
  vtkWriterHT.setDataOptions(gmds::N|gmds::R);
  vtkWriterHT.write("ModeleCAD1_HEX_TET0_BISBF.vtk");
}
