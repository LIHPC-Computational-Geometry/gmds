/*----------------------------------------------------------------------------*/
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
using namespace gmds;
using namespace hybrid;
using namespace operators;
using namespace simplicesNode;
using namespace simplicesTriangle;
using namespace simplicesCell;
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  std::cout << "==== Curve indexed correction ====" << std::endl;
  std::string fIn, fOut;
  if(argc != 2)
  {
      throw gmds::GMDSException("NO INPUT FILE : <mesh_file>");
  }
  fIn = std::string(argv[1]);
  if (fIn.find('.vtk') == std::string::npos) {
    throw gmds::GMDSException("<mesh_file> NOT A .vtk FILE");
  }

  std::string extansion(".vtk");
  std::size_t position = fIn.find(extansion);
  fOut = fIn.substr(0,position) +  "_CORRECTED.vtk";


  std::cout << "INPUT FILE: " << fIn << std::endl;
  std::cout << "OUTPUT FILE: " << fOut << std::endl;

  //==================================================================
  // MESH FILE READING
  //==================================================================
  SimplexMesh simplexMesh = SimplexMesh();
  gmds::ISimplexMeshIOService ioService(&simplexMesh);
  gmds::VTKReader vtkReader(&ioService);
  vtkReader.setCellOptions(gmds::R|gmds::N);
  vtkReader.setDataOptions(gmds::N);
  vtkReader.read(fIn);
  simplexMesh.buildAdjInfoGlobal();
  simplexMesh.initializeEdgeStructure();
  simplexMesh.buildSimplexHull();

  Variable<int>* BND_CURVE_COLOR = nullptr;
  Variable<int>* BND_SURFACE_COLOR = nullptr;
  Variable<int>* BND_VERTEX_COLOR = nullptr;

  try{
    BND_CURVE_COLOR  = simplexMesh.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
    BND_SURFACE_COLOR  = simplexMesh.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
    BND_VERTEX_COLOR  = simplexMesh.getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
  }catch(gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  std::set<TInt> wronglabeledCurve{};
  const gmds::BitVector& nodesBitvec = simplexMesh.getBitVectorNodes();

  for(unsigned int node_id = 0 ; node_id < nodesBitvec.capacity() ; node_id++)
  {
    if(nodesBitvec[node_id] != 0)
    {
      if((*BND_CURVE_COLOR)[node_id] != 0)
      {
        const std::vector<TInt>&& neighboorNodes = SimplicesNode(&simplexMesh, node_id).neighborNodes();
        std::vector<TInt> nodesInSameCurve{};
        for(auto const neighboorNode : neighboorNodes)
        {
          if((*BND_CURVE_COLOR)[neighboorNode] == (*BND_CURVE_COLOR)[node_id] || (*BND_VERTEX_COLOR)[neighboorNode] != 0)
          {
            nodesInSameCurve.push_back(neighboorNode);
          }
        }

        if(nodesInSameCurve.size() < 2)
        {
          wronglabeledCurve.insert((*BND_CURVE_COLOR)[node_id]);
        }
        else if(nodesInSameCurve.size() > 2)
        {
          //peut se produire meme pour des maillages valide exemple node 286 maillage S23.vtk
          /*std::cout << "nodesInSameCurve.size() -> " << nodesInSameCurve.size() << std::endl;
          std::cout << "problem with the node -> " << node_id << std::endl;
          SimplexMesh nodeMesh = SimplexMesh();
          nodeMesh.addNode(SimplicesNode(&simplexMesh, node_id).getCoords()[0], SimplicesNode(&simplexMesh, node_id).getCoords()[1], SimplicesNode(&simplexMesh, node_id).getCoords()[2]);
          nodeMesh.addTetraedre(0, 0, 0, 0);
          gmds::ISimplexMeshIOService ioServiceNode(&nodeMesh);
          gmds::VTKWriter vtkWriterNode(&ioServiceNode);
          vtkWriterNode.setCellOptions(gmds::N|gmds::R);
          vtkWriterNode.setDataOptions(gmds::N|gmds::R);
          vtkWriterNode.write("NODE_" + std::to_string(node_id) + ".vtk");
          throw gmds::GMDSException("nodesInSameCurve.size() > 2");*/
        }
      }
    }
  }

  const std::multimap<TInt, std::pair<TInt,TInt>> & edgesStructure = simplexMesh.getConstEdgeStructure();
  std::vector<TInt> newsLabel{};
  for(auto const wrongLabel : wronglabeledCurve)
  {
    bool flag = false;
    auto it = edgesStructure.equal_range(wrongLabel);
    for(auto itr = it.first ; itr != it.second ; itr++)
    {
      TInt node0 = itr->second.first;
      TInt node1 = itr->second.second;
      std::vector<TInt> edge{node0, node1};
      std::vector<TSimplexID> shell = SimplicesNode(&simplexMesh, node0).shell(SimplicesNode(&simplexMesh, node1));
      for(auto const simplex : shell)
      {
        if(simplex < 0) // triangle
        {
          std::vector<TInt> v = SimplicesTriangle(&simplexMesh, -simplex).getOtherNodeInSimplex(edge);
          if(v.size() == 1)
          {
            if((*BND_SURFACE_COLOR)[v.front()] != 0)
            {
              newsLabel.push_back((*BND_SURFACE_COLOR)[v.front()]);
              flag = true;
              break;
            }
          }
          else
          {
            throw gmds::GMDSException("v.size() != 1, in labeled correction algorithm");
          }
        }
      }
      if(flag)
      {
        break;
      }
    }
  }

  if(newsLabel.size() != wronglabeledCurve.size()){
    std::cout << "newsLabel.size() -> " << newsLabel.size() << std::endl;
    std::cout << "wronglabeledCurve.size() -> " << wronglabeledCurve.size() << std::endl;
    throw gmds::GMDSException("newsLabel.size() != wronglabeledCurve.size()");
  }

  unsigned int cpt = 0;
  for(auto const wrongLabel : wronglabeledCurve)
  {
    TInt newLabel = newsLabel[cpt];
    auto it = edgesStructure.equal_range(wrongLabel);
    for(auto itr = it.first ; itr != it.second ; itr++)
    {
      TInt node0 = itr->second.first;
      TInt node1 = itr->second.second;
      BND_CURVE_COLOR->set(node0, 0);
      BND_CURVE_COLOR->set(node1, 0);
      BND_SURFACE_COLOR->set(node0, newLabel);
      BND_SURFACE_COLOR->set(node1, newLabel);
    }
    cpt++;
  }

  gmds::VTKWriter vtkWriter(&ioService);
  vtkWriter.setCellOptions(gmds::N|gmds::R);
  vtkWriter.setDataOptions(gmds::N|gmds::R);
  vtkWriter.write(fOut);
}

/*----------------------------------------------------------------------------*/
