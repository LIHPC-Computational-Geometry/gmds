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
          if((*BND_CURVE_COLOR)[neighboorNode] == (*BND_CURVE_COLOR)[node_id] ||
              (*BND_VERTEX_COLOR)[neighboorNode] != 0)
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
          std::cout << "problem with the node -> " << node_id << std::endl;
          throw gmds::GMDSException("nodesInSameCurve.size() > 2");
        }
      }
    }
  }

  const std::multimap<TInt, std::pair<TInt,TInt>> & edgesStructure = simplexMesh.getConstEdgeStructure();
  for(auto const wrongLabel : wronglabeledCurve)
  {
    auto it = edgesStructure.equal_range(wrongLabel);
    for(auto itr = it.first ; itr != it.second ; itr++)
    {
      TInt node0 = itr->second.first;
      TInt node1 = itr->second.second;
      BND_CURVE_COLOR->set(node0, 0);
      BND_CURVE_COLOR->set(node1, 0);
    }
  }

  gmds::VTKWriter vtkWriter(&ioService);
  vtkWriter.setCellOptions(gmds::N|gmds::R);
  vtkWriter.setDataOptions(gmds::N|gmds::R);
  vtkWriter.write(fOut);
}

/*----------------------------------------------------------------------------*/
