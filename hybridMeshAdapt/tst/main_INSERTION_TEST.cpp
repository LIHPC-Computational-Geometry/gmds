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
using namespace gmds;
using namespace hybrid;
using namespace operators;
using namespace simplicesNode;
using namespace simplicesTriangle;
using namespace simplicesCell;
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  std::string fIn, pIn, fDI, fER, fHEX, fFF, fEI;
  if(argc != 3)
  {
      throw gmds::GMDSException("NO INPUT FILE : <mesh_file> <point_file>");
  }
  fIn = std::string(argv[1]);
  pIn = std::string(argv[2]);
  if (fIn.find('.vtk') == std::string::npos) {
    throw gmds::GMDSException("<mesh_file> NOT A .vtk FILE");
  }
  if (pIn.find('.vtk') == std::string::npos) {
    throw gmds::GMDSException("<point_file> NOT A .vtk FILE");
  }
  std::cout << "INPUT FILE: " << fIn << std::endl;
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

  Octree oc(&simplexMesh, 50);
  simplexMesh.setOctree(&oc);
  Variable<int>* BND_VERTEX_COLOR  = simplexMesh.getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
  Variable<int>* BND_CURVE_COLOR   = simplexMesh.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR = simplexMesh.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");

  //adding a metric to the mesh for the delaunay expansion ctriterion
  Variable<Eigen::Matrix3d>* var = simplexMesh.newVariable<Eigen::Matrix3d, SimplicesNode>("metric");
  Eigen::Matrix3d m =  Eigen::MatrixXd::Identity(3, 3);
  var->setValuesTo(m);
  //==================================================================
  // POINT FILE READING
  //==================================================================
  SimplexMesh simplexNodes = SimplexMesh();
  gmds::ISimplexMeshIOService ioServiceNodes(&simplexNodes);
  gmds::VTKReader vtkReaderNodes(&ioServiceNodes);
  vtkReaderNodes.setCellOptions(gmds::R|gmds::N|gmds::F);
  vtkReaderNodes.setDataOptions(gmds::N);
  vtkReaderNodes.read(pIn);
  Variable<int>* BND_CURVE_COLOR_NODES   = simplexNodes.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  Variable<int>* BND_SURFACE_COLOR_NODES = simplexNodes.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
  Variable<int>* BND_VERTEX_COLOR_NODES = simplexNodes.getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");

  const gmds::BitVector& nodesToAddIds = simplexNodes.getBitVectorNodes();
  const gmds::BitVector& nodePresentInMesh = simplexMesh.getBitVectorNodes();

  CriterionRAIS criterionRAIS(new VolumeCriterion());
  gmds::BitVector nodesAdded(simplexMesh.nodesCapacity());
  std::vector<TInt> nodes(nodesToAddIds.capacity(), -1);
  TInt border = std::numeric_limits<TInt>::min();

  //==================================================================
  // VOLUME POINT INSERTION
  //==================================================================
  std::cout << "VOLUME POINT INSERTION START" << std::endl;
  std::clock_t start;
  double duration;
  start = std::clock();
  unsigned int nodeCpt = 0;
  for(unsigned int idx = 0 ; idx < nodesToAddIds.capacity() ; idx++)
  {
    std::cout << "idx -> " << idx << std::endl;
    if(nodesToAddIds[idx] != 0)
    {
      const gmds::BitVector & nodesIds = simplexMesh.getBitVectorNodes();
      math::Point point = SimplicesNode(&simplexNodes, idx).getCoords();

      bool alreadyAdd = false;
      std::vector<TSimplexID> tetraContenaingPt{};
      TInt node = simplexMesh.addNodeAndcheck(point, tetraContenaingPt, alreadyAdd);
      std::cout << "node -> " << node << std::endl;
      if(!alreadyAdd)
      {
        if((*BND_CURVE_COLOR_NODES)[idx] != 0) {BND_CURVE_COLOR->set(node, (*BND_CURVE_COLOR_NODES)[idx]);}
        else if((*BND_SURFACE_COLOR_NODES)[idx] != 0) {BND_SURFACE_COLOR->set(node, (*BND_SURFACE_COLOR_NODES)[idx]);}
        else if((*BND_VERTEX_COLOR_NODES)[idx] != 0) {BND_VERTEX_COLOR->set(node, (*BND_VERTEX_COLOR_NODES)[idx]);}

        if((*BND_CURVE_COLOR)[node] == 0 && (*BND_SURFACE_COLOR)[node] == 0 && (*BND_SURFACE_COLOR)[node] == 0)
        {
          bool status = false;
          std::vector<TSimplexID> deletedSimplex{};
          const std::multimap<TInt, TInt> facesAlreadyBuilt{};
          std::cout << "DelaunayPointInsertion " << std::endl;
          DelaunayPointInsertion DI(&simplexMesh, SimplicesNode(&simplexMesh, node), criterionRAIS, tetraContenaingPt, status, nodesAdded, deletedSimplex, facesAlreadyBuilt);
          std::cout << "status -> " << status << std::endl;
          std::cout << std::endl;
          if(status)
          {
            nodeCpt++;
          }
        }
      }
    }
  }

  duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
  std::cout << "DELAUNAY VOLUME INSERTION DONE IN " << duration << std::endl;
  std::cout << "  INSERTED NODE -> "  << (double)nodeCpt / (double)nodesToAddIds.size() * 100.0 << "% " << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
}

/*----------------------------------------------------------------------------*/
