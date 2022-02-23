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

  std::string extansion(".vtk");
  std::size_t position = fIn.find(extansion);
  fDI = fIn.substr(0,position) +  "_DELAUNAY_INSERTION.vtk";
  fER = fIn.substr(0,position) +  "_EDGES_REMOVE.vtk";
  fHEX = fIn.substr(0,position) + "_HEX_DOMINANT.vtk";
  fFF = fIn.substr(0,position) + "_FORCE_FACE.vtk";
  fEI = fIn.substr(0,position) + "_EDGE_INSERTION.vtk";

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
  return 0;
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

  const gmds::BitVector& nodesToAddIds = simplexNodes.getBitVectorNodes();
  const gmds::BitVector& nodePresentInMesh = simplexMesh.getBitVectorNodes();

  CriterionRAIS criterionRAIS(new VolumeCriterion());
  gmds::BitVector nodesAdded(simplexMesh.nodesCapacity());
  std::vector<TInt> nodes(nodesToAddIds.capacity(), -1);
  TInt border = std::numeric_limits<TInt>::min();

  //==================================================================
  // DElAUNAY INSERTION OF THE POINT DATA
  //==================================================================
  std::cout << "DELAUNAY POINT INSERTION START" << std::endl;
  std::clock_t start;
  double duration;
  start = std::clock();
  for(unsigned int idx = 0 ; idx < nodesToAddIds.capacity() ; idx++)
  {
    if(nodesToAddIds[idx] != 0)
    {
      const gmds::BitVector & nodesIds = simplexMesh.getBitVectorNodes();
      math::Point point = SimplicesNode(&simplexNodes, idx).getCoords();

      bool alreadyAdd = false;
      std::vector<TSimplexID> tetraContenaingPt{};
      TInt node = simplexMesh.addNodeAndcheck(point, tetraContenaingPt, alreadyAdd);
      if(!alreadyAdd)
      {
        if((*BND_CURVE_COLOR_NODES)[idx] != 0) {BND_CURVE_COLOR->set(node, (*BND_CURVE_COLOR_NODES)[idx]);}
        else if((*BND_SURFACE_COLOR_NODES)[idx] != 0) {BND_SURFACE_COLOR->set(node, (*BND_SURFACE_COLOR_NODES)[idx]);}

        simplexMesh.getVariable<Eigen::Matrix3d, SimplicesNode>("metric")->value(node) = m;
        bool status = false;
        std::vector<TSimplexID> deletedSimplex{};
        const std::multimap<TInt, TInt> facesAlreadyBuilt{};
        DelaunayPointInsertion DI(&simplexMesh, SimplicesNode(&simplexMesh, node), criterionRAIS, tetraContenaingPt, status, nodesAdded, deletedSimplex, facesAlreadyBuilt);
        if(status)
        {
          nodes[idx] = node;
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
        nodes[idx] = node;
        if(nodesAdded.capacity() != nodesIds.capacity())
        {
          nodesAdded.resize(nodesIds.capacity());
        }
        nodesAdded.assign(node);
      }
    }
  }

  duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
  std::cout << "DELAUNAY INSERTION DONE IN " << duration << std::endl;
  std::cout << "  INSERTED NODE -> "  << (double)nodesAdded.size() / (double)nodesToAddIds.size() * 100.0 << "% " << std::endl;
  std::cout << std::endl;
  gmds::VTKWriter vtkWriterDI(&ioService);
  vtkWriterDI.setCellOptions(gmds::N|gmds::R|gmds::F);
  vtkWriterDI.setDataOptions(gmds::N|gmds::R|gmds::F);
  vtkWriterDI.write(fDI);

  start = std::clock();
  std::vector<TInt> deletedNodes{};
  unsigned int tmp = 0;
  for(;;)
  {
    std::cout << "  EDGE REMOVING....." << std::endl;
    unsigned int edgesRemoved = simplexMesh.edgesRemove(nodesAdded, deletedNodes);
    if(edgesRemoved == tmp || edgesRemoved == 0)
    {
      break;
    }
    tmp = edgesRemoved;
  }
  std::cout << "deletedNodes.size() --> " << deletedNodes.size() << std::endl;
  std::cout << "EDGE REMOVING REINSERTION " << std::endl;
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
      const std::multimap<TInt, TInt> facesAlreadyBuilt{};
      bool status = false;
      PointInsertion(&simplexMesh, SimplicesNode(&simplexMesh, deletedNode), criterionRAIS, status, tetraContenaingPt, nodesAdded, deletedNodes, facesAlreadyBuilt);
      if(status)
      {
        nodeReinsertedSize++;
      }
    }
  }

  duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
  std::cout << "EDGES REMOVED DONE IN " << duration << std::endl;
  std::cout << std::endl;
  gmds::VTKWriter vtkWriterER(&ioService);
  vtkWriterER.setCellOptions(gmds::N|gmds::R|gmds::F);
  vtkWriterER.setDataOptions(gmds::N|gmds::R|gmds::F);
  vtkWriterER.write(fER);


  ///////////////////////////EDGE BUILDER START HERE///////////////////////////
  std::multimap<TInt, TInt> edges{};
  const gmds::BitVector triIdx = simplexNodes.getBitVectorTri();
  unsigned int cptEdgeNotBuilt = 0;
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
        cptEdgeNotBuilt++;
        //std::cout << "e[0,1] --> " << e[0] << " | " << e[1] << std::endl;
      }
    }
  }
  std::cout << "EDGE NOT BUILT -> " << (double)cptEdgeNotBuilt / (double)triIdx.size() * 100.0 << " %" << std::endl;
  std::cout << std::endl;
  std::vector<std::vector<TSimplexID>> hexesNodes{};
  std::vector<std::vector<TInt>> nodesHex = simplexNodes.getHexadronData();
  gmds::BitVector markedTet;

  double hexBuiltCpt = 0;
  unsigned int faceBuildCpt = 0;
  unsigned int iter = 0;
  start = std::clock();
  std::cout << "HEX GENERATION START " << std::endl;
  std::multimap<TInt, TInt> edgeAlreadyBuilt{};
  for(;;)
  {
    hexBuiltCpt = 0;
    hexesNodes.clear();
    unsigned int tmp = 0;
    unsigned int ITER_MAX = 4;
    for(unsigned int iter = 0 ; iter < ITER_MAX ; iter++)
    {
      unsigned int edgeBuild = simplexMesh.buildEdges(edges, nodesAdded);
      if(edgeBuild == tmp || edgeBuild == 0)
      {
        break;
      }
      tmp = edgeBuild;
    }
    std::cout << "BUILD EDGE DONE " << std::endl;
    gmds::VTKWriter vtkWriterEI(&ioService);
    vtkWriterEI.setCellOptions(gmds::N|gmds::R);
    vtkWriterEI.setDataOptions(gmds::N|gmds::R);
    vtkWriterEI.write(fEI);
    /////////////////////HEXA'S FACES BUILDER START HERE /////////////////////////


    unsigned int faceBuiltTmp = 0;
    std::cout << "FACE BUILDING START " << std::endl;
    for(auto const h : nodesHex)
    {
      TInt n0 = nodes[h[0]]; TInt n1 = nodes[h[1]]; TInt n2 = nodes[h[2]]; TInt n3 = nodes[h[3]];
      TInt n4 = nodes[h[4]]; TInt n5 = nodes[h[5]]; TInt n6 = nodes[h[6]]; TInt n7 = nodes[h[7]];
      std::vector<TInt> hexeNodes{n0, n1, n2, n3, n4, n5, n6, n7};
      if(n0 != -1 && n1 != -1 && n2 != -1 && n3 != -1 && n4 != -1 && n5 != -1 && n6 != -1 && n7 != -1)
      {
        if(simplexMesh.buildFace(hexeNodes, nodesAdded, edgeAlreadyBuilt))
        {

          faceBuiltTmp++;
        }
      }
    }
    gmds::VTKWriter vtkWriterFF(&ioService);
    vtkWriterFF.setCellOptions(gmds::N|gmds::R);
    vtkWriterFF.setDataOptions(gmds::N|gmds::R);
    vtkWriterFF.write(fFF);
    std::cout << "BUILD FACE DONE " << std::endl;
    markedTet = simplexMesh.getBitVectorTet();

    for(auto const h : nodesHex)
    {
      TInt n0 = nodes[h[0]]; TInt n1 = nodes[h[1]]; TInt n2 = nodes[h[2]]; TInt n3 = nodes[h[3]];
      TInt n4 = nodes[h[4]]; TInt n5 = nodes[h[5]]; TInt n6 = nodes[h[6]]; TInt n7 = nodes[h[7]];
      std::vector<TInt> hexeNodes{n0, n1, n2, n3, n4, n5, n6, n7};
      if(n0 != -1 && n1 != -1 && n2 != -1 && n3 != -1 && n4 != -1 && n5 != -1 && n6 != -1 && n7 != -1)
      {
        std::vector<TSimplexID> simplices = simplexMesh.hex2tet(hexeNodes);
        if(simplices.size() != 0)
        {
          for(auto const simplex : simplices){ markedTet.unselect(simplex);}
          hexBuiltCpt = hexBuiltCpt + 1.0;
          hexesNodes.push_back(hexeNodes);
        }
      }
    }
    std::cout << "HEX 2 TET DONE " << std::endl;

    if(faceBuiltTmp == faceBuildCpt || faceBuiltTmp == 0)
    {
      break;
    }
    faceBuildCpt = faceBuiltTmp;
  }

  //RESULT DATA
  duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
  std::cout << "HEX GENERATION DONE " << duration << std::endl;
  std::cout << "  hex build % -> " << hexBuiltCpt / (double)nodesHex.size() << std::endl;
  std::cout << "  hexBuiltCpt -> " << hexBuiltCpt << std::endl;
  std::cout << "  hexesNodes.size() --> " << nodesHex.size() << std::endl;


  simplexMesh.setHexadronData(hexesNodes);
  simplexMesh.setMarkedTet(markedTet);

  gmds::VTKWriter vtkWriterHT(&ioService);
  vtkWriterHT.setCellOptions(gmds::N|gmds::R);
  vtkWriterHT.setDataOptions(gmds::N|gmds::R);
  vtkWriterHT.write(fHEX);
}

/*----------------------------------------------------------------------------*/
