/*----------------------------------------------------------------------------*/
#include <gmds/io/VTKReader.h>
#include <gmds/io/VTKWriter.h>
/*----------------------------------------------------------------------------*/
#include <gmds/hybridMeshAdapt/ICriterion.h>
#include <gmds/hybridMeshAdapt/SimplexMesh.h>
#include <gmds/hybridMeshAdapt/ISimplexMeshIOService.h>
#include <gmds/hybridMeshAdapt/DelaunayPointInsertion.h>
#include <gmds/hybridMeshAdapt/Octree.h>
#include <gmds/hybridMeshAdapt/MetricFFPointgeneration.h>
#include <gmds/math/Hexahedron.h>
#include <gmds/ig/Mesh.h>
#include <gmds/ig/MeshDoctor.h>
#include <gmds/io/IGMeshIOService.h>
/*----------------------------------------------------------------------------*/
#include <iostream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace math;
using namespace hybrid;
using namespace operators;
using namespace simplicesNode;
using namespace simplicesTriangle;
using namespace simplicesCell;
/*----------------------------------------------------------------------------*/
std::set<TInt> findNodeHullCpt(const std::vector<TSimplexID>& shell, SimplexMesh* simplexMesh, bool& status, std::set<TSimplexID>& seen_tet);
double computeTheScaledJacobianOfTheHull(const std::vector<TSimplexID>& shell, SimplexMesh* simplexMesh, std::vector<TSimplexID>& hexNode, std::multimap<TInt, TInt>& mmap_edges, double sJ);
std::multimap<TInt, TInt> sortEdgeNearSurface(SimplexMesh* simplexMesh);
/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    std::cout << "==== METRIC BASED ADAPTATION  ====" << std::endl;
    std::string fIn, fOut;
    if(argc != 2)
    {
        throw gmds::GMDSException("MISSING PARAMETERS : <mesh_in>");
    }
    fIn = std::string(argv[1]);
    if (fIn.find('.vtk') == std::string::npos) {
      throw gmds::GMDSException("<mesh_in> NOT A .vtk FILE");
    }
    std::string extansion(".vtk");
    std::size_t position = fIn.find(extansion);
    fOut = fIn.substr(0,position) + "_TET2HEX.vtk";
    std::cout << "INPUT FILE: " << fIn << std::endl;
    std::cout << "OUTPUT FILE: " << fOut << std::endl;

    ////////////////////////////////////////////////////////////////////////
    /*std::cout << "Reading " << std::endl;
    std::string extansion(".vtk");
    std::size_t position = fIn.find(extansion);
    fOut = fIn.substr(0,position) + "_TET2HEX.vtk";
    std::cout << "INPUT FILE: " << fIn << std::endl;
    std::cout << "OUTPUT FILE: " << fOut << std::endl;
    Mesh cube(MeshModel(DIM3 | R | F | E | N |
      R2N | F2N | E2N | R2F | F2R |
      F2E | E2F | R2E | N2R | N2F | N2E));

      const Node nodeHEX0 =  cube.newNode(math::Point({0.0,0.0,0.0}));
      const Node nodeHEX1 =  cube.newNode(math::Point({0.0,0.0,1.0}));
      const Node nodeHEX2 =  cube.newNode(math::Point({1.0,0.0,1.0}));
      const Node nodeHEX3 =  cube.newNode(math::Point({1.0,0.0,0.0}));


      const Node nodeHEX4 =  cube.newNode(math::Point({-0.2,1.0,-0.2}));
      const Node nodeHEX5 =  cube.newNode(math::Point({0.0,1.0,1.2}));
      const Node nodeHEX6 =  cube.newNode(math::Point({1.0,1.2,1.2}));
      const Node nodeHEX7 =  cube.newNode(math::Point({1.2,1.1,0.0}));

      cube.newHex(nodeHEX0, nodeHEX1, nodeHEX2, nodeHEX3,
              nodeHEX4, nodeHEX5, nodeHEX6, nodeHEX7);
      gmds::IGMeshIOService ioServiceCUBE(&cube);
      gmds::VTKWriter vtkWriterCUBE(&ioServiceCUBE);
      vtkWriterCUBE.setCellOptions(gmds::N|gmds::R);
      std::string extansionCUBE("HEX_0_");
      std::size_t positionCUBE = fIn.find(extansionCUBE);
      std::string fOutCube = fIn.substr(0,positionCUBE) + "CUBE.vtk";

      vtkWriterCUBE.write(fOutCube);
      return 0;
      std::cout << "HEX CREATED " << std::endl;*/
      ////////////////////////////////////////////////////////////////////////


    SimplexMesh simplexMesh = SimplexMesh();
    gmds::ISimplexMeshIOService ioService(&simplexMesh);
    gmds::VTKReader vtkReader(&ioService);
    vtkReader.setCellOptions(gmds::R|gmds::N);
    vtkReader.setDataOptions(gmds::N);
    std::cout << "Reading the input file" << std::endl;
    vtkReader.read(fIn);
    simplexMesh.buildAdjInfoGlobal();
    std::cout << "Adjacent info built" << std::endl;
    simplexMesh.initializeEdgeStructure();
    std::cout << "Edge structure built" << std::endl;

    SimplexMesh testMesh = SimplexMesh();
    const gmds::BitVector& tet_bit = simplexMesh.getBitVectorTet();
    unsigned int sizeTetra_init = tet_bit.size();
    std::vector<TInt> f0_{};
    std::vector<TInt> f1_{};
    std::vector<TInt> f2_{};
    std::vector<TInt> f3_{};
    std::unordered_set<TSimplexID> seen{};
    gmds::BitVector markedTet = simplexMesh.getBitVectorTet();
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // init edge i mesh
    const gmds::BitVector & nodeBitVector = simplexMesh.getBitVectorNodes();
    std::multimap<TInt, TInt> edges{};
    std::vector<TSimplexID> tetsShell{};
    std::set<std::pair<TInt,TInt>> edges_seen{};
    for(TInt n = 0 ; n < nodeBitVector.capacity() ; n++)
    {
      if(nodeBitVector[n] != 0)
      {
        const std::vector<TSimplexID> ballOf_n = SimplicesNode(&simplexMesh, n).ballOf();
        std::set<TInt> s{};
        for(auto const tet : ballOf_n)
        {
          const std::vector<TSimplexID> nodes = SimplicesCell(&simplexMesh, tet).getNodes();
          for(auto const nN : nodes)
          {
            s.insert(nN);
          }
        }

        for(auto const nN : s)
        {
          if(nN != n)
          {
            if(edges_seen.find(std::pair<TInt, TInt>{std::min(nN,n), std::max(nN,n)}) == edges_seen.end())
            {
              edges_seen.insert(std::pair<TInt, TInt>{std::min(nN,n), std::max(nN,n)});
              edges.insert(std::pair<TInt, TInt>{std::min(nN,n), std::max(nN,n)});
            }
          }
        }
      }
    }

    std::vector<std::vector<TInt>> sortedEdges_8{};
    std::vector<std::vector<TInt>> sortedEdges_7{};
    std::vector<std::vector<TInt>> sortedEdges_6{};
    std::vector<std::vector<TInt>> sortedEdges_5{};
    std::vector<std::vector<TInt>> sortedEdges_4{};
    std::vector<std::vector<TInt>> sortedEdges_3{};

    for(auto const edge : edges)
    {
      TInt nodeA = edge.first;
      TInt nodeB = edge.second;

      std::vector<TSimplexID> shell = SimplicesNode(&simplexMesh, nodeA).shell(SimplicesNode(&simplexMesh, nodeB));
      // parameters optional here :
      bool status = false;
      std::set<TSimplexID> seen_tet{};
      std::set<TInt> s = findNodeHullCpt(shell, &simplexMesh, status, seen_tet);
      if(s.size() == 8)
        sortedEdges_8.push_back({nodeA, nodeB});
      if(s.size() == 7)
        sortedEdges_7.push_back({nodeA, nodeB});
      if(s.size() == 6)
        sortedEdges_6.push_back({nodeA, nodeB});
    }

    std::vector<std::vector<TInt>> sortedEdges{};
    std::multimap<TInt, TInt> mmap_edges{};
    for(auto const & edge : sortedEdges_8){
      sortedEdges.push_back(edge);
      mmap_edges.insert(std::pair<TInt, TInt>{edge.front(), edge.back()});
    }

    for(auto const & edge : sortedEdges_7){
      sortedEdges.push_back(edge);
      mmap_edges.insert(std::pair<TInt, TInt>{edge.front(), edge.back()});
    }
    for(auto const & edge : sortedEdges_6){
      sortedEdges.push_back(edge);
      mmap_edges.insert(std::pair<TInt, TInt>{edge.front(), edge.back()});
    }
    for(auto const & edge : sortedEdges_5){
      sortedEdges.push_back(edge);
      mmap_edges.insert(std::pair<TInt, TInt>{edge.front(), edge.back()});
    }
    for(auto const & edge : sortedEdges_4){
      sortedEdges.push_back(edge);
      mmap_edges.insert(std::pair<TInt, TInt>{edge.front(), edge.back()});
    }
    for(auto const & edge : sortedEdges_3){
      sortedEdges.push_back(edge);
      mmap_edges.insert(std::pair<TInt, TInt>{edge.front(), edge.back()});
    }

    std::set<TSimplexID> seen_tet{};
    std::vector<std::vector<TInt>> hexs{};
    std::vector<double> sJ_values{0.90, 0.80, 0.70, 0.60, 0.50, 0.4, 0.3, 0.2};
    std::set<std::pair<TInt, TInt>> seen_Edges{};

    //std::multimap<TInt, TInt> sortedEdgeSurface = sortEdgeNearSurface(&simplexMesh);
    for(auto const sJ_ : sJ_values)
    {
      std::cout << "sJ -> " << sJ_ << std::endl;

      //for(auto const edge : sortedEdgeSurface)
      for(auto const edge : sortedEdges)
      {
        TInt nodeA = edge.front();
        TInt nodeB = edge.back();

        /*TInt nodeA = edge.first;
        TInt nodeB = edge.second;*/

        if(seen_Edges.find(std::pair<TInt, TInt>{nodeA, nodeB}) != seen_Edges.end())
          continue;

        //if(!(nodeA == 4842 && nodeB == 7444))
          //continue;
        std::vector<TSimplexID> shell = SimplicesNode(&simplexMesh, nodeA).shell(SimplicesNode(&simplexMesh, nodeB));
        shell.erase(std::remove_if(shell.begin(), shell.end(), [&](TSimplexID t){
          return (std::find(seen_tet.begin(), seen_tet.end(), t) != seen_tet.end());
        }), shell.end());
        bool status = false;
        std::set<TInt> s = findNodeHullCpt(shell, &simplexMesh, status, seen_tet);
        if(status)
          continue;

        if(s.size() == 8)
        {
          std::vector<TSimplexID> hexNode{};
          double sJ = computeTheScaledJacobianOfTheHull(shell, &simplexMesh, hexNode, mmap_edges, sJ_);
          if(hexNode.size() == 8)
          {
            hexs.push_back(hexNode);
            for(auto const t : shell)
              seen_tet.insert(t);
            seen_Edges.insert(std::pair<TInt, TInt>{nodeA, nodeB});
          }
        }
        else if(s.size() == 7)
        {
          double sJGlobal = std::numeric_limits<double>::lowest();
          TSimplexID bestAdjSimplex = std::numeric_limits<TSimplexID>::min();
          std::vector<TSimplexID> hexNode{};
          //add exterior tet to add a node to the hull
          //to build the hexahedron
          for(auto const tet : shell)
          {
            for(unsigned int f = 0 ; f < 4 ; f++)
            {
              TSimplexID adjTet = SimplicesCell(&simplexMesh, tet).oppositeTetraIdx(f);
              if(adjTet >= 0)
              {
                if(seen_tet.find(adjTet) != seen_tet.end())
                  continue;

                if(std::find(shell.begin(), shell.end(), adjTet) == shell.end())
                {
                  std::vector<TSimplexID> new_hull(shell.begin(), shell.end());
                  new_hull.push_back(adjTet);
                  bool status = false;
                  std::set<TInt> new_s = findNodeHullCpt(new_hull, &simplexMesh, status, seen_tet);
                  if(status)
                    continue;

                  if(new_s.size() == 8) // possible hex
                  {
                    std::vector<TInt> nodes = SimplicesCell(&simplexMesh, adjTet).getNodes();
                    std::vector<TSimplexID> h{};
                    //std::vector<TSimplexID> simplicesForFacesABCD{};

                    double sJ = computeTheScaledJacobianOfTheHull(new_hull, &simplexMesh, h, mmap_edges, sJ_);
                    if(sJGlobal < sJ && h.size() == 8)
                    {
                      sJGlobal = sJ;
                      bestAdjSimplex = adjTet;
                      hexNode = h;
                    }
                  }
                }
              }
            }
          }
          if(hexNode.size() == 8){
            hexs.push_back(hexNode);
            for(auto const t : shell)
              seen_tet.insert(t);
            seen_tet.insert(bestAdjSimplex);
            seen_Edges.insert(std::pair<TInt, TInt>{nodeA, nodeB});
          }
        }
        else if(s.size() == 6)
        {
          double sJGlobal = std::numeric_limits<double>::lowest();
          std::vector<TSimplexID> bestAdjSimplices{};
          std::vector<TSimplexID> hexNode{};
          //add exterior tet to add a node to the hull
          //to build the hexahedron
          for(auto const tet_first : shell)
          {
            for(unsigned int f0 = 0 ; f0 < 4 ; f0++)
            {
              TSimplexID adjTet_first = SimplicesCell(&simplexMesh, tet_first).oppositeTetraIdx(f0);
              if(adjTet_first >= 0)
              {
                if(seen_tet.find(adjTet_first) != seen_tet.end() ||
                    std::find(shell.begin(), shell.end(), adjTet_first) != shell.end())
                  continue;

                  for(auto const tet_second : shell)
                  {
                    if(tet_second != tet_first && tet_second >= 0)
                    {
                      for(unsigned int f1 = 0 ; f1 < 4 ; f1++)
                      {
                        TSimplexID adjTet_second = SimplicesCell(&simplexMesh, tet_second).oppositeTetraIdx(f1);
                        if(adjTet_second >= 0)
                        {
                          if(seen_tet.find(adjTet_second) != seen_tet.end() ||
                              std::find(shell.begin(), shell.end(), adjTet_second) != shell.end())
                              continue;

                          std::vector<TSimplexID> new_hull(shell.begin(), shell.end());
                          new_hull.push_back(adjTet_first);
                          new_hull.push_back(adjTet_second);
                          bool status = false;
                          std::set<TInt> new_s = findNodeHullCpt(new_hull, &simplexMesh, status, seen_tet);
                          if(status)
                            continue;

                          if(new_s.size() == 8) // possible hex
                          {
                            std::vector<TSimplexID> h{};
                            double sJ = computeTheScaledJacobianOfTheHull(new_hull, &simplexMesh, h, mmap_edges, sJ_);
                            if(sJGlobal < sJ && h.size() == 8)
                            {
                              sJGlobal = sJ;
                              bestAdjSimplices.clear();
                              bestAdjSimplices.push_back(adjTet_first);
                              bestAdjSimplices.push_back(adjTet_second);
                              hexNode = h;
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
          }
          if(hexNode.size() == 8){
            hexs.push_back(hexNode);
            for(auto const t : shell){
              seen_tet.insert(t);
            }
            seen_tet.insert(bestAdjSimplices.front());
            seen_tet.insert(bestAdjSimplices.back());
            seen_Edges.insert(std::pair<TInt, TInt>{nodeA, nodeB});
          }
        }
        else if(s.size() == 5)
        {
          double sJGlobal = std::numeric_limits<double>::lowest();
          std::vector<TSimplexID> bestAdjSimplices{};
          std::vector<TSimplexID> hexNode{};
          //add exterior tet to add a node to the hull
          //to build the hexahedron
          for(auto const tet_first : shell)
          {
            for(unsigned int f0 = 0 ; f0 < 4 ; f0++)
            {
              TSimplexID adjTet_first = SimplicesCell(&simplexMesh, tet_first).oppositeTetraIdx(f0);
              if(adjTet_first >= 0)
              {
                if(seen_tet.find(adjTet_first) != seen_tet.end() ||
                    std::find(shell.begin(), shell.end(), adjTet_first) != shell.end())
                  continue;

                  for(auto const tet_second : shell)
                  {
                    if(tet_second != tet_first && tet_second >= 0)
                    {
                      for(unsigned int f1 = 0 ; f1 < 4 ; f1++)
                      {
                        TSimplexID adjTet_second = SimplicesCell(&simplexMesh, tet_second).oppositeTetraIdx(f1);
                        bool flag = false;
                        for(unsigned int f = 0 ; f < 4 ; f++)
                        {
                          if(adjTet_second == SimplicesCell(&simplexMesh, adjTet_first).oppositeTetraIdx(f))
                          {
                            flag = true;
                            break;
                          }
                        }
                        if(flag)
                          continue;

                        if(adjTet_second >= 0)
                        {
                          if(seen_tet.find(adjTet_second) != seen_tet.end() ||
                              std::find(shell.begin(), shell.end(), adjTet_second) != shell.end())
                          continue;
                          for(auto const tet_third : shell)
                          {
                            if(tet_second != tet_third && tet_third >= 0)
                            {
                              for(unsigned int f2 = 0 ; f2 < 4 ; f2++)
                              {
                                TSimplexID adjTet_third = SimplicesCell(&simplexMesh, tet_third).oppositeTetraIdx(f2);

                                if(adjTet_third >= 0)
                                {
                                  if(seen_tet.find(adjTet_third) != seen_tet.end() ||
                                      std::find(shell.begin(), shell.end(), adjTet_third) != shell.end())
                                    continue;
                                  std::vector<TSimplexID> new_hull(shell.begin(), shell.end());
                                  new_hull.push_back(adjTet_first);
                                  new_hull.push_back(adjTet_second);
                                  new_hull.push_back(adjTet_third);
                                  bool status = false;
                                  std::set<TInt> new_s = findNodeHullCpt(new_hull, &simplexMesh, status, seen_tet);
                                  if(status)
                                    continue;

                                  if(new_s.size() == 8) // possible hex
                                  {
                                    std::vector<TSimplexID> h{};
                                    double sJ = computeTheScaledJacobianOfTheHull(new_hull, &simplexMesh, h, mmap_edges, sJ_);
                                    if(sJGlobal < sJ && h.size() == 8)
                                    {
                                      sJGlobal = sJ;
                                      bestAdjSimplices.clear();
                                      bestAdjSimplices.push_back(adjTet_first);
                                      bestAdjSimplices.push_back(adjTet_second);
                                      bestAdjSimplices.push_back(adjTet_third);
                                      hexNode = h;
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
            if(hexNode.size() == 8){
              hexs.push_back(hexNode);
              for(auto const t : shell)
                seen_tet.insert(t);

              seen_tet.insert(bestAdjSimplices.front());
              seen_tet.insert(bestAdjSimplices[1]);
              seen_tet.insert(bestAdjSimplices.back());
            }
          }
        }
      std::cout << "hexs.size() -> " << hexs.size() << std::endl;
    }


    std::cout << "hexs.size() -> " << hexs.size() << std::endl;
    std::cout << "tetra nbr initial : "  << sizeTetra_init << std::endl;
    std::cout << "marked tetrahedron  : "  << seen_tet.size() << std::endl;
    std::cout << "remaining tetrahedron  : "  << sizeTetra_init - seen_tet.size() << std::endl;
    std::cout << "remaining tetrahedron in % : "  << static_cast<double>(sizeTetra_init - seen_tet.size()) / static_cast<double>(sizeTetra_init) << std::endl;

    //look up for the sliver

    Mesh m(MeshModel(DIM3 | R | F | E | N |
      R2N | F2N | E2N | R2F | F2R |
      F2E | E2F | R2E | N2R | N2F | N2E));

      std::unordered_map<TInt, Node> map{};
      for(unsigned int n = 0 ; n < nodeBitVector.capacity() ; n++)
      {
        if(nodeBitVector[n] != 0)
        {
          const Node nodeHEX =  m.newNode(SimplicesNode(&simplexMesh, n).getCoords());
          map[n] = nodeHEX;
        }
      }
      for(auto const & hex : hexs){
        m.newHex(map[hex[0]], map[hex[1]], map[hex[2]], map[hex[3]], map[hex[4]], map[hex[5]], map[hex[6]], map[hex[7]]);
      }
      for(auto const & hex : simplexMesh.getHexadronData()){
        m.newHex(map[hex[0]], map[hex[1]], map[hex[2]], map[hex[3]], map[hex[4]], map[hex[5]], map[hex[6]], map[hex[7]]);
      }

      const gmds::BitVector & tets = simplexMesh.getBitVectorTet();
      for(unsigned int tet_id = 0 ; tet_id < tets.capacity() ; tet_id++){
        if(tets[tet_id] != 0)
        {
          std::vector<TSimplexID> adjCells = SimplicesCell(&simplexMesh, tet_id).adjacentTetra();
          if(adjCells[0] < 0 || adjCells[1] < 0 || adjCells[2] < 0 || adjCells[3] < 0)
          {
            //continue;
          }
          if(!SimplicesCell(&simplexMesh, tet_id).isSliver())
          {
            if(std::find(seen_tet.begin(), seen_tet.end(), tet_id) == seen_tet.end())
            {
              std::vector<TInt> nodes = SimplicesCell(&simplexMesh, tet_id).getNodes();
              m.newTet(map[nodes[0]], map[nodes[1]], map[nodes[2]], map[nodes[3]]);
            }
          }
        }
      }

      gmds::IGMeshIOService ioServiceMESH(&m);
      gmds::VTKWriter vtkWriterMESH(&ioServiceMESH);
      vtkWriterMESH.setCellOptions(gmds::N|gmds::R);
      vtkWriterMESH.setDataOptions(gmds::N|gmds::R);
      vtkWriterMESH.write(fOut);
      std::cout << "HEX CREATED " << std::endl;

      return 0;
}
/*----------------------------------------------------------------------------*/
std::set<TInt> findNodeHullCpt(const std::vector<TSimplexID>& shell, SimplexMesh* simplexMesh, bool& status, std::set<TSimplexID>& seen_tet)
{
  std::set<TInt> ans{};
  for(auto const tet : shell)
  {
    if(seen_tet.find(tet) != seen_tet.end())
    {
      status = true;
      break;
    }
    const std::vector<TSimplexID> nodes = SimplicesCell(simplexMesh, tet).getNodes();
    for(auto const n : nodes)
    {
      ans.insert(n);
    }
  }

  return ans;
}
/*----------------------------------------------------------------------------*/
std::multimap<TInt, TInt> sortEdgeNearSurface(SimplexMesh* simplexMesh)
{
  std::multimap<TInt, TInt> sortedEdges{};
  const gmds::BitVector& tet_bitVector = simplexMesh->getBitVectorTet();
  for(unsigned int t = 0 ; t < tet_bitVector.capacity() ; t++)
  {
    if(tet_bitVector[t] != 0)
    {
      for(unsigned int f = 0 ; f < 4 ; f++)
      {
        TSimplexID neighborSimplex = SimplicesCell(simplexMesh, t).oppositeTetraIdx(f);
        if(neighborSimplex < 0)
        {
          std::vector<TInt> face = SimplicesCell(simplexMesh, t).getOrderedFace(f);
          std::pair<TInt, TInt> p01{std::min(face[0], face[1]), std::max(face[0], face[1])};
          std::pair<TInt, TInt> p12{std::min(face[1], face[2]), std::max(face[1], face[2])};
          std::pair<TInt, TInt> p20{std::min(face[2], face[0]), std::max(face[2], face[0])};

          sortedEdges.insert(p01);
          sortedEdges.insert(p12);
          sortedEdges.insert(p20);
        }
      }
    }
  }
  return sortedEdges;
}
/*----------------------------------------------------------------------------*/
double computeTheScaledJacobianOfTheHull(const std::vector<TSimplexID>& shell, SimplexMesh* simplexMesh, std::vector<TSimplexID>& hexNode, std::multimap<TInt, TInt>& mmap_edges, double sJ)
{
  std::vector<std::vector<TInt>> faces{};
  //reorder the node in s to form an hex
  for(auto const tet : shell)
  {
    const std::vector<TSimplexID> adjTets = SimplicesCell(simplexMesh, tet).adjacentTetra();
    for(unsigned int lface = 0 ; lface < 4 ; lface++)
    {
      TSimplexID adjTet = adjTets[lface];
      if(std::find(shell.begin(), shell.end(), adjTet) == shell.end())
      {
        const std::vector<TInt> orderedFace = SimplicesCell(simplexMesh, tet).getOrderedFace(lface);
        faces.push_back(orderedFace);
      }
    }
  }

  TInt A = std::numeric_limits<TInt>::min();
  TInt B = std::numeric_limits<TInt>::min();
  TInt C = std::numeric_limits<TInt>::min();
  TInt D = std::numeric_limits<TInt>::min();

  TInt E = std::numeric_limits<TInt>::min();
  TInt F = std::numeric_limits<TInt>::min();
  TInt G = std::numeric_limits<TInt>::min();
  TInt H = std::numeric_limits<TInt>::min();

  //std::cout << "faces.size()" << " -> " << faces.size() << std::endl;
  std::vector<std::vector<TInt>> ABCD_ALL{};
  std::vector<std::vector<TInt>> EFGH_ALL{};
  if(faces.size() == 12)
  {
    std::set<std::vector<TInt>> seenFaces{};
    for(unsigned int firstFaceId = 0 ; firstFaceId < faces.size() ; firstFaceId++)
    {
      std::vector<TInt> firstFace = faces[firstFaceId];
      std::sort(firstFace.begin(), firstFace.end());
      std::vector<TInt> h_{};
      //std::cout << "    firstFace -> " << firstFace[0] << " | " << firstFace[1] << " | " << firstFace[2] << std::endl;
      //we find the first quad with the faces
      for(unsigned int lface = 0; lface < faces.size() ; lface++) // we ignore the first face
      {
        if(lface == firstFaceId)
          continue;
        std::vector<TInt> faceToCompare = faces[lface];
        std::sort(faceToCompare.begin(), faceToCompare.end());
        //std::cout << "    faceToCompare -> " << faceToCompare[0] << " | " << faceToCompare[1] << " | " << faceToCompare[2] << std::endl;
        std::vector<TInt> inter{};

        std::set_intersection(firstFace.begin(), firstFace.end(),
                              faceToCompare.begin(), faceToCompare.end(),
                              std::back_inserter(inter));

        //std::cout << "        inter.size() -> " << inter.size() << std::endl;
        if(inter.size() == 2)
        {
          //test pour voir si ce sont les 2 meilleurs faces à choisir
          std::set<TInt> s_test_face;
          s_test_face.insert(firstFace[0]) ; s_test_face.insert(firstFace[1]) ; s_test_face.insert(firstFace[2]);
          s_test_face.insert(faceToCompare[0]) ; s_test_face.insert(faceToCompare[1]) ; s_test_face.insert(faceToCompare[2]);
          std::vector<TInt> v_test_face(s_test_face.begin(), s_test_face.end());

          unsigned int cpt_complementary_face = 0;
          for(auto & face : faces)
          {
            std::sort(face.begin(), face.end());
            std::vector<TInt> inter{};
            std::set_intersection(face.begin(), face.end(),
                                  v_test_face.begin(), v_test_face.end(),
                                  std::back_inserter(inter));
            if(inter.size() == 0)
            {
              cpt_complementary_face++;
            }
          }

          if(cpt_complementary_face < 2)
            continue;
          firstFace.erase(std::remove(firstFace.begin(), firstFace.end(), inter.front()), firstFace.end());
          firstFace.erase(std::remove(firstFace.begin(), firstFace.end(), inter.back()), firstFace.end());
          faceToCompare.erase(std::remove(faceToCompare.begin(), faceToCompare.end(), inter.front()), faceToCompare.end());
          faceToCompare.erase(std::remove(faceToCompare.begin(), faceToCompare.end(), inter.back()), faceToCompare.end());

          if(firstFace.size() == 1 && faceToCompare.size() == 1)
          {
            A = inter.front();
            B = firstFace.front();
            C = inter.back();
            D = faceToCompare.front();
            std::vector<TSimplexID> sortedABCD{A,B,C,D};
            std::sort(sortedABCD.begin(), sortedABCD.end());
            if(seenFaces.find(sortedABCD) == seenFaces.end()){
              ABCD_ALL.push_back({A, B, C, D});
              seenFaces.insert(sortedABCD);
            }
          }
        }
      }
      if(ABCD_ALL.size() != 0)
        break;
    }

    for(auto const & ABCD_ : ABCD_ALL)
    {
      std::vector<unsigned int> compFacesId{};
      A = ABCD_[0] ; B = ABCD_[1];
      C = ABCD_[2] ; D = ABCD_[3];

      if(A != std::numeric_limits<TInt>::min() && B != std::numeric_limits<TInt>::min() &&
          C != std::numeric_limits<TInt>::min() && D != std::numeric_limits<TInt>::min())
      {
        std::vector<TInt> sortedABCD{A,B,C,D};
        std::sort(sortedABCD.begin(), sortedABCD.end());
        for(unsigned int lface = 0; lface < faces.size() ; lface++)
        {
          std::vector<TInt> face0 = faces[lface];
          std::sort(face0.begin(), face0.end());
          std::vector<TInt> inter{};

          std::set_intersection(face0.begin(), face0.end(),
                                sortedABCD.begin(), sortedABCD.end(),
                                std::back_inserter(inter));

          if(inter.size() == 0)
          {
            compFacesId.push_back(lface);
          }
        }

        if(compFacesId.size() == 2)
        {
          //compute edge intersection
          std::vector<TInt> f0 = faces[compFacesId.front()];
          std::vector<TInt> f1 = faces[compFacesId.back()];
          std::sort(f0.begin(), f0.end());
          std::sort(f1.begin(), f1.end());
          std::vector<TInt> inter{};

          std::set_intersection(f0.begin(), f0.end(),
                                f1.begin(), f1.end(),
                                std::back_inserter(inter));


          if(inter.size() == 2)
          {
            f0.erase(std::remove(f0.begin(), f0.end(), inter.front()), f0.end());
            f0.erase(std::remove(f0.begin(), f0.end(), inter.back()), f0.end());
            f1.erase(std::remove(f1.begin(), f1.end(), inter.front()), f1.end());
            f1.erase(std::remove(f1.begin(), f1.end(), inter.back()), f1.end());

            if(f0.size() == 1 && f1.size() == 1)
            {
              E = inter.front();
              F = f0.front();
              G = inter.back();
              H = f1.front();
              EFGH_ALL.push_back({E, F, G, H});
            }
          }
        }
      }
    }

    if(ABCD_ALL.size() == EFGH_ALL.size() /*&& ABCD_ALL.size() == 6*/)
    {
      for(unsigned int i = 0 ; i < ABCD_ALL.size() ; i++)
      {
        std::vector<TInt> ABCD = ABCD_ALL[i];
        std::vector<TInt> EFGH = EFGH_ALL[i];
        TInt A = ABCD[0] ; TInt B = ABCD[1] ;
        TInt C = ABCD[2] ; TInt D = ABCD[3] ;

        TInt E = EFGH[0] ; TInt F = EFGH[1] ;
        TInt G = EFGH[2] ; TInt H = EFGH[3] ;

        std::vector<TInt> ADCB{A,D,C,B};
        std::vector<TInt> EHGF{E,H,G,F};

        if(A != std::numeric_limits<TInt>::min() && B != std::numeric_limits<TInt>::min() &&
            C != std::numeric_limits<TInt>::min() && D != std::numeric_limits<TInt>::min() &&
              E != std::numeric_limits<TInt>::min() && F != std::numeric_limits<TInt>::min() &&
                G != std::numeric_limits<TInt>::min() && G != std::numeric_limits<TInt>::min())
        {
          std::vector<std::vector<TInt>> combinaisonFirstFaces{ABCD, ADCB};

          std::vector<std::vector<TInt>> combinaisonSecondFaces{};
          for(unsigned j = 0 ; j < 4 ; j++)
          {
              std::vector<TInt> h_{};
              h_.push_back(EFGH[j]); h_.push_back(EFGH[(j+1)%4]);
              h_.push_back(EFGH[(j+2)%4]); h_.push_back(EFGH[(j+3)%4]);
              combinaisonSecondFaces.push_back(h_);
          }
          for(unsigned j = 0 ; j < 4 ; j++)
          {
              std::vector<TInt> h_{};
              h_.push_back(EHGF[j]); h_.push_back(EHGF[(j+1)%4]);
              h_.push_back(EHGF[(j+2)%4]); h_.push_back(EHGF[(j+3)%4]);
              combinaisonSecondFaces.push_back(h_);
          }
          std::vector<TInt> goodhexahedron{};

          for(auto const & combiansionFirstFace : combinaisonFirstFaces)
          {
            for(auto const & combiansionSecondFace : combinaisonSecondFaces)
            {
                math::Point p0 = SimplicesNode(simplexMesh, combiansionFirstFace[0]).getCoords();
                math::Point p1 = SimplicesNode(simplexMesh, combiansionFirstFace[1]).getCoords();
                math::Point p2 = SimplicesNode(simplexMesh, combiansionFirstFace[2]).getCoords();
                math::Point p3 = SimplicesNode(simplexMesh, combiansionFirstFace[3]).getCoords();
                math::Point p4 = SimplicesNode(simplexMesh, combiansionSecondFace[0]).getCoords();
                math::Point p5 = SimplicesNode(simplexMesh, combiansionSecondFace[1]).getCoords();
                math::Point p6 = SimplicesNode(simplexMesh, combiansionSecondFace[2]).getCoords();
                math::Point p7 = SimplicesNode(simplexMesh, combiansionSecondFace[3]).getCoords();
                const Hexahedron hexahedron_(p0, p1, p2, p3, p4, p5, p6, p7);
                double sJ_tmp = hexahedron_.computeScaledJacobian();
                //std::cout << "    sJ_tmp -> " << sJ_tmp << std::endl;
                if(sJ_tmp >= sJ)
                {
                  /*
                  *			      7--------------6
                  * 			   /|		  				/|
                  * 			  / |	  				 / |
                  * 			 4--------------5  |
                  * 			 |  |						|  |
                  * 			 |  3-----------|--2
                  * 			 | /  					| /
                  * 			 |/							|/
                  * 			 0--------------1
                  */
                  //for paraview hew Representation
                  //check if those edges exist :
                  //0-5
                  //1-6
                  //2-7
                  //0-7

                  /*bool flag_edge05 = false;
                  bool flag_edge16 = false;
                  bool flag_edge27 = false;
                  bool flag_edge07 = false;
                  std::pair<TInt, TInt> p05{std::min(combiansionFirstFace[0], combiansionSecondFace[1]), std::max(combiansionFirstFace[0], combiansionSecondFace[1])};
                  std::pair<TInt, TInt> p16{std::min(combiansionFirstFace[1], combiansionSecondFace[2]), std::max(combiansionFirstFace[1], combiansionSecondFace[2])};
                  std::pair<TInt, TInt> p27{std::min(combiansionFirstFace[2], combiansionSecondFace[3]), std::max(combiansionFirstFace[2], combiansionSecondFace[3])};
                  std::pair<TInt, TInt> p07{std::min(combiansionFirstFace[0], combiansionSecondFace[3]), std::max(combiansionFirstFace[0], combiansionSecondFace[3])};



                  auto range = mmap_edges.equal_range(p05.first);
                  for (auto it = range.first; it != range.second; ++it){
                      if(it->second == p05.second)
                      {
                        flag_edge05 = true;
                        break;
                      }
                  }

                  range = mmap_edges.equal_range(p16.first);
                  for (auto it = range.first; it != range.second; ++it){
                      if(it->second == p16.second)
                      {
                        flag_edge16 = true;
                        break;
                      }
                  }

                  range = mmap_edges.equal_range(p27.first);
                  for (auto it = range.first; it != range.second; ++it){
                      if(it->second == p27.second)
                      {
                        flag_edge27 = true;
                        break;
                      }
                  }


                  range = mmap_edges.equal_range(p07.first);
                  for (auto it = range.first; it != range.second; ++it){
                      if(it->second == p07.second)
                      {
                        flag_edge07 = true;
                        break;
                      }
                  }*/

                  //if(flag_edge05 && flag_edge16 && flag_edge27 && flag_edge07)
                  {
                    /*std::cout << "ABCD -> " << A << " | " << B << " | " << C << " | " << D << std::endl;
                    std::cout << "EFGH -> " << E << " | " << F << " | " << G << " | " << H << std::endl;
                    std::cout << "p05 -> " << p05.first << " | " << p05.second << std::endl;
                    std::cout << "p16 -> " << p16.first << " | " << p16.second << std::endl;
                    std::cout << "p27 -> " << p27.first << " | " << p27.second << std::endl;
                    std::cout << "p07 -> " << p07.first << " | " << p07.second << std::endl;
                    std::cout << "flag_edge05 -> " << flag_edge05 << std::endl;
                    std::cout << "flag_edge16 -> " << flag_edge16 << std::endl;
                    std::cout << "flag_edge27 -> " << flag_edge27 << std::endl;
                    std::cout << "flag_edge07 -> " << flag_edge07 << std::endl;
                    std::cout << std::endl;*/
                    sJ = sJ_tmp;
                    goodhexahedron.clear();
                    //we add the first face
                    goodhexahedron.push_back(combiansionFirstFace[0]) ; goodhexahedron.push_back(combiansionFirstFace[1]) ;
                    goodhexahedron.push_back(combiansionFirstFace[2]) ; goodhexahedron.push_back(combiansionFirstFace[3]) ;
                    //we add the second face
                    goodhexahedron.push_back(combiansionSecondFace[0]) ; goodhexahedron.push_back(combiansionSecondFace[1]) ;
                    goodhexahedron.push_back(combiansionSecondFace[2]) ; goodhexahedron.push_back(combiansionSecondFace[3]) ;
                  }
                }
            }
          }
          if(goodhexahedron.size() == 8)
            hexNode = goodhexahedron;
        }
      }
    }
  }
  return sJ;
}
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/*for(unsigned int t = 0 ; t < tet_bit.capacity() ; t++)
{
  if(tet_bit[t] != 0 && seen.find(t) == seen.end())
  {
    std::vector<TInt> nodes = SimplicesCell(&simplexMesh, t).getNodes();
    std::vector<double> jacobians = SimplicesCell(&simplexMesh, t).getBestNodeJacobianNormalise();
    unsigned int lN_res = 0;
    double jacobian_tmp = std::numeric_limits<double>::min();
    for(unsigned int lN = 0 ; lN < 4 ; lN++)
    {
      double j = jacobians[lN];
      if(j > jacobian_tmp)
      {
        jacobian_tmp = j;
        lN_res = lN;
      }
    }

    if(jacobian_tmp < 0.70)
      continue;
    unsigned int lN1 = (lN_res+1)%4;
    unsigned int lN2 = (lN_res+2)%4;
    unsigned int lN3 = (lN_res+3)%4;

    TInt A = nodes[lN_res];
    TInt B = nodes[lN1];
    TInt C = nodes[lN2];
    TInt D = nodes[lN3];

    TInt F = std::numeric_limits<TSimplexID>::min();
    TInt E = std::numeric_limits<TSimplexID>::min();
    TInt H = std::numeric_limits<TSimplexID>::min();
    TInt G = std::numeric_limits<TSimplexID>::min();
    std::vector<TSimplexID> shelBC = SimplicesNode(&simplexMesh, B).shell(SimplicesNode(&simplexMesh, C));
    std::vector<TSimplexID> shelCD = SimplicesNode(&simplexMesh, C).shell(SimplicesNode(&simplexMesh, D));
    std::vector<TSimplexID> shelDB = SimplicesNode(&simplexMesh, D).shell(SimplicesNode(&simplexMesh, B));

    bool flag = false;
    TSimplexID t0_ = std::numeric_limits<TSimplexID>::min();
    TSimplexID t1_ = std::numeric_limits<TSimplexID>::min();
    TSimplexID t2_ = std::numeric_limits<TSimplexID>::min();
    std::vector<TInt> nn{};
    for(auto const t0 : shelBC)
    {
      if(seen.find(t0) != seen.end())
        continue;
      for(auto const t1 : shelCD)
      {
        if(seen.find(t1) != seen.end())
          continue;
        std::vector<TInt> commonNodes01 = SimplicesCell(&simplexMesh, t0).commonNode(SimplicesCell(&simplexMesh, t1));
        commonNodes01.erase(std::remove(commonNodes01.begin(), commonNodes01.end(), C), commonNodes01.end());
        if(commonNodes01.size() == 1)
        {
          for(auto const t2 : shelDB)
          {
            if(seen.find(t2) != seen.end())
              continue;
            std::vector<TInt> commonNodes02 = SimplicesCell(&simplexMesh, t0).commonNode(SimplicesCell(&simplexMesh, t2));
            std::vector<TInt> commonNodes12 = SimplicesCell(&simplexMesh, t1).commonNode(SimplicesCell(&simplexMesh, t2));
            commonNodes02.erase(std::remove(commonNodes02.begin(), commonNodes02.end(), B), commonNodes02.end());
            commonNodes12.erase(std::remove(commonNodes12.begin(), commonNodes12.end(), D), commonNodes12.end());
            if(commonNodes02.size() == 1 && commonNodes12.size() == 1)
            {
              if(commonNodes02.front() == commonNodes01.front())
              {
                G = commonNodes02.front();
                //compute the nodes of t0, t1 and t2 that are not shown in the parameters
                //to :
                std::vector<TInt> nodes0 = SimplicesCell(&simplexMesh, t0).getNodes();
                nodes0.erase(std::remove_if(nodes0.begin(), nodes0.end(), [=](TInt node){
                return (node == G || node == B || node == C);
                }), nodes0.end());
                if(nodes0.size() == 1)
                  F = nodes0.front();


                //t1
                std::vector<TInt> nodes1 = SimplicesCell(&simplexMesh, t1).getNodes();
                nodes1.erase(std::remove_if(nodes1.begin(), nodes1.end(), [=](TInt node){
                return (node == G || node == D || node == C);
                }), nodes1.end());
                if(nodes1.size() == 1)
                  E = nodes1.front();

                //t2
                std::vector<TInt> nodes2 = SimplicesCell(&simplexMesh, t2).getNodes();
                nodes2.erase(std::remove_if(nodes2.begin(), nodes2.end(), [=](TInt node){
                return (node == G || node == D || node == B);
                }), nodes2.end());
                if(nodes2.size() == 1)
                  H = nodes2.front();

              if(A != std::numeric_limits<TInt>::min() && B != std::numeric_limits<TInt>::min() && F != std::numeric_limits<TInt>::min() &&
                C != std::numeric_limits<TInt>::min() && D != std::numeric_limits<TInt>::min() && H != std::numeric_limits<TInt>::min() &&
                G != std::numeric_limits<TInt>::min() && E != std::numeric_limits<TInt>::min())
              {

                //check if the connection between face0 : (ABFC) and face1 (DHGE) exist
                std::vector<TSimplexID> shellAD = SimplicesNode(&simplexMesh, A).shell(SimplicesNode(&simplexMesh, D));
                if(shellAD.size() == 0)
                  continue;

                std::vector<TSimplexID> shellBH = SimplicesNode(&simplexMesh, B).shell(SimplicesNode(&simplexMesh, H));
                if(shellBH.size() == 0)
                  continue;

                std::vector<TSimplexID> shellFG = SimplicesNode(&simplexMesh, F).shell(SimplicesNode(&simplexMesh, G));
                if(shellFG.size() == 0)
                  continue;

                std::vector<TSimplexID> shellCE = SimplicesNode(&simplexMesh, C).shell(SimplicesNode(&simplexMesh, E));
                if(shellCE.size() == 0)
                  continue;

                  std::set<TInt> s{};
                  s.insert(A); s.insert(B); s.insert(C); s.insert(D);
                  s.insert(E); s.insert(F); s.insert(G); s.insert(H);
                  if(s.size() != 8)
                    continue;
                  std::vector<TInt> h{A, B, F,C, D, H, G, E};

                  math::Point p0 = SimplicesNode(&simplexMesh, h[0]).getCoords();
                  math::Point p1 = SimplicesNode(&simplexMesh, h[1]).getCoords();
                  math::Point p2 = SimplicesNode(&simplexMesh, h[2]).getCoords();
                  math::Point p3 = SimplicesNode(&simplexMesh, h[3]).getCoords();
                  math::Point p4 = SimplicesNode(&simplexMesh, h[4]).getCoords();
                  math::Point p5 = SimplicesNode(&simplexMesh, h[5]).getCoords();
                  math::Point p6 = SimplicesNode(&simplexMesh, h[6]).getCoords();
                  math::Point p7 = SimplicesNode(&simplexMesh, h[7]).getCoords();
                  const Hexahedron hexahedron_(p0, p1, p2, p3, p4, p5, p6, p7);
                  if(hexahedron_.computeScaledJacobian() < 0.5)
                    continue;

                 hexs.push_back(h);

                  //std::vector<TInt> hex{};
                  //for(auto const n : h)
                  //{
                      //const math::Point p = SimplicesNode(&simplexMesh, n).getCoords();
                      //TInt id = testMesh.addNode(p);
                      //hex.push_back(id);
                  //}
                  //hexs.push_back(hex);
                  seen.insert(t);
                  seen.insert(t0);
                  seen.insert(t1);
                  seen.insert(t2);
                  markedTet.unselect(t);
                  markedTet.unselect(t0);
                  markedTet.unselect(t1);
                  markedTet.unselect(t2);
                }
                flag = true;
                break;
              }
            }
          }
          if(flag)
            break;
        }
      }
      if(flag)
        break;
    }
    //if(hexs.size() == 10)
      //break;
  }
}

Mesh m(MeshModel(DIM3 | R | F | E | N |
  R2N | F2N | E2N | R2F | F2R |
  F2E | E2F | R2E | N2R | N2F | N2E));

  std::unordered_map<TInt, Node> map{};
  const gmds::BitVector & nodeBitVector = simplexMesh.getBitVectorNodes();
  for(unsigned int n = 0 ; n < nodeBitVector.capacity() ; n++)
  {
    if(nodeBitVector[n] != 0)
    {
      const Node nodeHEX =  m.newNode(SimplicesNode(&simplexMesh, n).getCoords());
      map[n] = nodeHEX;
    }
  }
for(auto const & hex : hexs){
  m.newHex(map[hex[0]], map[hex[1]], map[hex[2]], map[hex[3]], map[hex[4]], map[hex[5]], map[hex[6]], map[hex[7]]);
}

gmds::IGMeshIOService ioServiceMESH(&m);
gmds::VTKWriter vtkWriterMESH(&ioServiceMESH);
vtkWriterMESH.setCellOptions(gmds::N|gmds::R|gmds::F);
//vtkWriterMESH.setDataOptions(gmds::N|gmds::R|gmds::F);
vtkWriterMESH.write("HEX_ONLY_" + fOut);
std::cout << "HEX CREATED " << std::endl;

std::cout << "hexs size -> " << hexs.size() << std::endl;
std::cout << "market tet -> " << markedTet.size() << std::endl;
//testMesh.setHexadronData(hexs);
//simplexMesh.setMarkedTet(markedTet);
//gmds::ISimplexMeshIOService ioServiceTestMesh(&testMesh);
//gmds::VTKWriter vtkWriterTestMesh(&ioServiceTestMesh);
//vtkWriterTestMesh.setCellOptions(gmds::N|gmds::R);
//vtkWriterTestMesh.write(fOut);
simplexMesh.setHexadronData(hexs);
//simplexMesh.setMarkedTet(markedTet);
gmds::VTKWriter vtkWriter(&ioService);
vtkWriter.setCellOptions(gmds::N|gmds::R);
vtkWriter.write(fOut);
*/
