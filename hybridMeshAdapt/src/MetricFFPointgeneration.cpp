#include "gmds/hybridMeshAdapt/MetricFFPointgeneration.h"
#include "gmds/hybridMeshAdapt/SimplexMesh.h"
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
#include <gmds/io/VTKWriter.h>
/*----------------------------------------------------------------------------*/
#include <unordered_set>
#include <deque>
#include <ctime>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace math;
using namespace hybrid;
using namespace operators;
using namespace simplicesNode;
using namespace simplicesTriangle;
using namespace simplicesCell;
/*----------------------------------------------------------------------------*/
MetricFFPointgeneration::MetricFFPointgeneration(SimplexMesh* simplexMesh):m_simplexMesh(simplexMesh),m_oc(Octree(simplexMesh, 10))
{
  //m_nodesMesh.newVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  m_nodesMesh.newVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  m_nodesMesh.newVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
  m_nodesMesh.newVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
}
/*----------------------------------------------------------------------------*/
MetricFFPointgeneration::~MetricFFPointgeneration()
{

}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::execute()
{
  std::clock_t start = std::clock();
  gmds::Variable<int>* BND_SURFACE_COLOR = nullptr;

  std::vector<double> simplexMesh_Borders = m_oc.getBorderOctree();
  Octree* simplexNodes_Octree = new Octree(&m_nodesMesh, 3,
                                    simplexMesh_Borders[0],simplexMesh_Borders[1],
                                    simplexMesh_Borders[2],simplexMesh_Borders[3],
                                    simplexMesh_Borders[4],simplexMesh_Borders[5]);

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  Mesh m0(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N | R2F
                        | F2R |F2E | E2F | R2E | N2R | N2F | N2E));

  m_nodesMesh.setOctree(simplexNodes_Octree);
  std::vector<std::vector<Node>> nodes{};
  m_nodesMesh.getOctree()->writeOctree(m0, nodes);
  for(auto const & node : nodes)
    m0.newHex(node[0], node[1], node[2], node[3], node[4], node[5], node[6], node[7]);


  gmds::IGMeshIOService ioServiceM0(&m0);
  gmds::VTKWriter vtkWriterM0(&ioServiceM0);
	vtkWriterM0.setCellOptions(gmds::N|gmds::R);
	vtkWriterM0.setDataOptions(gmds::N|gmds::R);
	vtkWriterM0.write("Octree_Nodes.vtk");
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  try{
    BND_SURFACE_COLOR = m_nodesMesh.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  std::vector<double> edges_length{};
  const std::map<unsigned int, std::vector<TInt>> sortedEdges = buildSortedEdges();
  const std::vector<std::vector<double>> edgesU = buildParamEdgeU(sortedEdges, edges_length);


  std::vector<TInt> nodeAdded{};
  unsigned int i = 0;
  for(auto const & sortedEdge : sortedEdges)
  {
    //if(i == 8)
    {
      const unsigned int edgeId = sortedEdge.first;
      std::vector<TInt> edge = sortedEdge.second;
      std::vector<double> edgeU = edgesU[i];
      double edge_length = edges_length[i];
      subdivideEdgeUsingMetric_Relaxation(nodeAdded, edge, edgeU, edge_length, edgeId);
    }
    i++;
  }

  std::cout << "EDGE NODES CREATED " << std::endl;
  for(unsigned int n = 0 ; n < m_nodesMesh.getBitVectorNodes().capacity() ; n++)
  {
    if(m_nodesMesh.getBitVectorNodes()[n] != 0)
    {
      m_nodesMesh.addTetraedre(n,n,n,n);
    }
  }
  gmds::ISimplexMeshIOService ioServiceEDGE(&m_nodesMesh);
  gmds::VTKWriter vtkWriterEGE(&ioServiceEDGE);
  vtkWriterEGE.setCellOptions(gmds::N|gmds::R);
  vtkWriterEGE.setDataOptions(gmds::N|gmds::R);
  vtkWriterEGE.write("metricFF_EDGE.vtk");

  while(nodeAdded.size() != 0)
  {
    std::cout << "    nodeAdded.size() -> " << nodeAdded.size() << std::endl;
    nodesSpreading(nodeAdded, true);
    //break;
  }
  for(unsigned int n = 0 ; n < m_nodesMesh.getBitVectorNodes().capacity() ; n++)
  {
    if(m_nodesMesh.getBitVectorNodes()[n] != 0)
    {
      m_nodesMesh.addTetraedre(n,n,n,n);
    }
  }
  gmds::ISimplexMeshIOService ioServiceSURFACE(&m_nodesMesh);
  gmds::VTKWriter vtkWriterSURFACE(&ioServiceSURFACE);
  vtkWriterSURFACE.setCellOptions(gmds::N|gmds::R);
  vtkWriterSURFACE.setDataOptions(gmds::N|gmds::R);
  vtkWriterSURFACE.write("metricFF_SURFACE.vtk");
  std::cout << "SURFACE NODES CREATED " << std::endl;
  throw gmds::GMDSException("END");
  const gmds::BitVector& m_nodes = m_nodesMesh.getBitVectorNodes();
  for(unsigned int n = 0 ; n < m_nodes.capacity() ; n++)
  {
    if(m_nodes[n] != 0)
    {
      if((*BND_SURFACE_COLOR)[n] != 0)
      {
        nodeAdded.push_back(n);
      }
    }
  }

  while(nodeAdded.size() != 0)
  {
    std::cout << "    nodeAdded.size() -> " << nodeAdded.size() << std::endl;
    nodesSpreading(nodeAdded);
  }
  gmds::ISimplexMeshIOService ioServiceVOLUME(&m_nodesMesh);
  gmds::VTKWriter vtkWriterVOLUME(&ioServiceVOLUME);
  vtkWriterVOLUME.setCellOptions(gmds::N|gmds::R);
  vtkWriterVOLUME.setDataOptions(gmds::N|gmds::R);
  vtkWriterVOLUME.write("metricFF_VOLUME.vtk");
  std::cout << "VOLUME NODES CREATED " << std::endl;

  processNodesStructure();
  std::set<std::vector<TInt>> hexs{};
  computeHexa(hexs);

  std::cout << "hex size -> " << hexs.size() << std::endl;

  Mesh m(MeshModel(DIM3 | R | F | E | N |
                   R2N | F2N | E2N | R2F | F2R |
                   F2E | E2F | R2E | N2R | N2F | N2E));


  const gmds::BitVector& nodesMeshBitVector = m_nodesMesh.getBitVectorNodes();
  std::vector<Node> seenNodes(nodesMeshBitVector.capacity(), Node());
  for(auto const & hex : hexs)
  {
    std::vector<Node> nodesHex{};
    for(unsigned int n = 0 ; n < 8 ; n++)
    {
      if(seenNodes[hex[n]] == Node())
      {
        const Node n0 = m.newNode(SimplicesNode(&m_nodesMesh, hex[n]).getCoords());
        seenNodes[hex[n]] = n0;
        nodesHex.push_back(n0);
      }
      else
        nodesHex.push_back(seenNodes[hex[n]]);
    }
    m.newHex(nodesHex[0], nodesHex[1], nodesHex[2], nodesHex[3], nodesHex[4], nodesHex[5], nodesHex[6], nodesHex[7]);
  }

  gmds::IGMeshIOService ioService(&m);
  gmds::VTKWriter vtkWriter(&ioService);
  vtkWriter.setCellOptions(gmds::N|gmds::R);
  vtkWriter.setDataOptions(gmds::N|gmds::R);
  vtkWriter.write("HEX.vtk");

  for(auto const & d : m_nodeStructure)
  {
    TInt node = d.first;
    for(auto const n : d.second)
    {
      if(n != -1)
      {
        m_nodesMesh.addTriangle(d.first, n, n);
      }
    }
  }

  gmds::ISimplexMeshIOService ioServiceGRID(&m_nodesMesh);
  gmds::VTKWriter vtkWriterGRID(&ioServiceGRID);
  vtkWriterGRID.setCellOptions(gmds::N|gmds::F);
  vtkWriterGRID.setDataOptions(gmds::N|gmds::F);
  vtkWriterGRID.write("metricFF_Grid.vtk");

  double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  std::cout << "DURATION -> " << duration << std::endl;
}
/*----------------------------------------------------------------------------*/
std::vector<TInt> MetricFFPointgeneration::findBoundedNode(const double t, const std::vector<TInt>& edgeNodes) const
{
  if(t > 1 || t < 0)
    throw gmds::GMDSException("t > 1 || t < 0");

  unsigned int n = edgeNodes.size();
  std::vector<TInt> ans{};
  for(unsigned int i = 0 ; i < n-1; i++)
  {
    double boundDown = static_cast<double>(i) / static_cast<double>(n-1);
    double boudUp = static_cast<double>(i+1) / static_cast<double>(n-1);

    if(t >= boundDown && t <= boudUp)
    {
      if(i != 0)
        ans.push_back(i-1);

      ans.push_back(i);
      ans.push_back(i+1);

      if(i+1 != edgeNodes.size() - 1)
        ans.push_back(i+2);

      break;
    }
  }
  return ans;
}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::correctNodeLabel()
{
  gmds::Variable<int>* BND_SURFACE_COLOR    = nullptr;
  gmds::Variable<int>* BND_CURVE_COLOR    = nullptr;

  try{
    BND_SURFACE_COLOR = m_nodesMesh.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
    BND_CURVE_COLOR = m_nodesMesh.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  const gmds::BitVector& m_nodes = m_nodesMesh.getBitVectorNodes();
  for(unsigned int n = 0 ; n < m_nodes.capacity() ; n++)
  {
    if(m_nodes[n] != 0)
    {
      if((*BND_SURFACE_COLOR)[n] == 0 || (*BND_CURVE_COLOR)[n] == 0 /*add vertex color*/)
      {
        int surfaceLabel = 0;
        m_simplexMesh->onSurface(SimplicesNode(&m_nodesMesh,n).getCoords(), surfaceLabel);
        if(surfaceLabel != 0)
          m_nodesMesh.deleteNode(n);
      }
    }
  }
}
/*----------------------------------------------------------------------------*/
/*void MetricFFPointgeneration::computeHexa(std::set<std::vector<TInt>> & hexa)
{
  const gmds::BitVector& nodeBitVector = m_nodesMesh.getBitVectorNodes();
  std::set<std::vector<TInt>> faces{};
  unsigned int sizeFACE = 4;
  computeQuadFaces(faces);


  //compute the hull of a node in faces
  std::multimap<TInt, std::vector<TInt>> mm{};
  std::unordered_map<TInt, std::set<TInt>> um{};

  for(auto const face : faces)
  {
    for(unsigned int n = 0 ; n < sizeFACE ; n++)
    {
      std::pair<TInt, std::vector<TInt>> p{face[n], std::vector<TInt>{face[(n + 1) % sizeFACE], face[(n + 2) % sizeFACE], face[(n + 3) % sizeFACE]}};
      mm.insert(p);
    }
  }

  for(unsigned int n = 0 ; n < nodeBitVector.capacity() ; n++)
  {
    if(nodeBitVector[n] != 0)
    {
      std::set<TInt> s{};
      auto r = mm.equal_range(n);
      for(auto it = r.first ; it != r.second ; it++)
      {
        for(auto const neighborNode : it->second)
        {
          s.insert(neighborNode);
        }
      }
      um.insert(std::pair<TInt, std::set<TInt>>{n, s});
    }
  }

  std::vector<TInt> v{};
  std::set<std::vector<TInt>> seenHexa{};
  for(auto const & p0 : um)
  {
    TInt nodeA = p0.first;
    for(auto const adjNode : p0.second)
    {
      std::unordered_set<TInt> seen{};
      for(auto const nodeB : m_nodeStructure[adjNode])
      {
        if(std::find(m_nodeStructure[nodeA].begin(), m_nodeStructure[nodeA].end(), adjNode) == m_nodeStructure[nodeA].end())
        {
          if(seen.find(nodeB) == seen.end())
          {
            std::set<TInt> p1 = um[nodeB];
            seen.insert(nodeB);
            if(std::find(p1.begin(), p1.end(), nodeA) == p1.end())
            {
              v.clear();
              if(p1 != p0.second)
              {
                for(auto const & v0 : p0.second)
                {
                  for(auto const & v1 : p1)
                  {
                    if(v0 == v1 && v0 != nodeA && v0 != nodeB)
                    {
                      v.push_back(v0);
                    }
                  }
                }
              }

              if(v.size() == 6)
              {
                v.push_back(nodeA);
                v.push_back(nodeB);
                std::sort(v.begin(), v.end());
                if(seenHexa.find(v) == seenHexa.end())
                {
                  seenHexa.insert(v);
                  hexa.insert(v);
                }
              }
            }
          }
        }
      }
    }
  }
}*/
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::computeHexa(std::set<std::vector<TInt>> & hexas)
{
  const gmds::BitVector& nodeBitVector = m_nodesMesh.getBitVectorNodes();
  std::set<std::vector<TInt>> faces{};
  unsigned int sizeFACE = 4;

  computeQuadFaces(faces);

  std::unordered_map<TInt, std::vector<std::set<TInt>>> um{};
  std::unordered_map<std::string, std::vector<TInt>> indirect_um{};
  std::unordered_set<std::string> seen{};

  for(auto const face : faces)
  {
    std::set<TInt> s(face.begin(), face.end());
    std::string str{};
    for(auto const nbr : s)
      str += std::to_string(nbr) + "|";

    indirect_um[str] = face;
    for(unsigned int n = 0 ; n < sizeFACE ; n++)
      um[face[n]].push_back(s);
  }

  unsigned int i = 0;
  for(auto const face : faces)
  {
    std::unordered_map<std::string, unsigned int> h{};
    const std::vector<TInt>& neiboor0 = m_nodeStructure[face[0]];
    const std::vector<TInt>& neiboor1 = m_nodeStructure[face[1]];
    const std::vector<TInt>& neiboor2 = m_nodeStructure[face[2]];
    const std::vector<TInt>& neiboor3 = m_nodeStructure[face[3]];

    std::vector<std::vector<TInt>> neigboors{neiboor0, neiboor1, neiboor2, neiboor3};
    for(auto const & neigboor : neigboors)
    {
      for(auto const & n : neigboor)
      {
        if(n != face[0] && n != face[1] && n != face[2] && n != face[3] && n >= 0)
        {
          for(auto const & vec : um[n])
          {
            std::string str{};
            for(auto const nbr : vec)
              str += std::to_string(nbr) + "|";
            h[str]++;
          }
        }
      }
    }

    for(auto const & d : h)
    {
      if(d.second == 4)
      {
        std::vector<TInt> hexa(face.begin(), face.end());
        std::vector<TInt> face1{};
        for(auto const node : face)
        {
          for(auto const neigboor : m_nodeStructure[node])
          {
            if(std::find(indirect_um[d.first].begin(), indirect_um[d.first].end(), neigboor) != indirect_um[d.first].end())
              if(std::find(face1.begin(), face1.end(), neigboor) == face1.end())
                face1.push_back(neigboor);
          }
        }

        if(face1.size() == sizeFACE)
        {
          std::vector<TInt> test_seen(face.begin(), face.end());
          std::copy(face1.begin(), face1.end(), std::back_inserter(test_seen));
          std::sort(test_seen.begin(), test_seen.end());

          std::string str{};
          for(auto const var : test_seen)
            str += std::to_string(var);

          if(seen.find(str) == seen.end())
          {
            std::copy(face1.begin(), face1.end(), std::back_inserter(hexa));
            std::set<TInt> redunancy(hexa.begin(), hexa.end());
            if(redunancy.size() == 8)
            {
              const math::Point p0 = SimplicesNode(&m_nodesMesh, hexa[0]).getCoords();
              const math::Point p1 = SimplicesNode(&m_nodesMesh, hexa[1]).getCoords();
              const math::Point p2 = SimplicesNode(&m_nodesMesh, hexa[2]).getCoords();
              const math::Point p3 = SimplicesNode(&m_nodesMesh, hexa[3]).getCoords();
              const math::Point p4 = SimplicesNode(&m_nodesMesh, hexa[4]).getCoords();
              const math::Point p5 = SimplicesNode(&m_nodesMesh, hexa[5]).getCoords();
              const math::Point p6 = SimplicesNode(&m_nodesMesh, hexa[6]).getCoords();
              const math::Point p7 = SimplicesNode(&m_nodesMesh, hexa[7]).getCoords();

              const Hexahedron h(p0, p1, p2, p3, p4, p5, p6, p7);
              if(h.isValid()/* && h.computeScaledJacobian() > -0.7*/)
                hexas.insert(hexa);

              seen.insert(str);
            }
          }
        }
      }
    }
  }
}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::computeQuadFaces(std::set<std::vector<TInt>> & faces) const
{
  faces.clear();
  std::set<std::set<TInt>> seen{};
  std::multimap<TInt, TInt> edges{};

  for(auto const & s : m_nodeStructure)
  {
    TInt node = s.first;
    for(auto const & n : s.second)
    {
      if(n != -1)
      {
        std::pair<TInt, TInt> p(node,n);
        edges.insert(p);
      }
    }
  }

  for(auto const & edge : edges)
  {
    TInt nodeA = edge.first;
    TInt nodeB = edge.second;

    auto rangeA = edges.equal_range(nodeA);
    auto rangeB = edges.equal_range(nodeB);

    gmds::BitVector bitVectorNodeRangeA(m_nodesMesh.getBitVectorNodes().capacity());
    for (auto i = rangeA.first; i != rangeA.second; ++i)
    {
      if(i->second != nodeB)
      {
        bitVectorNodeRangeA.assign(i->second);
      }
    }

    for (auto j = rangeB.first; j != rangeB.second; ++j)
    {
      if(j->second != nodeA)
      {
        auto rangeC = edges.equal_range(j->second);
        for (auto k = rangeC.first; k != rangeC.second; ++k)
        {
          if(bitVectorNodeRangeA[k->second] != 0)
          {
            std::vector<TInt> v{nodeA, nodeB, j->second, k->second};
            std::set<TInt> s(v.begin(), v.end());
            if(seen.find(s) == seen.end())
            {
              faces.insert(v);
              seen.insert(s);
            }
          }
        }
      }
    }
  }
}
/*----------------------------------------------------------------------------*/
Point MetricFFPointgeneration::computeTheEdgeNodeCoordinate(const double u, const std::vector<TInt>& edge, const std::vector<double>& edgeU) const
{
  Point pt;
  if(edgeU.size() < 2)
  {
    throw gmds::GMDSException("edgeU.size() < 2");
  }
  if(!(u >= edgeU.front() && u <= edgeU.back()))
  {
    std::cout << "edgeU.front() -> " << edgeU.front() << std::endl;
    std::cout << "edgeU.back() -> " << edgeU.back() << std::endl;
    std::cout << "u -> " << u << std::endl;
    throw gmds::GMDSException("!(u >= edgeU.front() && u <= edgeU.back())");
  }

  for(unsigned int i = 0 ; i < edgeU.size() - 1 ; i++)
  {
    const double uA = edgeU[i];
    const double uB = edgeU[i + 1];
    if(uA <= u && uB >= u)
    {
      //interpolation (linear) of the middle position using the interval ua, ub and u = 0.5
      if(uB != uA)
      {
        const TInt nodeA = edge[i];
        const TInt nodeB = edge[i + 1];

        const Point ptA = SimplicesNode(m_simplexMesh, nodeA).getCoords();
        const Point ptB = SimplicesNode(m_simplexMesh, nodeB).getCoords();

        const double t = (u-uA) / (uB-uA);

        std::vector<TInt> newEdge0{};
        std::vector<TInt> newEdge1{};

        std::vector<double> edgeU_0 {};
        std::vector<double> edgeU_1 {};

        if(t == 0.0)
        {
          pt = ptA;
        }
        else if(t == 1.0)
        {
          pt = ptB;
        }
        else
        {
          pt = ptA * (1.0 - t) + ptB * t;
        }
        break;
      }
      else
      {
        throw gmds::GMDSException("uB == uA");
      }
    }
  }

  return pt;
}
/*----------------------------------------------------------------------------*/
Eigen::Matrix3d MetricFFPointgeneration::metricInterpolationWithDistorsion(const Eigen::Matrix3d  & metric, const Eigen::Matrix3d  & frameField) const
{
  auto signLambda = [=](const double a){
    double ans = (a < 0.0)?-1.0:1.0;
    return ans;
  };

  Eigen::Matrix3d ans = Eigen::Matrix3d::Zero();

  double distorsion = computeDistorsionMetric(metric);
  //std::cout << "distorsion -> " << distorsion << std::endl;
  distorsion = (distorsion < 1E-4)? 0.0 : distorsion;
  //distorsion = 0.5;//1E-4;

  /*Eigen::EigenSolver<Eigen::Matrix3d> es_metric(metric);
  Eigen::EigenSolver<Eigen::Matrix3d> es_ff(frameField);

  Eigen::MatrixXcd P_metric = es_metric.eigenvectors();
  Eigen::MatrixXcd D_metric = es_metric.eigenvalues().asDiagonal();
  std::cout << "P_metric : " << std::endl << P_metric << std::endl;
  std::cout << "P_metric.inverse() : " << std::endl << P_metric.inverse() << std::endl;
  std::cout << std::endl;

  Eigen::MatrixXcd P_ff = es_ff.eigenvectors();
  Eigen::MatrixXcd D_ff = es_ff.eigenvalues().asDiagonal();
  std::cout << "P_ff : " << std::endl << P_ff << std::endl;
  std::cout << "P_ff.inverse() : " << std::endl << P_ff.inverse() << std::endl;
  std::cout << std::endl;

  D_metric(0,0) = std::complex<double>(signLambda(D_metric(0,0).real()) * std::pow(std::abs(D_metric(0,0).real()), 1.0 - distorsion), signLambda(D_metric(0,0).imag()) * std::pow(std::abs(D_metric(0,0).imag()), 1.0 - distorsion));
  D_metric(1,1) = std::complex<double>(signLambda(D_metric(1,1).real()) * std::pow(std::abs(D_metric(1,1).real()), 1.0 - distorsion), signLambda(D_metric(1,1).imag()) * std::pow(std::abs(D_metric(1,1).imag()), 1.0 - distorsion));
  D_metric(2,2) = std::complex<double>(signLambda(D_metric(2,2).real()) * std::pow(std::abs(D_metric(2,2).real()), 1.0 - distorsion), signLambda(D_metric(2,2).imag()) * std::pow(std::abs(D_metric(2,2).imag()), 1.0 - distorsion));
  std::cout << "D_metric after: " << std::endl << D_metric << std::endl;
  std::cout << std::endl;

  D_ff(0,0) = std::complex<double>(signLambda(D_ff(0,0).real()) * std::pow(std::abs(D_ff(0,0).real()), distorsion), signLambda(D_ff(0,0).imag()) * std::pow(std::abs(D_ff(0,0).imag()), distorsion));
  D_ff(1,1) = std::complex<double>(signLambda(D_ff(1,1).real()) * std::pow(std::abs(D_ff(1,1).real()), distorsion), signLambda(D_ff(1,1).imag()) * std::pow(std::abs(D_ff(1,1).imag()), distorsion));
  D_ff(2,2) = std::complex<double>(signLambda(D_ff(2,2).real()) * std::pow(std::abs(D_ff(2,2).real()), distorsion), signLambda(D_ff(2,2).imag()) * std::pow(std::abs(D_ff(2,2).imag()), distorsion));
  std::cout << "D_ff after: " << std::endl << D_ff << std::endl;


  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "(P_metric.inverse() * D_metric * P_metric) : " << std::endl << (P_metric.inverse() * D_metric * P_metric) << std::endl;
  std::cout << std::endl;
  std::cout << "(P_ff * D_ff * P_ff.inverse()) : " << std::endl << (P_ff * D_ff * P_ff.inverse()) << std::endl;
  //Eigen::Matrix3cd t = (P_metric.inverse() * D_metric * P_metric) * (P_ff.inverse() * D_ff * P_ff);*/

  Eigen::Vector3d col0 = frameField.col(0).real();
  Eigen::Vector3d col1 = frameField.col(1).real();
  Eigen::Vector3d col2 = frameField.col(2).real();

  double a = sqrt(col0.transpose() * metric * col0);
  double b = sqrt(col1.transpose() * metric * col1);
  double c = sqrt(col2.transpose() * metric * col2);

  Eigen::Matrix3d scalarReduction = Eigen::Matrix3d::Identity();
  scalarReduction(0,0) = 1.0/a;
  scalarReduction(1,1) = 1.0/b;
  scalarReduction(2,2) = 1.0/c;

  Eigen::Matrix3cd t = (1-distorsion)*metric + distorsion*scalarReduction*frameField;
  //std::cout << "t : " << std::endl << t << std::endl;
  //std::cout << std::endl;
  //std::cout << "t.real() : " << std::endl << t.real() << std::endl;

  return t.real();

}
/*----------------------------------------------------------------------------*/
double MetricFFPointgeneration::computeDistorsionMetric(const Eigen::Matrix3d & m) const
{
  Eigen::Matrix3cd m_;
  m_ << std::complex<double>(m(0,0), 0), std::complex<double>(m(0,1), 0), std::complex<double>(m(0,2), 0),
  std::complex<double>(m(1,0), 0), std::complex<double>(m(1,1), 0), std::complex<double>(m(1,2), 0),
  std::complex<double>(m(2,0), 0), std::complex<double>(m(2,1), 0), std::complex<double>(m(2,2), 0);


  Eigen::EigenSolver<Eigen::Matrix3d> es(m);
  Eigen::MatrixXcd L_ = es.eigenvalues();

  const double lambda0 = L_(0).real();
  const double lambda1 = L_(1).real();
  const double lambda2 = L_(2).real();

  const double h0 = 1.0/sqrt(lambda0);
  const double h1 = 1.0/sqrt(lambda1);
  const double h2 = 1.0/sqrt(lambda2);

  const double distorsion0 = std::pow(lambda0/lambda1 - 1, 2);
  const double distorsion1 = std::pow(lambda0/lambda2 - 1, 2);
  const double distorsion2 = std::pow(lambda1/lambda2 - 1, 2);

  double distorsion = exp(-(distorsion0 + distorsion1 + distorsion2));

  return distorsion;
}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::nodesSpreading(std::vector<TInt>& nodesAdded, bool surfaceFlag)
{
  static int j = 0;
  static SimplexMesh test_mesh;
  std::vector<TInt> newNodes{};

  Variable<Eigen::Matrix3d>* metric  = nullptr;
  Variable<Eigen::Matrix3d>* metricNodes  = nullptr;
  gmds::Variable<int>* BND_SURFACE_COLOR    = nullptr;

  try{
    BND_SURFACE_COLOR    = m_nodesMesh.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
    metric = m_simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
    metricNodes = m_nodesMesh.getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  std::vector<nodeSamplingData> nodesSamplingData{};
  //std::clock_t start = std::clock();
  for(auto const node : nodesAdded)
  {
    Eigen::Matrix3d M = metricNodes->value(node);
    Eigen::Matrix3d FF = Eigen::Matrix3d::Zero();

    std::vector<math::Vector3d> frames{};
    const math::Point pt = SimplicesNode(&m_nodesMesh, node).getCoords();

    m_simplexMesh->getFrameAt(pt, frames);
    FF << frames[0].X() , frames[1].X() , frames[2].X(),
    frames[0].Y() , frames[1].Y() , frames[2].Y(),
    frames[0].Z() , frames[1].Z() , frames[2].Z();

    //std::cout << "FF : " << std::endl << FF << std::endl;
    //std::cout << "M : " << std::endl << M << std::endl;
    Eigen::Matrix3d t = metricInterpolationWithDistorsion(M, FF);
    std::vector<Eigen::Matrix3d> M_result{t, -t};

    //for(auto & f : frames)
    for(auto & t : M_result)
    {
      for(unsigned i = 0 ; i < 3 ; i++)
      {
        //Eigen::Vector3d d(f.X(), f.Y(), f.Z());
        Eigen::Vector3d d = t.col(i);
        //d = d.normalized();
        //M(0,0) = 1.0 /sqrt(M(0,0)); M(1,1) = 1.0 /sqrt(M(1,1)); M(2,2) = 1.0 /sqrt(M(2,2));
        //const Eigen::Vector3d vec = M * d;
        //const double m = vec.norm();
        //math::Point newCoord = pt + m*math::Point(d.x(), d.y(), d.z());
        math::Point newCoord = pt + math::Point(d.x(), d.y(), d.z());
        nodeSamplingData samplingData{};

        //if(m != 0.0)
        {
          const math::Point pt = SimplicesNode(&m_nodesMesh, node).getCoords();
          TInt newNodeId = -1;
          //compute the metric lentgh based on newNodeId and the metric at node
          findOptimimalPosition(node, newCoord);
          //compute the nearest Simplices of newCoord with octree
          std::vector<TSimplexID> simplices = m_oc.findSimplicesInOc(newCoord);
          bool inCell = false;
          for(auto const simplex : simplices)
          {
            if(simplex >= 0)
            {
              std::vector<double> UVWT{};
              //inCell = SimplicesCell(m_simplexMesh, simplex).isCellClose(newCoord, UVWT, 0.001);
              inCell = SimplicesCell(m_simplexMesh, simplex).isInCell(newCoord);
              if(inCell)
              {
                break;
              }
            }
          }

          //we project initiale node on surface
          //the good projection is determine by the minimum distance
          if(simplices.size() != 0)
          {
            double minDistance = std::numeric_limits<int>::max();
            if(!inCell)
            {
              double distance;
              math::Point goodProjection;
              for(auto const simplex : simplices)
              {
                if(simplex < 0)
                {
                  math::Point projectedPoint;
                  const SimplicesTriangle t(m_simplexMesh, simplex);
                  const std::vector<TInt> nodes = t.getNodes();
                  inCell = m_simplexMesh->pointInTriangle(newCoord,
                                       SimplicesNode(m_simplexMesh, nodes[0]).getCoords(),
                                       SimplicesNode(m_simplexMesh, nodes[1]).getCoords(),
                                       SimplicesNode(m_simplexMesh, nodes[2]).getCoords(),
                                       distance,projectedPoint);
                  if(inCell && distance < minDistance)
                  {
                    goodProjection = projectedPoint;
                    minDistance = distance;
                  }
                }
              }

              if(minDistance != std::numeric_limits<int>::max())
                newCoord = goodProjection;
            }

            bool onSurface = true;
            int surfaceLabel = 0;
            m_simplexMesh->onSurface(newCoord, surfaceLabel);

            if(surfaceLabel == 0)
              onSurface = false;
            if(surfaceFlag == onSurface)
            {
              if(minDistance != std::numeric_limits<int>::max() || inCell)
              {
                std::vector<TInt> neighboorNodes;
                nodeFiltering(newCoord, neighboorNodes);
                if(neighboorNodes.size() == 0)
                {
                  newNodeId = m_nodesMesh.addNode(newCoord);
                  if(surfaceFlag == true)
                    BND_SURFACE_COLOR->set(newNodeId, surfaceLabel);

                  m_nodesMesh.setAnalyticMetric(newNodeId, m_simplexMesh->getOctree());
                  newNodes.push_back(newNodeId);
                  std::unordered_set<TInt> seen{};
                  m_nodesMesh.getOctree()->addNode(newNodeId, seen);

                  TInt nnn = test_mesh.addNode(newCoord);
                  test_mesh.addTetraedre(nnn, nnn, nnn, nnn);
                  gmds::ISimplexMeshIOService ioServiceVOLUME(&test_mesh);
                  gmds::VTKWriter vtkWriterVOLUME(&ioServiceVOLUME);
                  vtkWriterVOLUME.setCellOptions(gmds::N|gmds::R);
                  vtkWriterVOLUME.setDataOptions(gmds::N|gmds::R);
                  vtkWriterVOLUME.write("testSurfaceNode1_" + std::to_string(j) + ".vtk");
                  ++j;
                }
              }
            }
          }
        }
      }
    }
  }

  //double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  //std::cout << "DURATION -> " << duration << std::endl;;

  nodesAdded.clear();
  nodesAdded = newNodes;
}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::processNodesStructure()
{
  const gmds::BitVector& nodesMesh = m_nodesMesh.getBitVectorNodes();
  for(unsigned int nodeId = 0 ; nodeId < nodesMesh.capacity() ; nodeId++)
  {
    if(nodesMesh[nodeId] != 0)
    {
      std::vector<TInt> neighboorNodes{};
      const math::Point pt = SimplicesNode(&m_nodesMesh, nodeId).getCoords();
      const double k = 1.2 * sqrt(2.0) ;
      //const double k = 0.7 * sqrt(2.0) * 0.5;
      nodeFiltering(pt, neighboorNodes, k);

      //comparing the distance and the angle divergence with the all the vector from frames
      //at current position
      std::vector<math::Vector3d> frames{};
      m_simplexMesh->getFrameAt(pt, frames);
      for(auto vec : frames)
      {
        TInt newNodeId = -1;
        vec.normalize();
        double angle_diff = 0.85;
        double distMin = std::numeric_limits<double>::max();
        for(auto const n : neighboorNodes)
        {
          math::Point pid = SimplicesNode(&m_nodesMesh, n).getCoords();
          math::Vector3d v = math::Vector3d({pid.X() - pt.X(), pid.Y() - pt.Y(), pid.Z() - pt.Z()});
          double dist = v.norm();
          v.normalize();
          double angle_dot = vec.dot(v);
          if(angle_dot >= angle_diff && dist <= distMin && n != nodeId)
          {
            distMin = dist;
            angle_diff = angle_dot;
            newNodeId = n;
          }
        }
        if(newNodeId != -1)
          m_nodeStructure[nodeId].push_back(newNodeId);
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::nodeFiltering(const math::Point& pt, std::vector<TInt> & neighboorNode, double k, bool flag)
{
  Variable<Eigen::Matrix3d>* metric  = nullptr;
  gmds::Variable<int>* BND_CURVE_COLOR_NODE = nullptr;

  try{
    metric = m_simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
    BND_CURVE_COLOR_NODE = m_nodesMesh.getVariable<int, SimplicesNode>("BND_CURVE_COLOR");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  const gmds::BitVector& nodesMesh = m_nodesMesh.getBitVectorNodes();
  //for(unsigned int nodeId = 0 ; nodeId < nodesMesh.capacity() ; nodeId++)

  std::vector<TInt> nodes = m_nodesMesh.getOctree()->findNodesNextTo(pt);

  for(auto const nodeId : nodes)
  {
    if(nodesMesh[nodeId] != 0)
    {
      const math::Point nodeCoord = SimplicesNode(&m_nodesMesh, nodeId).getCoords();
      Eigen::Vector3d v(pt.X() - nodeCoord.X(), pt.Y() - nodeCoord.Y(), pt.Z() - nodeCoord.Z());
      const double euclidianNorm = v.norm();
      //the metric being analytique we do not have to interpolate the metric at point
      const Eigen::Matrix3d m0 = m_nodesMesh.getAnalyticMetric(pt, m_simplexMesh->getOctree());
      const Eigen::Matrix3d m1 = m_nodesMesh.getAnalyticMetric(nodeCoord, m_simplexMesh->getOctree());
      //compute the current length based on the metric atached to the mesh
      const double metricLenght = 0.5 * sqrt(v.dot(m0*v)) + 0.5 * sqrt(v.dot(m1*v));

      if(metricLenght <= k)
      {
        neighboorNode.push_back(nodeId);
      }
    }
  }
}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::findOptimimalPosition(const TInt node, math::Point& initialCoord)
{
  double epsilon = 0.01;
  Variable<Eigen::Matrix3d>* metricNodes  = nullptr;
  try{
    metricNodes = m_nodesMesh.getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }


  //the node is the layer node so we can not move this once
  const math::Point pt = SimplicesNode(&m_nodesMesh, node).getCoords();
  Eigen::Vector3d dir(initialCoord.X() - pt.X(), initialCoord.Y() - pt.Y(), initialCoord.Z() - pt.Z());

  //the metric being analytique we do not have to interpolate the metric at point
  const Eigen::Matrix3d m0 = (*metricNodes)[node];
  const Eigen::Matrix3d m1 = m_simplexMesh->getAnalyticMetric(initialCoord, m_simplexMesh->getOctree());

  //compute the current length based on the metric atached to the mesh
  const double metricLenght = 0.5 * sqrt(dir.dot(m0*dir)) +  0.5 * sqrt(dir.dot(m1*dir));

  if(metricLenght < 1.0 - epsilon)
  {
    initialCoord = pt + 1.5 * math::Vector3d({dir.x(), dir.y(), dir.z()});
    findOptimimalPosition(node, initialCoord);
  }
  else if(metricLenght > 1.0 + epsilon)
  {
    initialCoord = 0.5 * (pt + initialCoord);
    findOptimimalPosition(node, initialCoord);
  }
}
/*----------------------------------------------------------------------------*/
std::vector<gmds::hybrid::MetricFFPointgeneration::DataEdges> MetricFFPointgeneration::subdivideEdge(const std::vector<TInt>& edge, const std::vector<double>& edgeU, const double sizeEdge) const
{
  //This function will subdivide the current edge if the metric is to much discontinue because of the curvature
  std::vector<DataEdges> dataEdges{};
  std::list<int> ts{};
  double epsilon = 1E-4;
  bool flag = true;
  /*for(int i = 1 ; i < edge.size() - 1 ; i++)
  {
    double tA = static_cast<double>(i) / static_cast<double>(edge.size()-1);
    double curvatureValueA = std::abs(curvature(tA, edge));
    if((curvatureValueA > epsilon) == flag)// the edge presente curvature at tA
    {
      ts.push_back(i-1);
      flag = !flag;
    }
  }*/

  for(int i = 1 ; i < edge.size() - 1 ; i++)
  {
    double tA = static_cast<double>(i) / static_cast<double>(edge.size()-1);
    double tB = static_cast<double>(i-1) / static_cast<double>(edge.size()-1);
    double curvatureValueA = std::abs(curvature(tA, edge));
    double curvatureValueB = std::abs(curvature(tB, edge));
    if(std::abs(curvatureValueA - curvatureValueB) > 1E-4)// the edge presente curvature at tA
    {
      ts.push_back(i);
    }
  }

  if(ts.size() % 2 != 0 && ts.size() != 1)
    throw gmds::GMDSException("ts.size() % 2 != 0 && ts.size() != 1");



  if(ts.front() != 0 || ts.size() == 0)
      ts.push_front(0);

  if(ts.back() != edge.size() - 1)
      ts.push_back(edge.size() - 1);

  for(std::list<int>::iterator it = ts.begin() ; it != ts.end() ; ++it)
  {
    if(it != ts.begin())
    {
      int val = *it;
      int p_val = *(std::next(it,-1));
      if(p_val + 1 == val)
      {
        ts.erase(it);
      }
    }
  }

  for(std::list<int>::iterator it = ts.begin() ; it != --(ts.end()) ; ++it)
  {
    int t0 = *it;
    int t1 = *(std::next(it,1));
    double ratio = (t1 - t0) / static_cast<double>(edge.size());
    double sizeSubEdge = ratio * sizeEdge;
    std::vector<int> SubEdge {};
    std::vector<double> SubEdgeU {};
    int cpt = 0;
    for(unsigned int j = t0 ; j <= t1; ++j)
    {
      SubEdgeU.push_back(static_cast<double>(cpt) / static_cast<double>(t1 - t0));
      SubEdge.push_back(edge[j]);
      ++cpt;
    }
    DataEdges dEdges;
    dEdges.subEdge = SubEdge;
    dEdges.subEdgeU = SubEdgeU;
    dEdges.sizeEdge = sizeSubEdge;
    dataEdges.push_back(dEdges);
  }

  return dataEdges;
}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::subdivideEdgeUsingMetric_Relaxation(std::vector<TInt>& nodesAdded, const std::vector<TInt>& edge_, const std::vector<double>& edgeU_, const double sizeEdge_, const unsigned int edgeId)
{
  double epsilon = 1E-3;
  Variable<Eigen::Matrix3d>* metric  = nullptr;
  gmds::Variable<int>* BND_CURVE_COLOR_NODE = nullptr;
  gmds::Variable<int>* BND_CURVE_VERTEX_NODE = nullptr;

  try{
    metric = m_simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
    BND_CURVE_COLOR_NODE = m_nodesMesh.getVariable<int, SimplicesNode>("BND_CURVE_COLOR");
    BND_CURVE_VERTEX_NODE = m_nodesMesh.getVariable<int, SimplicesNode>("BND_VERTEX_COLOR");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  //subdivide the edge if nedded ..
  std::vector<DataEdges> dataEdges = subdivideEdge(edge_, edgeU_, sizeEdge_);
  //compute the how many node should be in the edge  (depending on the Metric define on the mesh 's edges')
  for(auto const &d : dataEdges)
  {
    double den = 0.0;
    std::vector<double> edgeU = d.subEdgeU;
    std::vector<TInt> edge = d.subEdge;
    double sizeEdge = d.sizeEdge;

    for(unsigned int i = 0 ; i < edge.size() - 1 ; i++)
    {
      double tA = static_cast<double>(i) / static_cast<double>(edge.size()-1);
      double tB = static_cast<double>(i+1) / static_cast<double>(edge.size()-1);

      const TInt nodeA = edge[i];
      const TInt nodeB = edge[i + 1];


      const math::Point ptA = SimplicesNode(m_simplexMesh, nodeA).getCoords();
      const math::Point ptB = SimplicesNode(m_simplexMesh, nodeB).getCoords();

      Eigen::Vector3d dir = Eigen::Vector3d(ptB.X() - ptA.X(), ptB.Y() - ptA.Y(), ptB.Z() - ptA.Z());
      dir.normalize();

      //find the metric of the space in A and B
      Eigen::Matrix3d MA = metric->value(nodeA);
      Eigen::Matrix3d MB = metric->value(nodeB);

      //find the metric associate to the node (because of the curvature of the edge)
      double curvatureValueA = std::abs(curvature(tA, edge));
      double curvatureValueB = std::abs(curvature(tB, edge));
      double alpha = 0.1;

      if(curvatureValueA > epsilon)
      {
        //we built the isotrope metric on the canonical basis
        //on the adge at t with the curvatureValue (1/rayon de courbure)
        double curvatureRadiusA =  alpha * (1.0 / curvatureValueA);
        Eigen::Matrix3d curveMetricA;
        curveMetricA(0,0) = 1.0 / std::pow(curvatureRadiusA,2);
        curveMetricA(1,1) = 1.0 / std::pow(curvatureRadiusA,2);
        curveMetricA(2,2) = 1.0 / std::pow(curvatureRadiusA,2);

        MA = computeIntersectionMetric(curveMetricA, MA);
      }
      if(curvatureValueB > epsilon)
      {
        //we built the isotrope metric on the canonical basis
        //on the adge at t with the curvatureValue (1/rayon de courbure)
        double curvatureRadiusB =  alpha * (1.0 / curvatureValueB);
        Eigen::Matrix3d curveMetricB;
        curveMetricB(0,0) = 1.0 / std::pow(curvatureRadiusB,2);
        curveMetricB(1,1) = 1.0 / std::pow(curvatureRadiusB,2);
        curveMetricB(2,2) = 1.0 / std::pow(curvatureRadiusB,2);

        MB = computeIntersectionMetric(curveMetricB, MB);
      }
      Eigen::Matrix3d MA_modified;
      Eigen::Matrix3d MB_modified;

      MA_modified(0,0) = 1.0 /sqrt(MA(0,0)); MA_modified(1,1) = 1.0 /sqrt(MA(1,1)); MA_modified(2,2) = 1.0 /sqrt(MA(2,2));
      MB_modified(0,0) = 1.0 /sqrt(MB(0,0)); MB_modified(1,1) = 1.0 /sqrt(MB(1,1)); MB_modified(2,2) = 1.0 /sqrt(MB(2,2));



      const Eigen::Vector3d vecA = MA_modified * dir;
      const Eigen::Vector3d vecB = MB_modified * dir;

      const double mA = vecA.norm();
      const double mB = vecB.norm();

      const double sizeInterval = edgeU[i + 1] - edgeU[i];
      //const double sizeInterval = 1.0;

      //https://fr.wikipedia.org/wiki/Calcul_num%C3%A9rique_d%27une_int%C3%A9grale
      //basic discret square integral
      den += sizeInterval *  (0.5 * mA + 0.5 * mB);
    }

    if(den == 0.0)
    {
      throw gmds::GMDSException("den == 0.0");
    }

    //std::cout << "sizeEdge -> " << sizeEdge << std::endl;
    //std::cout << "den -> " << den << std::endl;
    unsigned int n = std::ceil(static_cast<double>(sizeEdge) / den);
    //std::cout << "n -> " << n << std::endl;

    //throw gmds::GMDSException("END");
    std::vector<double>res{};
    //std::cout << "init n -> " << n << std::endl;
    while(!metricSamplingEdge(n, res, edge, edgeU)){
      ++n;
      //std::cout << "n -> " << n << std::endl;
    }

    for(unsigned int i = 0 ; i < res.size()  ; i++)
    {
      const double u = res[i];
      Point pt = computeTheEdgeNodeCoordinate(u, edge, edgeU);
      std::vector<TSimplexID> tetraContenaingPt{};
      bool status = false;
      TInt newNodeId = m_nodesMesh.addNodeAndcheck(pt, tetraContenaingPt, status);
      if(!status)
      {
        (*BND_CURVE_COLOR_NODE)[newNodeId] = edgeId;
        m_nodesMesh.setAnalyticMetric(newNodeId, m_simplexMesh->getOctree());
        nodesAdded.push_back(newNodeId);
        std::unordered_set<TInt> seen{};
        m_nodesMesh.getOctree()->addNode(newNodeId, seen);
      }
    }
  }
  //throw gmds::GMDSException("END");

}
/*----------------------------------------------------------------------------*/
bool MetricFFPointgeneration::metricSamplingEdge(const unsigned int n, std::vector<double>& res, const std::vector<TInt>& edge, const std::vector<double>& edgeU) const
{
  if(edge.size() < 2)
  {
    throw gmds::GMDSException("edge.size < 2");
  }
  double epsilon = 1E-4;
  double minimumPercentageError = 0.05;
  //double epsilon = 0.0001;
  //first we subdivide the edge uniformly (whitout using any random generator for the position)
  std::vector<double> params_u{};
  const double e = 1.0 / static_cast<double>(n);
  for(unsigned int i = 0 ; i < n + 1 ; i++)
  {
    params_u.push_back(e * static_cast<double>(i));
  }
  //Compute the right sample node position in u coordinate with Loyd methodes in CVT -> Fast Methods for Computing Centroidal Voronoi Tessellations
  //compute V(zk) & the corresponding zk -> [0, 1]
  static unsigned int cpt = 0;
  double error_max_prev = 0.0;
  double error_max = error_max_prev;
  double error_moy = 0.0;
  do
  {
    error_max = 0.0;
    error_moy = 0.0;
    std::vector<std::vector<double>> V{};
    std::vector<double> p{0.0};
    for(unsigned int i = 1 ; i < params_u.size() - 1; i++)
    {
      const double uA = params_u[i - 1];
      const double uB = params_u[i    ];
      const double uC = params_u[i + 1];
      const double V0 = 0.5*(uA + uB);
      const double V1 = 0.5*(uC + uB);
      V.push_back(std::vector<double>{V0, V1});
    }

    //compute the new zk
    for(unsigned int i = 0; i < V.size() ; i++)
    {
      const double u_min = V[i][0];
      const double u_max = V[i][1];
      const double zk_1 = (u_min + u_max) * 0.5;

      const double tA = u_min;
      const double tB = u_max;

      const Point pt_min = computeTheEdgeNodeCoordinate(u_min, edge, edgeU);
      const Point pt_max = computeTheEdgeNodeCoordinate(u_max, edge, edgeU);

      Eigen::Vector3d dir = Eigen::Vector3d(pt_max.X() - pt_min.X(), pt_max.Y() - pt_min.Y(), pt_max.Z() - pt_min.Z());
      dir.normalize();

      Eigen::Matrix3d M_min  = m_simplexMesh->getAnalyticMetric(pt_min, m_simplexMesh->getOctree());
      Eigen::Matrix3d M_max  = m_simplexMesh->getAnalyticMetric(pt_max, m_simplexMesh->getOctree());

      //////////////////////////////////////////////////////////////////////////
      double curvatureValueA = std::abs(curvature(tA, edge));
      double curvatureValueB = std::abs(curvature(tB, edge));
      double alpha = 0.1;
      if(curvatureValueA > epsilon)
      {
        //we built the isotrope metric on the canonical basis
        //on the adge at t with the curvatureValue (1/rayon de courbure)
        double curvatureRadiusA =  alpha * (1.0 / curvatureValueA);
        Eigen::Matrix3d curveMetricA = Eigen::Matrix3d::Zero();
        curveMetricA(0,0) = 1.0 / std::pow(curvatureRadiusA,2);
        curveMetricA(1,1) = 1.0 / std::pow(curvatureRadiusA,2);
        curveMetricA(2,2) = 1.0 / std::pow(curvatureRadiusA,2);

        M_min = computeIntersectionMetric(curveMetricA, M_min);
      }
      if(curvatureValueB > epsilon)
      {
        //we built the isotrope metric on the canonical basis
        //on the adge at t with the curvatureValue (1/rayon de courbure)
        double curvatureRadiusB =  alpha * (1.0 / curvatureValueB);
        Eigen::Matrix3d curveMetricB = Eigen::Matrix3d::Zero();
        curveMetricB(0,0) = 1.0 / std::pow(curvatureRadiusB,2);
        curveMetricB(1,1) = 1.0 / std::pow(curvatureRadiusB,2);
        curveMetricB(2,2) = 1.0 / std::pow(curvatureRadiusB,2);

        M_max = computeIntersectionMetric(curveMetricB, M_max);
      }
      //////////////////////////////////////////////////////////////////////////

      M_min(0,0) = 1.0 /sqrt(M_min(0,0)); M_min(1,1) = 1.0 /sqrt(M_min(1,1)); M_min(2,2) = 1.0 /sqrt(M_min(2,2));
      M_max(0,0) = 1.0 /sqrt(M_max(0,0)); M_max(1,1) = 1.0 /sqrt(M_max(1,1)); M_max(2,2) = 1.0 /sqrt(M_max(2,2));


      const Eigen::Vector3d vecA = M_max * dir;
      const Eigen::Vector3d vecB = M_min * dir;

      const double r_min = vecA.norm();
      const double r_max = vecB.norm();

      const double num = u_min*r_min + u_max*r_max;
      const double den = r_min + r_max;

      if(den == 0.0)
      {
        throw gmds::GMDSException("den == 0");
      }

      const double zk = num / den;
      p.push_back(zk);

      error_max = std::max(error_max, std::abs(zk - zk_1));
      error_moy += std::abs(zk - zk_1) / params_u.size();
    }
    if(std::abs(error_max_prev - error_moy) / std::abs(error_moy) < minimumPercentageError)
      return false;

    //error_max_prev = error_max;
    error_max_prev = error_moy;
    p.push_back(1.0);
    params_u = p;

    /**************************************************************************/
    /*SimplexMesh sm;
    for(auto const u : params_u)
    {
      Point pt = computeTheEdgeNodeCoordinate(u, edge, edgeU);
      std::vector<TSimplexID> tetraContenaingPt{};
      bool status = false;
      TInt newNodeId = sm.addNode(pt);
      sm.addTetraedre(newNodeId, newNodeId, newNodeId, newNodeId);
    }
    gmds::ISimplexMeshIOService ioServiceM0(&sm);
    gmds::VTKWriter vtkWriterM0(&ioServiceM0);
    vtkWriterM0.setCellOptions(gmds::N|gmds::R);
    vtkWriterM0.write("Loyd_Relaxation_" + std::to_string(cpt)+ ".vtk");
    ++cpt;*/
    /**************************************************************************/

  }while(error_moy > epsilon);//while(error_max > epsilon);

      /**************************************************************************/
      /*static SimplexMesh sm;
      for(auto const u : params_u)
      {
        Point pt = computeTheEdgeNodeCoordinate(u, edge, edgeU);
        std::vector<TSimplexID> tetraContenaingPt{};
        bool status = false;
        TInt newNodeId = sm.addNode(pt);
        sm.addTetraedre(newNodeId, newNodeId, newNodeId, newNodeId);
      }
      static gmds::ISimplexMeshIOService ioServiceM0(&sm);
      gmds::VTKWriter vtkWriterM0(&ioServiceM0);
      vtkWriterM0.setCellOptions(gmds::N|gmds::R);
      vtkWriterM0.write("Loyd_Relaxation_" + std::to_string(cpt)+ ".vtk");
      ++cpt;*/
      /**************************************************************************/

  res = params_u;
  return true;
}
/*----------------------------------------------------------------------------*/
Eigen::VectorXd MetricFFPointgeneration::findCubicInterpolation(const double t, const std::vector<TInt>& edgeNodes, std::vector<int>& subEdgeIdx) const
{
  std::vector<TInt> nodesSubEdge = findBoundedNode(t, edgeNodes);
  if(nodesSubEdge.size() == 0)
    throw gmds::GMDSException("nodesSubEdge.size() == 0");

  subEdgeIdx = nodesSubEdge;
  int n_edges = edgeNodes.size();
  int n = nodesSubEdge.size();
  // Construction de la matrice des coefficients
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3*n, 3*n);
  for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        A(i, j) = std::pow(static_cast<double>(i) / static_cast<double>(n_edges-1), j);
        A(n+i, n+j) = std::pow(static_cast<double>(i) / static_cast<double>(n_edges-1), j);
        A(2*n+i, 2*n+j) = std::pow(static_cast<double>(i) / static_cast<double>(n_edges-1), j);
      }
  }

  //Calcul du vecteur des coefficients de la courbe polynomiale
  Eigen::VectorXd b = Eigen::VectorXd::Zero(3*n);
  for (int i = 0; i < n; i++) {
      b(i) = SimplicesNode(m_simplexMesh, edgeNodes[nodesSubEdge[i]]).getCoords()[0];
      b(n+i) = SimplicesNode(m_simplexMesh, edgeNodes[nodesSubEdge[i]]).getCoords()[1];
      b(2*n+i) = SimplicesNode(m_simplexMesh, edgeNodes[nodesSubEdge[i]]).getCoords()[2];
  }

  /*std::cout << A << std::endl;
  std::cout << std::endl;
  std::cout << b << std::endl;
  std::cout << std::endl;
  std::cout << "x : " << std::endl;
  std::cout << x << std::endl;
  std::cout <<  std::endl;*/
  Eigen::VectorXd x = A.inverse() * b;
  return x;
}
/*----------------------------------------------------------------------------*/
double MetricFFPointgeneration::curvature(const double t, const std::vector<TInt>& edgeNodes) const
{
  std::vector<int> subEdgesIdx{};
  Eigen::VectorXd coefficients = findCubicInterpolation(t, edgeNodes, subEdgesIdx);
  double new_t = t - subEdgesIdx[0] / static_cast<double>(edgeNodes.size() - 1);

  unsigned int n = coefficients.size() / 3;

  double courbureValue = 0.0;
  if(coefficients.size() == 12) //interpolation cubique
  {
    double c00 = coefficients[1]; double c01 = coefficients[2]; double c02 = coefficients[3];
    double c10 = coefficients[5]; double c11 = coefficients[6]; double c12 = coefficients[7];
    double c20 = coefficients[9]; double c21 = coefficients[10]; double c22 = coefficients[11];

    //dérivée premiere
    double dX = c00 +  2.0*c01*new_t + 3.0*c02*new_t*new_t;
    double dY = c10 +  2.0*c11*new_t + 3.0*c12*new_t*new_t;
    double dZ = c20 +  2.0*c21*new_t + 3.0*c22*new_t*new_t;

    //dérivée seconde
    double ddX = 2.0*c01 + 3.0*2.0*c02*new_t;
    double ddY = 2.0*c11 + 3.0*2.0*c12*new_t;
    double ddZ = 2.0*c21 + 3.0*2.0*c22*new_t;

    Eigen::Vector3d dXYZ(dX, dY, dZ);
    Eigen::Vector3d ddXYZ(ddX, ddY, ddZ);

    Eigen::Vector3d crossProduct = dXYZ.cross(ddXYZ);
    courbureValue = crossProduct.norm() / std::pow(dXYZ.norm(), 3);
  }
  else if (coefficients.size() == 9) //interpolation quadratique because we are on boder edge
  {
    double c00 = coefficients[1]; double c01 = coefficients[2];
    double c10 = coefficients[4]; double c11 = coefficients[5];
    double c20 = coefficients[7]; double c21 = coefficients[8];

    //dérivée premiere
    double dX = c00 + 2*c01*new_t;
    double dY = c10 + 2*c11*new_t;
    double dZ = c20 + 2*c21*new_t;

    //dérivée seconde
    double ddX = 2*c01;
    double ddY = 2*c11;
    double ddZ = 2*c21;

    Eigen::Vector3d dXYZ(dX, dY, dZ);
    Eigen::Vector3d ddXYZ(ddX, ddY, ddZ);

    Eigen::Vector3d crossProduct = dXYZ.cross(ddXYZ);
    courbureValue = crossProduct.norm() / std::pow(dXYZ.norm(), 3);
  }
  else
  {
    std::cout << "coefficients.size() -> " << coefficients.size() << std::endl;
    throw gmds::GMDSException("wrong coefficients.size()");
  }

  return courbureValue;
}
/*----------------------------------------------------------------------------*/
Eigen::Matrix3d MetricFFPointgeneration::computeIntersectionMetric(const Eigen::Matrix3d& m1, const Eigen::Matrix3d& m2) const
{
  Eigen::Matrix3cd m1_;
  m1_ << std::complex<double>(m1(0,0), 0), std::complex<double>(m1(0,1), 0), std::complex<double>(m1(0,2), 0),
  std::complex<double>(m1(1,0), 0), std::complex<double>(m1(1,1), 0), std::complex<double>(m1(1,2), 0),
  std::complex<double>(m1(2,0), 0), std::complex<double>(m1(2,1), 0), std::complex<double>(m1(2,2), 0);

  Eigen::Matrix3cd m2_;
  m2_ << std::complex<double>(m2(0,0), 0), std::complex<double>(m2(0,1), 0), std::complex<double>(m2(0,2), 0),
  std::complex<double>(m2(1,0), 0), std::complex<double>(m2(1,1), 0), std::complex<double>(m2(1,2), 0),
  std::complex<double>(m2(2,0), 0), std::complex<double>(m2(2,1), 0), std::complex<double>(m2(2,2), 0);


  Eigen::Matrix3d N = m1.inverse()*m2;

  Eigen::EigenSolver<Eigen::Matrix3d> es(N);
  Eigen::MatrixXcd V_ = es.eigenvectors();

  std::complex<double> lambda0_m1_ = (V_.col(0).transpose() * m1_ * V_.col(0));
  std::complex<double> lambda1_m1_ = (V_.col(1).transpose() * m1_ * V_.col(1));
  std::complex<double> lambda2_m1_ = (V_.col(2).transpose() * m1_ * V_.col(2));

  std::complex<double> lambda0_m2_ = (V_.col(0).transpose() * m2_ * V_.col(0));
  std::complex<double> lambda1_m2_ = (V_.col(1).transpose() * m2_ * V_.col(1));
  std::complex<double> lambda2_m2_ = (V_.col(2).transpose() * m2_ * V_.col(2));

  Eigen::Matrix3d d_;
  d_(0,0) = std::max(lambda0_m1_.real(), lambda0_m2_.real());
  d_(1,1) = std::max(lambda1_m1_.real(), lambda1_m2_.real());
  d_(2,2) = std::max(lambda2_m1_.real(), lambda2_m2_.real());

  //std::cout << "d -> " << std::endl << d_ << std::endl;
  return V_.real().inverse().transpose() * d_.real() * V_.real().inverse();
}
/*----------------------------------------------------------------------------*/
std::map<unsigned int, std::vector<TInt>> MetricFFPointgeneration::buildSortedEdges() const
{
  //sort the edge using the edge Id
  std::map<unsigned int, std::vector<TInt>> res{};
  const gmds::BitVector& meshNode = m_simplexMesh->getBitVectorNodes();
  Variable<int>* BND_VERTEX_COLOR    = nullptr;
  Variable<int>* BND_CURVE_COLOR     = nullptr;

  try{
    BND_VERTEX_COLOR  = m_simplexMesh->getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
    BND_CURVE_COLOR  = m_simplexMesh->getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  std::unordered_set<unsigned int> s{};
  std::vector<std::list<unsigned int>> linkedLists{};
  for(unsigned int node = 0 ; node < meshNode.capacity() ; node++)
  {
    if((*BND_CURVE_COLOR)[node] != 0 && s.find((*BND_CURVE_COLOR)[node]) == s.end())
    {
      s.insert((*BND_CURVE_COLOR)[node]);
      std::list<unsigned int> linkedList = {node};
      linkedLists.push_back(linkedList);
    }
  }



  std::vector<unsigned int> indices{};
  for(auto & linkedList : linkedLists)
  {
    s.clear();
    unsigned int currentIndice = (*BND_CURVE_COLOR)[linkedList.front()];
    indices.push_back(currentIndice);
    std::vector<TInt> neighborNodesFront = SimplicesNode(m_simplexMesh, linkedList.front()).neighborNodes();;
    std::vector<TInt> neighborNodesBack = SimplicesNode(m_simplexMesh, linkedList.back()).neighborNodes();;
    s.insert(linkedList.front()); // linkedList.front() = linkedList.back()
    bool nothingToAddFrontaly = false;
    bool nothingToAddBackely = false;
    while(!nothingToAddBackely || !nothingToAddFrontaly)
    {
      nothingToAddFrontaly = true;
      nothingToAddBackely = true;
      for(auto const neighborNode : neighborNodesFront)
      {
        if(((*BND_CURVE_COLOR)[neighborNode] == currentIndice || (*BND_VERTEX_COLOR)[neighborNode] != 0) && s.find(neighborNode) == s.end())
        {
          nothingToAddFrontaly = false;
          s.insert(neighborNode);
          linkedList.push_front(neighborNode);
          break;
        }
      }

      for(auto const neighborNode : neighborNodesBack)
      {
        if(((*BND_CURVE_COLOR)[neighborNode] == currentIndice || (*BND_VERTEX_COLOR)[neighborNode] != 0) && s.find(neighborNode) == s.end())
        {
          nothingToAddBackely = false;
          s.insert(neighborNode);
          linkedList.push_back(neighborNode);
          break;
        }
      }
      neighborNodesFront = SimplicesNode(m_simplexMesh, linkedList.front()).neighborNodes();
      neighborNodesBack = SimplicesNode(m_simplexMesh, linkedList.back()).neighborNodes();
    }
  }

  unsigned int i = 0;
  for(auto const & likedList : linkedLists)
  {
    std::vector<TInt> nodes{};
    for(auto const node : likedList)
    {
      nodes.push_back(node);
    }
    res[indices[i]] = nodes;
    i++;
  }

  return res;
}
/*----------------------------------------------------------------------------*/
std::vector<std::vector<double>> MetricFFPointgeneration::buildParamEdgeU(const std::map<unsigned int, std::vector<TInt>>& sortedEdge, std::vector<double> & length_edges) const
{
  std::vector<double> sizeEdges{};
  std::vector<std::vector<double>> sizeSubEdges{};

  for(auto const & edge : sortedEdge)
  {
    double AB_length = 0.0;
    std::vector<double> sizeSubEdge{0.0};
    for(unsigned int nodeIdx = 0 ; nodeIdx < edge.second.size() - 1; nodeIdx++)
    {
      TInt nodeA = edge.second[nodeIdx];
      TInt nodeB = edge.second[nodeIdx + 1];

      const Point ptA = SimplicesNode(m_simplexMesh, nodeA).getCoords();
      const Point ptB = SimplicesNode(m_simplexMesh, nodeB).getCoords();

      const Vector3d AB = ptB - ptA;
      AB_length += AB.norm();
      sizeSubEdge.push_back(AB_length);
    }
    length_edges.push_back(AB_length);
    sizeSubEdges.push_back(sizeSubEdge);
    sizeEdges.push_back(AB_length);
  }

  for(unsigned int subEdgeIdx = 0 ; subEdgeIdx < sizeSubEdges.size() ; subEdgeIdx++)
  {
    for(unsigned int i = 0 ; i < sizeSubEdges[subEdgeIdx].size() ; i++)
    {
      sizeSubEdges[subEdgeIdx][i] /= sizeEdges[subEdgeIdx];
    }
  }


  return sizeSubEdges;
}
