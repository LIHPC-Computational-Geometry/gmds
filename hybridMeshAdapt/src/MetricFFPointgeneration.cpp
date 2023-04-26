#include "gmds/hybridMeshAdapt/MetricFFPointgeneration.h"
#include "gmds/hybridMeshAdapt/SimplexMesh.h"
/*----------------------------------------------------------------------------*/
#include <gmds/ig/Mesh.h>
/*----------------------------------------------------------------------------*/
#include <gmds/io/VTKWriter.h>
/*----------------------------------------------------------------------------*/
#include <unordered_set>
#include <deque>
#include <queue>
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
MetricFFPointgeneration::MetricFFPointgeneration(SimplexMesh* simplexMesh, const std::string & name):m_simplexMesh(simplexMesh),m_oc(Octree(simplexMesh, 10)),m_layerNbr(0),m_minDistance(1.0 * (sqrt(2.0) * 0.5)), m_name(name)
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
void MetricFFPointgeneration::connectionWithNeighbor(const std::vector<TInt>& nodesAdded)
{
  gmds::Variable<int>* BND_VERTEX_COLOR = nullptr;
  gmds::Variable<int>* BND_CURVE_COLOR = nullptr;
  try{
    BND_CURVE_COLOR = m_nodesMesh.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  const gmds::BitVector& nodesBitVector = m_nodesMesh.getBitVectorNodes();
  for(auto const & n : nodesAdded)
  {
    if(nodesBitVector[n] != 0)
    {
      //look the predecessor who creat n
      TInt pred = m_nodeGeneratedBy[n].front();
      std::vector<TInt> connected_neigbors = m_nodeStructure[pred];
      for(auto const connected_neigbor : connected_neigbors)
      {
        const math::Point p = SimplicesNode(&m_nodesMesh, n).getCoords();
        auto comp = [=](TInt nodeA, TInt nodeB){
          const math::Point pA = SimplicesNode(&m_nodesMesh, nodeA).getCoords();
          const math::Point pB = SimplicesNode(&m_nodesMesh, nodeB).getCoords();
          const math::Vector3d vA = pA - p;
          const math::Vector3d vB = pB - p;
          double distA = vA.norm();
          double distB = vB.norm();
          return distA < distB;
        };
        std::set<int,decltype(comp)>  s = std::set<int,decltype(comp)>( comp );

        if(connected_neigbor != n)
        {
          std::vector<TInt> neigbors = m_nodeStructure[connected_neigbor];
          for(auto const neigbor : neigbors)
          {
            if(neigbor != pred)
            s.insert(neigbor);
          }
        }

        if(!s.empty())
        {
          if(m_nodeLayerNbr[n] == m_nodeLayerNbr[(*s.begin())])
          {
            m_nodeStructure[n].push_back((*s.begin()));
            m_nodeStructure[(*s.begin())].push_back(n);
          }
        }
      }
    }
  }
}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::addNodeToLayer(const TInt nodeId, const TInt fromNode, bool surfaceFlag)
{
  gmds::Variable<int>* BND_SURFACE_COLOR = nullptr;
  try{
    BND_SURFACE_COLOR = m_nodesMesh.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }


  m_layers[m_layerNbr].push_back(nodeId);
  m_nodeLayerNbr[nodeId] = m_layerNbr;

  if(m_layerNbr != 0) //recombine node for the grid structure
  {
    m_nodeStructure[nodeId].push_back(fromNode);
    m_nodeStructure[fromNode].push_back(nodeId);
  }
}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::incrementLayer()
{
  ++m_layerNbr;
}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::initializeGridWithEdge()
{
  for(auto const l : m_listEdge)
  {
    for(unsigned int i = 1 ; i < l.second.size() ; i++)
    {
      m_nodeStructure[l.second[i]].push_back(l.second[i-1]);
      m_nodeStructure[l.second[i-1]].push_back(l.second[i]);
    }
  }
}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::correctionVertexNode()
{
  double epsilon = 1E-6;
  gmds::Variable<int>* BND_CURVE_COLOR = nullptr;
  gmds::Variable<int>* BND_VERTEX_COLOR = nullptr;
  gmds::Variable<int>* BND_CURVE_VERTEX_COLOR = nullptr;


  try{
    BND_VERTEX_COLOR = m_simplexMesh->getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
    BND_CURVE_COLOR = m_nodesMesh.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
    BND_CURVE_VERTEX_COLOR = m_nodesMesh.getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  const gmds::BitVector& nodesBitVector = m_simplexMesh->getBitVectorNodes();
  for(auto const & l : m_layers[0])
  {
    for(unsigned int n = 0 ; n < nodesBitVector.capacity() ; n++)
    {
      if((*BND_VERTEX_COLOR)[n] != 0)
      {
        math::Point l_pt = SimplicesNode(&m_nodesMesh, l).getCoords();
        math::Point n_pt = SimplicesNode(m_simplexMesh, n).getCoords();
        math::Vector3d v = l_pt - n_pt;
        double dist = v.norm();
        if(dist <= epsilon)
        {
          BND_CURVE_VERTEX_COLOR->set(l, (*BND_VERTEX_COLOR)[n]);
          BND_CURVE_COLOR->set(l, 0);
        }
      }
    }
  }
}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::correctionNodeStructure()
{
  gmds::Variable<int>* BND_SURFACE_COLOR = nullptr;
  gmds::Variable<int>* BND_CURVE_COLOR = nullptr;
  gmds::Variable<int>* BND_VERTEX_COLOR = nullptr;
  gmds::Variable<Eigen::Matrix3d>* metric = nullptr;

  try{
    BND_SURFACE_COLOR = m_nodesMesh.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
    BND_CURVE_COLOR = m_nodesMesh.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
    BND_VERTEX_COLOR = m_nodesMesh.getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
    metric = m_nodesMesh.getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  //correction edge
  correctionVertexNode();

  const gmds::BitVector& nodeBitVector = m_nodesMesh.getBitVectorNodes();

  for(unsigned int n = 0 ; n < nodeBitVector.capacity() ; n++)
  {
    if(nodeBitVector[n] != 0)
    {
      if((*BND_VERTEX_COLOR)[n] != 0)
      {
      }
      else if((*BND_CURVE_COLOR)[n] != 0)
      {
        //check if n have 4 neigbors
        const std::map<unsigned int, std::pair<unsigned int, unsigned int>>&  e2s = m_simplexMesh->getEdgeTianglesIndices();
        if(m_nodeStructure[n].size() != 4)
        {
          //edge node n need 2 neighbors edge's node and 2 neighbors surface's node (from different surface)
          //edge is conneted for sur because of the initializeGridWithEdge() function
          //check the surfaces node
          std::vector<unsigned int> connectedSurface{e2s.at((*BND_CURVE_COLOR)[n]).first, e2s.at((*BND_CURVE_COLOR)[n]).second};
          connectedSurface.erase(std::remove_if(connectedSurface.begin(), connectedSurface.end(), [=](const unsigned int s){

            return ( std::find_if(m_nodeStructure[n].begin(), m_nodeStructure[n].end(), [=](const TInt neighbor){
              if((*BND_SURFACE_COLOR)[neighbor] != 0 && (*BND_SURFACE_COLOR)[neighbor] == s)
              {
                return true;
              }
              return false;
            }) != m_nodeStructure[n].end() );



          }), connectedSurface.end());
          for(auto const indiceSurface : connectedSurface)
          {
            TInt nodeToConnect = -1;
            double minDistance = std::numeric_limits<double>::max();
            for(auto const m : m_layers[1]) // first layer after edge is layer : 1
            {
              if(indiceSurface == (*BND_SURFACE_COLOR)[m] && std::find(m_nodeStructure[n].begin(), m_nodeStructure[n].end(), m) == m_nodeStructure[n].end())
              {
                math::Point p0 = SimplicesNode(&m_nodesMesh, m).getCoords();
                math::Point p1 = SimplicesNode(&m_nodesMesh, n).getCoords();
                math::Vector3d v = p0 - p1;
                double dist = v.norm();
                if(dist < minDistance)
                {
                  nodeToConnect = m;
                  minDistance = dist;
                }
              }
            }

            m_nodeStructure[n].push_back(nodeToConnect);
            if(std::find(m_nodeStructure[nodeToConnect].begin(), m_nodeStructure[nodeToConnect].end(), n) == m_nodeStructure[nodeToConnect].end())
            m_nodeStructure[nodeToConnect].push_back(n);
          }
        }
      }
      else if((*BND_SURFACE_COLOR)[n] != 0)
      {
        while(m_nodeStructure[n].size() < 4)
        {
          TInt nodeToConnect = -1;
          double minDistance = 2.1*m_minDistance;
          math::Point p0 = SimplicesNode(&m_nodesMesh, n).getCoords();
          Eigen::Matrix3d m0 = (*metric)[n];
          std::vector<TInt> neighborNodes = m_nodesMesh.getOctree()->findNodesNextTo(p0);
          for(auto const m : neighborNodes) // first layer after edge is layer : 1
          {
            if(n != m && (*BND_SURFACE_COLOR)[n] == (*BND_SURFACE_COLOR)[m] && std::find(m_nodeStructure[n].begin(), m_nodeStructure[n].end(), m) == m_nodeStructure[n].end()&&
              m_nodeLayerNbr[n] == m_nodeLayerNbr[m] /*m_nodeLayerNbr[n] >= m_nodeLayerNbr[m] */)
            {
              math::Point p1 = SimplicesNode(&m_nodesMesh, m).getCoords();
              Eigen::Matrix3d m1 = (*metric)[m];
              math::Vector3d v_ = p0 - p1;
              Eigen::Vector3d v(v_.X(), v_.Y(), v_.Z());
              double dist = 0.5 * sqrt(v.dot(m0*v)) + 0.5 * sqrt(v.dot(m1*v));
              if(dist < minDistance)
              {
                nodeToConnect = m;
                minDistance = dist;
              }
            }
          }
          if(nodeToConnect != -1)
          {
            m_nodeStructure[n].push_back(nodeToConnect);
            if(std::find(m_nodeStructure[nodeToConnect].begin(), m_nodeStructure[nodeToConnect].end(), n) == m_nodeStructure[nodeToConnect].end())
            m_nodeStructure[nodeToConnect].push_back(n);
          }
          else
            break;
        }
      }
      else //volume node
      {
        //check if n have 4 neigbors
        while(m_nodeStructure[n].size() < 5)
        {
          TInt nodeToConnect = -1;
          double minDistance = 2.1*m_minDistance;
          math::Point p0 = SimplicesNode(&m_nodesMesh, n).getCoords();
          Eigen::Matrix3d m0 = (*metric)[n];
          std::vector<TInt> neighborNodes = m_nodesMesh.getOctree()->findNodesNextTo(p0);
          for(auto const m : neighborNodes) // first layer after edge is layer : 1
          {
            if(n != m && std::find(m_nodeStructure[n].begin(), m_nodeStructure[n].end(), m) == m_nodeStructure[n].end() && m_nodeLayerNbr[n] == m_nodeLayerNbr[m])
            {
              math::Point p1 = SimplicesNode(&m_nodesMesh, m).getCoords();
              Eigen::Matrix3d m1 = (*metric)[m];
              math::Vector3d v_ = p0 - p1;
              Eigen::Vector3d v(v_.X(), v_.Y(), v_.Z());
              double dist = 0.5 * sqrt(v.dot(m0*v)) + 0.5 * sqrt(v.dot(m1*v));
              if(dist < minDistance)
              {
                nodeToConnect = m;
                minDistance = dist;
              }
            }
          }
          if(nodeToConnect != -1)
          {
            m_nodeStructure[n].push_back(nodeToConnect);
            if(std::find(m_nodeStructure[nodeToConnect].begin(), m_nodeStructure[nodeToConnect].end(), n) == m_nodeStructure[nodeToConnect].end())
            m_nodeStructure[nodeToConnect].push_back(n);
          }
          else
            break;
        }
      }
    }
  }
}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::correctUnwantedConnectionVOLUME()
{
  gmds::Variable<int>* BND_SURFACE_COLOR = nullptr;
  gmds::Variable<int>* BND_CURVE_COLOR = nullptr;
  gmds::Variable<int>* BND_VERTEX_COLOR = nullptr;
  gmds::Variable<Eigen::Matrix3d>* metric = nullptr;

  try{
    BND_SURFACE_COLOR = m_nodesMesh.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
    BND_CURVE_COLOR = m_nodesMesh.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
    BND_VERTEX_COLOR = m_nodesMesh.getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
    metric = m_nodesMesh.getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  const gmds::BitVector& nodeBitVector = m_nodesMesh.getBitVectorNodes();

  for(unsigned int n = 0 ; n < nodeBitVector.capacity() ; n++)
  {
    if(nodeBitVector[n] != 0)
    {
      if((*BND_SURFACE_COLOR)[n] == 0 && (*BND_CURVE_COLOR)[n] == 0 && (*BND_VERTEX_COLOR)[n] == 0)
      {
        //delete the next node not directly created by the nodesSpreading function
        //because of the correctionNodeStructure
        if(m_nodeStructure[n].size() > 5)
        {
          int currentLayer = m_nodeLayerNbr[n];
          std::vector<TInt> nodesConnectedToNextLayer{};
          std::vector<TInt> nodeToDelete{};

          m_nodeStructure[n].erase(std::remove_if(m_nodeStructure[n].begin(), m_nodeStructure[n].end(), [&](const TInt node){
            if(m_nodeLayerNbr[node] > currentLayer && std::find(m_nodeGeneratedBy[node].begin(), m_nodeGeneratedBy[node].end(), n) == m_nodeGeneratedBy[node].end() ){
              nodeToDelete.push_back(node);
              return true;
            }
            return false;
          }), m_nodeStructure[n].end());

          for(auto const next : nodeToDelete)
          {
            m_nodeStructure[next].erase(std::remove_if(m_nodeStructure[next].begin(), m_nodeStructure[next].end(), [&](const TInt beforeNode){
              if(beforeNode == n)
              return true;
              return false;
            }), m_nodeStructure[next].end());
          }
        }
      }
    }
  }
}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::correctUnwantedConnectionSURFACE()
{
  gmds::Variable<int>* BND_SURFACE_COLOR = nullptr;
  gmds::Variable<int>* BND_CURVE_COLOR = nullptr;
  gmds::Variable<int>* BND_VERTEX_COLOR = nullptr;
  gmds::Variable<Eigen::Matrix3d>* metric = nullptr;

  try{
    BND_SURFACE_COLOR = m_nodesMesh.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
    BND_CURVE_COLOR = m_nodesMesh.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
    BND_VERTEX_COLOR = m_nodesMesh.getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
    metric = m_nodesMesh.getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  const gmds::BitVector& nodeBitVector = m_nodesMesh.getBitVectorNodes();

  for(unsigned int n = 0 ; n < nodeBitVector.capacity() ; n++)
  {
    if(nodeBitVector[n] != 0)
    {
      if((*BND_SURFACE_COLOR)[n] == 0)
        continue;
      //delete the next node not directly created by the nodesSpreading function
      //because of the correctionNodeStructure
      if(m_nodeStructure[n].size() > 3)
      {
        int currentLayer = m_nodeLayerNbr[n];
        std::vector<TInt> nodesConnectedToNextLayer{};
        std::vector<TInt> nodeToDelete{};

        m_nodeStructure[n].erase(std::remove_if(m_nodeStructure[n].begin(), m_nodeStructure[n].end(), [&](const TInt node){
          if(m_nodeLayerNbr[node] > currentLayer && std::find(m_nodeGeneratedBy[node].begin(), m_nodeGeneratedBy[node].end(), n) == m_nodeGeneratedBy[node].end() ){
            nodeToDelete.push_back(node);
            return true;
          }
          return false;
        }), m_nodeStructure[n].end());

        for(auto const next : nodeToDelete)
        {
          m_nodeStructure[next].erase(std::remove_if(m_nodeStructure[next].begin(), m_nodeStructure[next].end(), [&](const TInt beforeNode){
            if(beforeNode == n)
            return true;
            return false;
          }), m_nodeStructure[next].end());
        }
      }
    }
  }
}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::sortBySurfaceNodeAdded(std::vector<TInt>& nodesAdded)
{
  gmds::Variable<int>* BND_SURFACE_COLOR = nullptr;

  try{
    BND_SURFACE_COLOR = m_nodesMesh.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  //first : we sort nodesAdded by surface color
  std::sort(nodesAdded.begin(), nodesAdded.end(), [=](const TInt nodeA, const TInt nodeB){
    return ((*BND_SURFACE_COLOR)[nodeA] != 0 && ((*BND_SURFACE_COLOR)[nodeA] < (*BND_SURFACE_COLOR)[nodeB]));
  });

  //second : we sort nodesAdded by neigbor node
  std::sort(nodesAdded.begin(), nodesAdded.end(), [=](const TInt nodeA, const TInt nodeB){
    if(m_nodeStructure.find(nodeA) != m_nodeStructure.end()){
      std::vector<TInt> v = m_nodeStructure[nodeA];
      return (std::find(v.begin(), v.end(), nodeB) != v.end());
    }
    return false;
  });
}
/*----------------------------------------------------------------------------*/
bool MetricFFPointgeneration::belongToEdge(const math::Point & nodeCoord)
{
  double epsilon = 1E-4;
  for(auto const & l : m_listEdge)
  {
    for(unsigned int i = 0 ; i < l.second.size() -1 ; i++)
    {
      const math::Point p0 = SimplicesNode(&m_nodesMesh, l.second[i]).getCoords();
      const math::Point p1 = SimplicesNode(&m_nodesMesh, l.second[i+1]).getCoords();


      const math::Vector3d v10 = p1-p0;
      const math::Vector3d v = nodeCoord-p0;

      if(std::abs((v10.cross(v)).norm()) > epsilon)
      continue;

      if(v10.dot(v) < 0.0)
      continue;

      double normEDGE = v10.norm();
      double d = v10.dot(v);
      if(d > normEDGE)
      continue;

      return true;
    }
  }
  return false;
}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::execute()
{
  std::clock_t start = std::clock();
  gmds::Variable<int>* BND_VERTEX_COLOR = nullptr;
  gmds::Variable<int>* BND_CURVE_COLOR = nullptr;
  gmds::Variable<int>* BND_SURFACE_COLOR = nullptr;

  std::vector<double> simplexMesh_Borders = m_oc.getBorderOctree();
  Octree* simplexNodes_Octree = new Octree(&m_nodesMesh, 3,
    simplexMesh_Borders[0],simplexMesh_Borders[1],
    simplexMesh_Borders[2],simplexMesh_Borders[3],
    simplexMesh_Borders[4],simplexMesh_Borders[5]);
    m_nodesMesh.setOctree(simplexNodes_Octree);


    try{
      BND_VERTEX_COLOR = m_nodesMesh.getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
      BND_CURVE_COLOR = m_nodesMesh.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
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
      const unsigned int edgeId = sortedEdge.first;
      std::vector<TInt> edge = sortedEdge.second;
      std::vector<double> edgeU = edgesU[i];
      double edge_length = edges_length[i];
      subdivideEdgeUsingMetric_Relaxation(nodeAdded, edge, edgeU, edge_length, edgeId);
      i++;
    }
    initializeGridWithEdge();
    incrementLayer();

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
    vtkWriterEGE.setCellOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterEGE.setDataOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterEGE.write("metricFF_EDGE_" + m_name + ".vtk");

    std::queue<TInt> q{};
    while(nodeAdded.size() != 0)
    {
      //fill the stack (change after the nodesSpreading function)
      for(unsigned int i = 0 ; i < nodeAdded.size() ; i++)
      q.push(nodeAdded[i]);
      std::cout << "    nodeAdded.size() -> " << nodeAdded.size() << std::endl;
      nodesSpreading(nodeAdded, true);
      incrementLayer();

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
      gmds::ISimplexMeshIOService ioServiceGRIDTEST(&m_nodesMesh);
      gmds::VTKWriter vtkWriterGRIDTEST(&ioServiceGRIDTEST);
      vtkWriterGRIDTEST.setCellOptions(gmds::N|gmds::R|gmds::F);
      vtkWriterGRIDTEST.setDataOptions(gmds::N|gmds::R|gmds::F);
      vtkWriterGRIDTEST.write("metricFF_Grid_LAYER_SORTING_COLOR_NEIGHBOR_" + m_name + "_" + std::to_string(m_layerNbr) + ".vtk");
    }

    correctionNodeStructure();
    std::cout << "correctionNodeStructure() for SURFACE END" << std::endl;
    correctUnwantedConnectionSURFACE();
    std::cout << "correctUnwantedConnection() for SURFACE END " << std::endl;


    for(auto const & d : m_nodeStructure)
    {
      TInt node = d.first;
      if(node != -1)
      {
        for(auto const n : d.second)
        {
          if(n != -1)
          {
            m_nodesMesh.addTriangle(d.first, n, n);
          }
        }
      }
    }

    gmds::ISimplexMeshIOService ioServiceGRIDTEST(&m_nodesMesh);
    gmds::VTKWriter vtkWriterGRIDTEST(&ioServiceGRIDTEST);
    vtkWriterGRIDTEST.setCellOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterGRIDTEST.setDataOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterGRIDTEST.write("metricFF_Grid_LAYER_SORTING_COLOR_NEIGHBOR_" + m_name + "_" + std::to_string(m_layerNbr) + ".vtk");

    gmds::ISimplexMeshIOService ioServiceSURFACE(&m_nodesMesh);
    gmds::VTKWriter vtkWriterSURFACE(&ioServiceSURFACE);
    vtkWriterSURFACE.setCellOptions(gmds::N|gmds::R);
    vtkWriterSURFACE.setDataOptions(gmds::N|gmds::R);
    vtkWriterSURFACE.write("metricFF_SURFACE_" + m_name + ".vtk");
    std::cout << "SURFACE NODES CREATED " << std::endl;

    const gmds::BitVector& m_nodes = m_nodesMesh.getBitVectorNodes();
    for(unsigned int n = 0 ; n < m_nodes.capacity() ; n++)
    {
      if(m_nodes[n] != 0)
      {
        if((*BND_SURFACE_COLOR)[n] == 2/*!= 0*/)
        {
          nodeAdded.push_back(n);
        }
      }
    }
    sortBySurfaceNodeAdded(nodeAdded);

    //fill the stack (change after the nodesSpreading function)
    nodeAdded.clear();
    unsigned int sizeQ = q.size();
    for(unsigned int i = 0 ; i < sizeQ ; i++)
    {
      nodeAdded.push_back(q.front());
      q.pop();
    }

    while(nodeAdded.size() != 0)
    {
      std::cout << "    nodeAdded.size() -> " << nodeAdded.size() << std::endl;
      nodesSpreading(nodeAdded);
      incrementLayer();


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
      gmds::ISimplexMeshIOService ioServiceGRIDTEST(&m_nodesMesh);
      gmds::VTKWriter vtkWriterGRIDTEST(&ioServiceGRIDTEST);
      vtkWriterGRIDTEST.setCellOptions(gmds::N|gmds::R|gmds::F);
      vtkWriterGRIDTEST.setDataOptions(gmds::N|gmds::R|gmds::F);
      vtkWriterGRIDTEST.write("metricFF_Grid_LAYER_SORTING_COLOR_NEIGHBOR_" + m_name + "_" + std::to_string(m_layerNbr) + ".vtk");
    }

    correctionNodeStructure();
    std::cout << "correctionNodeStructure() for VOLUME END" << std::endl;
    correctUnwantedConnectionVOLUME();
    std::cout << "correctUnwantedConnection() for VOLUME END " << std::endl;


    std::cout << "END COMPUTING " << std::endl;
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


    std::set<std::vector<TInt>> hexs{};
    computeHexa(hexs);

    std::cout << "hex size -> " << hexs.size() << std::endl;

    Mesh m(MeshModel(DIM3 | R | F | E | N |
      R2N | F2N | E2N | R2F | F2R |
      F2E | E2F | R2E | N2R | N2F | N2E));

      Variable<int>* BND_VERTEX_COLOR_END = m.newVariable<int, GMDS_NODE>("BND_VERTEX_COLOR"  );
      Variable<int>* BND_CURVE_COLOR_END = m.newVariable<int, GMDS_NODE>("BND_CURVE_COLOR"  );
      Variable<int>* BND_SURFACE_COLOR_END = m.newVariable<int, GMDS_NODE>("BND_SURFACE_COLOR");

      BND_VERTEX_COLOR_END->setValuesTo(-1);
      BND_CURVE_COLOR_END->setValuesTo(-1);
      BND_SURFACE_COLOR_END->setValuesTo(-1);


      const gmds::BitVector& m_nodesEND = m_nodesMesh.getBitVectorNodes();
      unsigned int cpt = 0;
      std::unordered_map<TInt, Node> map{};
      for(unsigned int n = 0 ; n < m_nodesEND.capacity() ; n++)
      {
        if(m_nodesEND[n] != 0)
        {
          const Node nodeHEX =  m.newNode(SimplicesNode(&m_nodesMesh, n).getCoords());
          m.newTet(nodeHEX,nodeHEX,nodeHEX,nodeHEX);
          map[n] = nodeHEX;
          if((*BND_VERTEX_COLOR)[n] != 0)
          {
            BND_VERTEX_COLOR_END->set(cpt, (*BND_VERTEX_COLOR)[n]);
          }
          else if((*BND_CURVE_COLOR)[n] != 0)
          {
            BND_CURVE_COLOR_END->set(cpt, (*BND_CURVE_COLOR)[n]);
          }
          else if((*BND_SURFACE_COLOR)[n] != 0)
          {
            BND_SURFACE_COLOR_END->set(cpt, (*BND_SURFACE_COLOR)[n]);
          }
          ++cpt;
        }
      }


      for(auto const & hex : hexs)
      {
        m.newHex(map[hex[0]], map[hex[1]], map[hex[2]], map[hex[3]], map[hex[4]], map[hex[5]], map[hex[6]], map[hex[7]]);
      }


      gmds::IGMeshIOService ioService(&m);
      gmds::VTKWriter vtkWriter(&ioService);
      vtkWriter.setCellOptions(gmds::N|gmds::R);
      vtkWriter.setDataOptions(gmds::N|gmds::R);
      vtkWriter.write(m_name+".vtk");
      std::cout << "HEX CREATED " << std::endl;


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
            //m_simplexMesh->onSurface(SimplicesNode(&m_nodesMesh,n).getCoords(), surfaceLabel);
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
Point MetricFFPointgeneration::computeTheEdgeNodeCoordinate(const double u, const std::vector<TInt>& edge, const std::vector<double>& edgeU, TInt& nodeA_, TInt& nodeB_) const
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
  TInt nodeA;
  TInt nodeB;
  for(unsigned int i = 0 ; i < edgeU.size() - 1 ; i++)
  {
    const double uA = edgeU[i];
    const double uB = edgeU[i + 1];
    if(uA <= u && uB >= u)
    {
      //interpolation (linear) of the middle position using the interval ua, ub and u = 0.5
      if(uB != uA)
      {
        nodeA = edge[i];
        nodeB = edge[i + 1];

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

  nodeA_ = nodeA;
  nodeB_ = nodeB;
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
  distorsion = 1.0;//1E-4;

  /*Eigen::EigenSolver<Eigen::Matrix3d> es_metric(metric);
  Eigen::MatrixXd P_metric = es_metric.eigenvectors().real();
  Eigen::VectorXd D_metric = es_metric.eigenvalues().real();

  Eigen::Vector3d t0 = distorsion*frameField.col(0) + (1.0-distorsion) * 1.0/sqrt(D_metric(0))*P_metric.col(0);
  Eigen::Vector3d t1 = distorsion*frameField.col(1) + (1.0-distorsion) * 1.0/sqrt(D_metric(1))*P_metric.col(1);
  Eigen::Vector3d t2 = distorsion*frameField.col(2) + (1.0-distorsion) * 1.0/sqrt(D_metric(2))*P_metric.col(2);*/

  Eigen::Vector3d t0 = distorsion*frameField.col(0);
  Eigen::Vector3d t1 = distorsion*frameField.col(1);
  Eigen::Vector3d t2 = distorsion*frameField.col(2);


  Eigen::Matrix3d t = Eigen::Matrix3d::Zero();
  t.col(0) << t0 / (/*sqrt*/(t0.transpose() * metric * t0));
  t.col(1) << t1 / (/*sqrt*/(t1.transpose() * metric * t1));
  t.col(2) << t2 / (/*sqrt*/(t2.transpose() * metric * t2));

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  /*Eigen::Vector3d col0 = frameField.col(0).real();
  Eigen::Vector3d col1 = frameField.col(1).real();
  Eigen::Vector3d col2 = frameField.col(2).real();

  double a = sqrt(col0.transpose() * metric * col0);
  double b = sqrt(col1.transpose() * metric * col1);
  double c = sqrt(col2.transpose() * metric * col2);

  Eigen::Matrix3d scalarReduction = Eigen::Matrix3d::Identity();
  scalarReduction(0,0) = 1.0/a;
  scalarReduction(1,1) = 1.0/b;
  scalarReduction(2,2) = 1.0/c;

  Eigen::Matrix3cd t = (1-distorsion)*metric + distorsion*scalarReduction*frameField;*/
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  return t;

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
  std::clock_t c_start = std::clock();
  double durationTOT;
  double durationRec;
  double durationInSim;
  double durationFiltering;
  double durationAddNode;
  double durationSurface;
  double durationExistingNode;
  double durationGetFrame;
  double durationFOR;
  double duration;

  const double epsilon = 5E-2;
  std::vector<TInt> newNodes{};
  Variable<Eigen::Matrix3d>* metric  = nullptr;
  Variable<Eigen::Matrix3d>* metricNodes  = nullptr;
  gmds::Variable<int>* BND_SURFACE_COLOR    = nullptr;
  gmds::Variable<int>* BND_CURVE_COLOR    = nullptr;

  try{
    BND_SURFACE_COLOR    = m_nodesMesh.getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
    BND_CURVE_COLOR    = m_nodesMesh.getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
    metric = m_simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
    metricNodes = m_nodesMesh.getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  SimplexMesh testMesh;
  SimplexMesh testMesh0;
  for(auto const node : nodesAdded)
  {

    Eigen::Matrix3d M = metricNodes->value(node);
    Eigen::Matrix3d FF = Eigen::Matrix3d::Zero();

    std::vector<math::Vector3d> frames{};
    const math::Point pt = SimplicesNode(&m_nodesMesh, node).getCoords();

    std::clock_t c_start7 = std::clock();
    m_simplexMesh->getFrameAt(pt, frames);
    durationGetFrame +=  std::clock() - c_start7;
    FF << frames[0].X() , frames[1].X() , frames[2].X(),
    frames[0].Y() , frames[1].Y() , frames[2].Y(),
    frames[0].Z() , frames[1].Z() , frames[2].Z();

    Eigen::Matrix3d t = metricInterpolationWithDistorsion(M, FF);
    std::vector<Eigen::Matrix3d> M_result{t, -t};

    for(auto & t : M_result)
    {
      std::clock_t c_start8 = std::clock();
      for(unsigned i = 0 ; i < 3 ; i++)
      {
        Eigen::Vector3d d = t.col(i);
        math::Point newCoord = pt + math::Point(d.x(), d.y(), d.z());
        nodeSamplingData samplingData{};
        const math::Point pt = SimplicesNode(&m_nodesMesh, node).getCoords();
        TInt newNodeId = -1;
        std::clock_t c_start1 = std::clock();
        if(!findOptimimalPosition(node, newCoord, surfaceFlag))
          continue;
        durationRec += std::clock() - c_start1;
        //check if the node is on the mesh during the volume propagation
        std::clock_t c_start2 = std::clock();
        bool flag = false;
        TSimplexID simplex = std::numeric_limits<TSimplexID>::min();
        if(!surfaceFlag)
        {
          std::vector<TSimplexID> simplices = m_oc.findSimplicesInOc(newCoord);
          if(simplices.size() == 0) //outside the octree
          flag = true;
          for(auto const & s : simplices)
          {
            if(s >= 0)
            {
              if(SimplicesCell(m_simplexMesh, s).isInCell(newCoord))
              {
                flag = true;
                break;
              }
            }
          }
          if(!flag)
          continue;
        }
        durationInSim += std::clock() - c_start2;

        //if(belongToEdge(newCoord))
          //continue;


        bool onSurface = true;
        int surfaceLabel = 0;

        std::clock_t c_start5 = std::clock();
        if(surfaceFlag)
          m_simplexMesh->onSurface(newCoord, surfaceLabel, simplex);
        durationSurface += std::clock() - c_start5;

        if(surfaceLabel == 0)
          onSurface = false;
        if(surfaceFlag == onSurface)
        {
          std::clock_t c_start3 = std::clock();
          std::vector<TInt> neighboorNodes;
          if(!nodeFiltering(newCoord, node, simplex, neighboorNodes))
            continue;
          durationFiltering += std::clock() - c_start3;

          if(neighboorNodes.size() == 0)
          {
            std::clock_t c_start4 = std::clock();
            newNodeId = m_nodesMesh.addNode(newCoord);
            nodeBelongingTO[newNodeId] = simplex;
            m_nodeGeneratedBy[newNodeId].push_back(node);
            if(surfaceFlag == true)
            BND_SURFACE_COLOR->set(newNodeId, surfaceLabel);

            newNodes.push_back(newNodeId);
            addNodeToLayer(newNodeId, node, surfaceFlag);
            m_nodesMesh.setAnalyticMetricFromMesh(newNodeId, m_simplexMesh, nodeBelongingTO);
            std::unordered_set<TInt> seen{};
            m_nodesMesh.getOctree()->addNode(newNodeId, seen);
            durationAddNode += std::clock() - c_start4;
          }
          else
          {
            std::clock_t c_start6 = std::clock();
            //if(!surfaceFlag)
            {
              int neigborNode = -1;
              double distanceMin = std::numeric_limits<double>::max();
              for(auto const n : neighboorNodes)
              {
                const math::Point pt_n = SimplicesNode(&m_nodesMesh, n).getCoords();
                math::Vector3d v = pt_n - newCoord;
                double dist = v.norm();
                if(dist < distanceMin)
                {
                  distanceMin = dist;
                  neigborNode = n;
                }
              }

              if(distanceMin <= epsilon &&  neigborNode!= -1 &&
                std::find(m_nodeStructure[node].begin(), m_nodeStructure[node].end(), neigborNode) == m_nodeStructure[node].end()
                && node != neigborNode && m_nodeLayerNbr[node] == m_nodeLayerNbr[neigborNode])
                {

                  m_nodeStructure[node].push_back(neigborNode);
                  m_nodeGeneratedBy[neigborNode].push_back(node);
                  if(std::find(m_nodeStructure[neigborNode].begin(), m_nodeStructure[neigborNode].end(), node) == m_nodeStructure[neigborNode].end())
                  m_nodeStructure[neigborNode].push_back(node);
                }
              }
              durationExistingNode += std::clock() - c_start6;
            }
          }
        }
        durationFOR += std::clock() - c_start8;
      }
    }

    durationTOT = ( std::clock() - c_start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"      durationTOT : "<< durationTOT << std::endl;;
    std::cout<<"        durationRec : "<< durationRec / (double) CLOCKS_PER_SEC << std::endl;;
    std::cout<<"        durationInSim : "<< durationInSim / (double) CLOCKS_PER_SEC << std::endl;;
    std::cout<<"        durationFiltering : "<< durationFiltering / (double) CLOCKS_PER_SEC << std::endl;;
    std::cout<<"        durationAddNode : "<< durationAddNode / (double) CLOCKS_PER_SEC << std::endl;
    std::cout<<"        durationExistingNode : "<< durationExistingNode / (double) CLOCKS_PER_SEC << std::endl;;
    std::cout<<"        durationGetFrame : "<< durationGetFrame / (double) CLOCKS_PER_SEC << std::endl;;
    std::cout<<"        durationFOR : "<< durationFOR / (double) CLOCKS_PER_SEC << std::endl;;
    std::cout<<"        durationSurface : "<< durationSurface / (double) CLOCKS_PER_SEC << std::endl;;
    std::cout << std::endl;
    std::cout << std::endl;

    nodesAdded.clear();
    std::copy(newNodes.begin(), newNodes.end(), std::back_inserter(nodesAdded));
  }
  /*----------------------------------------------------------------------------*/
  bool MetricFFPointgeneration::nodeFiltering(const math::Point& pt, const TInt fromNode, const TSimplexID simplex, std::vector<TInt> & neighboorNode)
  {
    Variable<Eigen::Matrix3d>* metric  = nullptr;
    gmds::Variable<int>* BND_CURVE_COLOR_NODE = nullptr;

    try{
      //metric = m_simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
      metric = m_nodesMesh.getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
      BND_CURVE_COLOR_NODE = m_nodesMesh.getVariable<int, SimplicesNode>("BND_CURVE_COLOR");
    }catch (gmds::GMDSException e)
    {
      throw gmds::GMDSException(e);
    }

    const gmds::BitVector& nodesMesh = m_nodesMesh.getBitVectorNodes();
    //for(unsigned int nodeId = 0 ; nodeId < nodesMesh.capacity() ; nodeId++)

    std::vector<TInt> nodes = m_nodesMesh.getOctree()->findNodesNextTo(pt);
    Eigen::Matrix3d m0;
    if(simplex != std::numeric_limits<TSimplexID>::min()){
      m0 = m_nodesMesh.getAnalyticMetricFromSimplex(pt, m_simplexMesh, simplex);
    }
    else{
      bool status = true;
      m0 = m_nodesMesh.getAnalyticMetric(pt, m_simplexMesh, status);
      if(!status)
      return status;
    }


    for(auto const nodeId : nodes)
    {
      if(nodesMesh[nodeId] != 0/* && nodeId != fromNode*/)
      {
        const math::Point nodeCoord = SimplicesNode(&m_nodesMesh, nodeId).getCoords();
        Eigen::Vector3d v(pt.X() - nodeCoord.X(), pt.Y() - nodeCoord.Y(), pt.Z() - nodeCoord.Z());
        const double euclidianNorm = v.norm();
        //the metric being analytique we do not have to interpolate the metric at point
        //const Eigen::Matrix3d m1 = m_nodesMesh.getAnalyticMetricFromSimplex(nodeCoord, m_simplexMesh, );
        const Eigen::Matrix3d m1 = (*metric)[nodeId];
        //compute the current length based on the metric atached to the mesh
        const double metricLenght = 0.5 * sqrt(v.dot(m0*v)) + 0.5 * sqrt(v.dot(m1*v));

        if(metricLenght <= m_minDistance)
        {
          neighboorNode.push_back(nodeId);
        }
      }
    }
    return true;
  }
  /*----------------------------------------------------------------------------*/
  bool MetricFFPointgeneration::findOptimimalPosition(const TInt node, math::Point& initialCoord, bool surfaceFlag, int cpt, double epsilon)
  {
    if(cpt == 0){
      epsilon *= 2.0;
      cpt = 10;
    }


    double DistanceMin = 1E-6;
    Variable<Eigen::Matrix3d>* metricNodes  = nullptr;
    try{
      metricNodes = m_nodesMesh.getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
    }catch (gmds::GMDSException e)
    {
      throw gmds::GMDSException(e);
    }

    //////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////
    math::Point newCoord = initialCoord;
    TSimplexID triangleProj;
    if(surfaceFlag)
    {
      double minDistance = std::numeric_limits<int>::max();
      double distance;
      math::Point goodProjection;
      //const gmds::BitVector& triangleVec = m_simplexMesh->getBitVectorTri();
      //for(unsigned int triangle = 1 ; triangle < triangleVec.capacity() ; triangle++)
      std::vector<TInt> triangles = m_simplexMesh->getOctree()->findTriangleInOc(newCoord);
      for(auto const triangle : triangles)
      {
        //if(triangleVec[triangle] != 0)
        if(triangle < 0)
        {
          math::Point projectedPoint;
          const SimplicesTriangle t(m_simplexMesh, triangle);
          const std::vector<TInt> nodes = t.getNodes();
          bool inCell = m_simplexMesh->pointInTriangle(newCoord,
            SimplicesNode(m_simplexMesh, nodes[0]).getCoords(),
            SimplicesNode(m_simplexMesh, nodes[1]).getCoords(),
            SimplicesNode(m_simplexMesh, nodes[2]).getCoords(),
            distance,projectedPoint);

            if(inCell && std::abs(distance) < minDistance)
            {
              goodProjection = projectedPoint;
              minDistance = std::abs(distance);
              triangleProj = triangle;
            }
          }
        }

        if(minDistance != std::numeric_limits<int>::max()){
          newCoord = goodProjection;
        }
        else{
          return false;
        }
      }
      //////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////

      //the node is the layer node so we can not move this once
      const math::Point pt = SimplicesNode(&m_nodesMesh, node).getCoords();
      Eigen::Vector3d dir(newCoord.X() - pt.X(), newCoord.Y() - pt.Y(), newCoord.Z() - pt.Z());

      if(dir.norm() <= DistanceMin){
        initialCoord = newCoord;
        return false;
      }

      //the metric being analytique we do not have to interpolate the metric at point
      const Eigen::Matrix3d m0 = (*metricNodes)[node];
      Eigen::Matrix3d m1;

      if(surfaceFlag){
        m1 = m_simplexMesh->getAnalyticMetricFromSimplex(initialCoord, m_simplexMesh, SimplicesTriangle(m_simplexMesh, triangleProj).neighborSimplex()[0]);
      }
      else{
        bool status = true;
        m1 = m_simplexMesh->getAnalyticMetric(initialCoord, m_simplexMesh, status);

        if(!status) // the node is probably outside the mesh so there is no metric
        return status;
      }
      /*bool status = true;
      const Eigen::Matrix3d m1 = m_simplexMesh->getAnalyticMetric(initialCoord, m_simplexMesh->getOctree(), status);

      if(!status) // the node is probably outside the mesh so there is no metric
      return status;*/

      //compute the current length based on the metric atached to the mesh
      const double metricLenght = 0.5 * sqrt(dir.dot(m0*dir)) +  0.5 * sqrt(dir.dot(m1*dir));

      if(metricLenght < (1.0 - epsilon))
      {
        initialCoord = pt + 1.5 * math::Vector3d({dir.x(), dir.y(), dir.z()});
        return findOptimimalPosition(node, initialCoord, surfaceFlag, --cpt, epsilon);
      }
      else if(metricLenght > (1.0 + epsilon))
      {
        initialCoord = 0.5 * (pt + initialCoord);
        return findOptimimalPosition(node, initialCoord, surfaceFlag, --cpt, epsilon);
      }
      else
      {
        initialCoord = newCoord;
        return true;
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
void MetricFFPointgeneration::subdivideEdgeUsingMetric_Relaxation(std::vector<TInt>& nodesAdded, const std::vector<TInt>& edge, const std::vector<double>& edgeU, const double sizeEdge, const unsigned int edgeId)
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
  //std::vector<DataEdges> dataEdges = subdivideEdge(edge_, edgeU_, sizeEdge_);
  //compute the how many node should be in the edge  (depending on the Metric define on the mesh 's edges')
  //for(auto const &d : dataEdges)
  {
    double den = 0.0;
    /*std::vector<double> edgeU = d.subEdgeU;
    std::vector<TInt> edge = d.subEdge;
    double sizeEdge = d.sizeEdge;*/

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
      //std::cout << "MA_modified -> " << std::endl << MA << std::endl;
      //std::cout << "MB_modified -> " << std::endl << MB << std::endl;

      //use it if you want to add curvature information for the metric computation
      //for the edge subdivision
      //find the metric associate to the node (because of the curvature of the edge)
      /*double curvatureValueA = std::abs(curvature(tA, edge));
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
  }*/
  Eigen::Matrix3d MA_modified;
  Eigen::Matrix3d MB_modified;

  MA_modified(0,0) = 1.0 /sqrt(MA(0,0)); MA_modified(1,1) = 1.0 /sqrt(MA(1,1)); MA_modified(2,2) = 1.0 /sqrt(MA(2,2));
  MB_modified(0,0) = 1.0 /sqrt(MB(0,0)); MB_modified(1,1) = 1.0 /sqrt(MB(1,1)); MB_modified(2,2) = 1.0 /sqrt(MB(2,2));
  //std::cout << "MA_modified -> " << std::endl << MA_modified << std::endl;
  //std::cout << "MB_modified -> " << std::endl << MB_modified << std::endl;


  const Eigen::Vector3d vecA = MA_modified * dir;
  const Eigen::Vector3d vecB = MB_modified * dir;

  const double mA = vecA.norm();
  const double mB = vecB.norm();

  const double sizeInterval = edgeU[i + 1] - edgeU[i];
  //const double sizeInterval = 1.0;

  //https://fr.wikipedia.org/wiki/Calcul_num%C3%A9rique_d%27une_int%C3%A9grale
  //basic discret square integral
  //den += sizeInterval  *  (0.5 * mA + 0.5 * mB);
  den += (tB - tA)  *  (0.5 * mA + 0.5 * mB);
}
if(den == 0.0)
{
  throw gmds::GMDSException("den == 0.0");
}

unsigned int n = std::ceil(sizeEdge / den) + 1;
std::vector<double>res{};
/*while(!metricSamplingEdge(n, res, edge, edgeU)){
++n;
}*/
metricSamplingEdge(n, res, edge, edgeU);
for(unsigned int i = 0 ; i < res.size()  ; i++)
{
  const double u = res[i];
  TInt nodeA; TInt nodeB;
  Point pt = computeTheEdgeNodeCoordinate(u, edge, edgeU, nodeA, nodeB);
  std::vector<TSimplexID> tetraContenaingPt{};
  bool status = false;
  TInt newNodeId = m_nodesMesh.addNodeAndcheck(pt, tetraContenaingPt, status);
  if(!status){
    addNodeToLayer(newNodeId);
    std::vector<TSimplexID> shell = SimplicesNode(m_simplexMesh, nodeA).shell(SimplicesNode(m_simplexMesh, nodeB));
    if(shell.size() == 0)
    throw gmds::GMDSException("shell.size() == 0 in subdivideEdgeUsingMetric_Relaxation function");

    for(auto const s : shell){
      if(s >= 0)
      {
        nodeBelongingTO[newNodeId] = s;
        break;
      }
    }
    (*BND_CURVE_COLOR_NODE)[newNodeId] = edgeId;
    m_nodesMesh.setAnalyticMetricFromMesh(newNodeId, m_simplexMesh, nodeBelongingTO);
    nodesAdded.push_back(newNodeId);
    std::unordered_set<TInt> seen{};
    m_nodesMesh.getOctree()->addNode(newNodeId, seen);
  }
  m_listEdge[edgeId].push_back(newNodeId);
}
}
}
/*----------------------------------------------------------------------------*/
bool MetricFFPointgeneration::metricSamplingEdge(const unsigned int n, std::vector<double>& res, const std::vector<TInt>& edge, const std::vector<double>& edgeU) const
{
  if(edge.size() < 2)
  {
    throw gmds::GMDSException("edge.size < 2");
  }
  double epsilon = 1E-2;
  double minimumPercentageError = 0.1;
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
  double error_moy = 0.0;
  do
  {
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

      TInt nodeA_; TInt nodeB_; //useless here
      const Point pt_min = computeTheEdgeNodeCoordinate(u_min, edge, edgeU, nodeA_, nodeB_);
      const Point pt_max = computeTheEdgeNodeCoordinate(u_max, edge, edgeU, nodeA_, nodeB_);

      Eigen::Vector3d dir = Eigen::Vector3d(pt_max.X() - pt_min.X(), pt_max.Y() - pt_min.Y(), pt_max.Z() - pt_min.Z());
      dir.normalize();

      //std::cout << "OK3" << std::endl;
      bool status0, status1;
      Eigen::Matrix3d M_min  = m_simplexMesh->getAnalyticMetric(pt_min, m_simplexMesh, status0);
      //std::cout << "OK4" << std::endl;
      Eigen::Matrix3d M_max  = m_simplexMesh->getAnalyticMetric(pt_max, m_simplexMesh,status1);

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

      error_moy += std::abs(zk - zk_1) / params_u.size();
    }
    //a checker plus tard
    if(std::abs(error_max_prev - error_moy) / std::abs(error_moy) < minimumPercentageError)
    return false;

    error_max_prev = error_moy;
    p.push_back(1.0);
    params_u = p;
  }while(error_moy > epsilon);
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
