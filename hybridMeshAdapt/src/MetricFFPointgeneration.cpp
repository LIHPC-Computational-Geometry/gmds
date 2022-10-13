#include "gmds/hybridMeshAdapt/MetricFFPointgeneration.h"
#include "gmds/hybridMeshAdapt/SimplexMesh.h"
/*----------------------------------------------------------------------------*/
#include <gmds/io/VTKWriter.h>
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
  m_nodesMesh.newVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  m_nodesMesh.newVariable<int,SimplicesNode>("BND_CURVE_COLOR");
}
/*----------------------------------------------------------------------------*/
MetricFFPointgeneration::~MetricFFPointgeneration()
{

}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::execute()
{
  CriterionRAIS criterionRAIS(new VolumeCriterion());

  Variable<int>* BND_VERTEX_COLOR    = nullptr;
  Variable<int>* BND_SURFACE_COLOR   = nullptr;
  Variable<int>* BND_CURVE_COLOR     = nullptr;
  Variable<Eigen::Matrix3d>* metric  = nullptr;

  try{
    metric = m_simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
    BND_VERTEX_COLOR  = m_simplexMesh->getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
    BND_SURFACE_COLOR  = m_simplexMesh->getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
    BND_CURVE_COLOR  = m_simplexMesh->getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
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
    //dichotomie algo reduction edge to see if the edge is sample enough
    //subdivideEdgeUsingMetric_Dichotomie(nodeAdded, edge, edgeU, edge_length);
    subdivideEdgeUsingMetric_Relaxation(nodeAdded, edge, edgeU, edge_length, edgeId);
    i++;
  }

  while(nodeAdded.size() != 0)
  {
    std::cout << "nodeAdded.size() BEFORE -> " << nodeAdded.size() << std::endl;
    nodesSpreading(nodeAdded);
    std::cout << "nodeAdded.size() AFTER  -> " << nodeAdded.size() << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
  }

  for(auto const & d : m_nodeStructure)
  {
    std::vector<TInt> nodes{};
    for(auto const & n : d.second)
    {
      if(n != -1)
      {
        nodes.push_back(n);
      }
    }
    for(auto const n : nodes)
    {
      m_nodesMesh.addTriangle(d.first,n,n);
    }
  }


  const gmds::BitVector& nodesMesh = m_nodesMesh.getBitVectorNodes();
  for(unsigned int n = 0 ; n < nodesMesh.capacity() ; n++)
  {
    if(nodesMesh[n] != 0)
    {
      m_nodesMesh.addTetraedre(n,n,n,n);
    }
  }

  gmds::ISimplexMeshIOService ioService(&m_nodesMesh);
  gmds::VTKWriter vtkWriterMA(&ioService);
  vtkWriterMA.setCellOptions(gmds::N|gmds::R|gmds::F);
  vtkWriterMA.setDataOptions(gmds::N|gmds::R|gmds::F);
  vtkWriterMA.write("metricFF_Node.vtk");
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
void MetricFFPointgeneration::nodesSpreading(std::vector<TInt>& nodesAdded)
{
  static int cpt = 0;
  std::vector<TInt> newNodes{};
  //test direction, we force the direction tospread the nodes to i,j,k
  /*Eigen::Vector3d x_dir(1.0, 0.0, 0.0);
  Eigen::Vector3d y_dir(0.0, 1.0, 1.0);
  Eigen::Vector3d z_dir(0.0, -1.0, 1.0);
  x_dir.normalize(); y_dir.normalize(); z_dir.normalize();*/
  const Eigen::Vector3d x_dir(1.0, 0.0, 0.0);
  const Eigen::Vector3d y_dir(0.0, 1.0, 0.0);
  const Eigen::Vector3d z_dir(0.0, 0.0, 1.0);
  std::vector<double> spins{1.0, -1.0};

  std::vector<Eigen::Vector3d> dirs{x_dir, y_dir, z_dir};
  Variable<Eigen::Matrix3d>* metric  = nullptr;
  try{
    metric = m_simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  for(auto const node : nodesAdded)
  {
    unsigned int id = 0;
    std::vector<TInt> v{-1, -1, -1, -1, -1, -1};
    const math::Point pt = SimplicesNode(&m_nodesMesh, node).getCoords();
    for(auto const & dir : dirs)
    {
      for(auto const spin : spins)
      {
        Eigen::Matrix3d M = metric->value(node);
        M(0,0) = 1.0 /sqrt(M(0,0)); M(1,1) = 1.0 /sqrt(M(1,1)); M(2,2) = 1.0 /sqrt(M(2,2));
        Eigen::Vector3d d = spin * dir;
        d.normalize();
        const Eigen::Vector3d vec = M * d;
        const double m = vec.norm();

        math::Point newCoord = pt + m*math::Point(d.x(), d.y(), d.z());
        //check if the pt is in the mesh
        if(m_oc.belongToOc(newCoord))
        {
          //std::cout << "node belong to oc " << std::endl;
          //compute the metric lentgh based on newNodeId and the metric at node
          findOptimimalPosition(node, newCoord);
          std::vector<TInt> neighboorNodes;
          nodeFiltering(node, newCoord, neighboorNodes);

          if(!neighboorNodes.size())
          {
            TInt newNodeId = m_nodesMesh.addNode(newCoord);
            v[id] = newNodeId;
            newNodes.push_back(newNodeId);
            gmds::ISimplexMeshIOService ioService(&m_nodesMesh);
            gmds::VTKWriter vtkWriterMA(&ioService);
            vtkWriterMA.setCellOptions(gmds::N);
            vtkWriterMA.setDataOptions(gmds::N);
            vtkWriterMA.write("metricFF_Node_CPT_" + std::to_string(cpt) + ".vtk");
            cpt++;
          }
          else
          {
            double distMin = std::numeric_limits<double>::max();
            TInt neighboorNode;
            for(auto const n : neighboorNodes)
            {
              math::Point pid = SimplicesNode(&m_nodesMesh, n).getCoords();
              math::Vector3d v = math::Vector3d(pid.X() - pt.X(), pid.Y() - pt.Y(), pid.Z() - pt.Z());
              double dist = v.norm();
              if(dist < distMin)
              {
                distMin = dist;
                neighboorNode = n;
              }
            }
            v[id] = neighboorNode;
          }
        }
        id++;
      }
    }
    m_nodeStructure.insert( std::pair<TInt, std::vector<TInt>>(node, v) );
  }

  nodesAdded.clear();
  nodesAdded = newNodes;
}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::nodeFiltering(const TInt node, const math::Point& pt, std::vector<TInt> & neighboorNode)
{
  const double k = 0.8;
  const double epsilon = 0.01;
  Variable<Eigen::Matrix3d>* metric  = nullptr;
  gmds::Variable<int>* BND_CURVE_COLOR_NODE = nullptr;

  try{
    metric = m_simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
    BND_CURVE_COLOR_NODE = m_nodesMesh.getVariable<int, SimplicesNode>("BND_CURVE_COLOR");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  bool flag = false;
  std::vector<math::Point> pts{};
  const gmds::BitVector& nodesMesh = m_nodesMesh.getBitVectorNodes();
  for(unsigned int nodeId = 0 ; nodeId < nodesMesh.capacity() ; nodeId++)
  {
    if(nodesMesh[nodeId] != 0 && nodeId != node)
    {
      const math::Point nodeCoord = SimplicesNode(&m_nodesMesh, nodeId).getCoords();
      Eigen::Vector3d v(pt.X() - nodeCoord.X(), pt.Y() - nodeCoord.Y(), pt.Z() - nodeCoord.Z());
      //the metric being analytique we do not have to interpolate the metric at point
      const Eigen::Matrix3d m0 = (*metric)[nodeId];
      const Eigen::Matrix3d m1 = m_nodesMesh.getAnalyticMetric(pt);
      //compute the current length based on the metric atached to the mesh
      const double metricLenght = 0.5 * sqrt(v.dot(m0*v)) +  0.5 * sqrt(v.dot(m1*v));
      if(metricLenght <= k*(1.0 - epsilon))
      {
        neighboorNode.push_back(nodeId);
      }
    }
  }
}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::findOptimimalPosition(const TInt node, math::Point& initialCoord)
{
  double epsilon = 0.001;
  Variable<Eigen::Matrix3d>* metric  = nullptr;
  try{
    metric = m_simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  //the node is the layer node so we can not move this once
  const math::Point pt = SimplicesNode(&m_nodesMesh, node).getCoords();
  Eigen::Vector3d dir(initialCoord.X() - pt.X(), initialCoord.Y() - pt.Y(), initialCoord.Z() - pt.Z());

  //the metric being analytique we do not have to interpolate the metric at point
  const Eigen::Matrix3d m0 = (*metric)[node];
  const Eigen::Matrix3d m1 = m_nodesMesh.getAnalyticMetric(initialCoord);

  //compute the current length based on the metric atached to the mesh
  const double metricLenght = 0.5 * sqrt(dir.dot(m0*dir)) +  0.5 * sqrt(dir.dot(m1*dir));

  if(metricLenght < 1.0 - epsilon)
  {
    initialCoord = pt + 1.5 * math::Vector3d(dir.x(), dir.y(), dir.z());
    findOptimimalPosition(node, initialCoord);
  }
  else if(metricLenght > 1.0 + epsilon)
  {
    initialCoord = 0.5 * (pt + initialCoord);
    findOptimimalPosition(node, initialCoord);
  }
}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::subdivideEdgeUsingMetric_Relaxation(std::vector<TInt>& nodesAdded, const std::vector<TInt>& edge, const std::vector<double>& edgeU, const double sizeEdge, const unsigned int edgeId)
{
  double den = 0.0;
  Variable<Eigen::Matrix3d>* metric  = nullptr;
  gmds::Variable<int>* BND_CURVE_COLOR_NODE = nullptr;

  try{
    metric = m_simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
    BND_CURVE_COLOR_NODE = m_nodesMesh.getVariable<int, SimplicesNode>("BND_CURVE_COLOR");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  //compute the how many node should be in the edge  (depending on the Metric define on the mesh 's edges')
  for(unsigned int i = 0 ; i < edge.size() - 1 ; i++)
  {
    const TInt nodeA = edge[i];
    const TInt nodeB = edge[i + 1];

    const math::Point ptA = SimplicesNode(m_simplexMesh, nodeA).getCoords();
    const math::Point ptB = SimplicesNode(m_simplexMesh, nodeB).getCoords();

    Eigen::Vector3d dir = Eigen::Vector3d(ptB.X() - ptA.X(), ptB.Y() - ptA.Y(), ptB.Z() - ptA.Z());
    dir.normalize();

    const Eigen::Matrix3d MA = metric->value(nodeA);
    const Eigen::Matrix3d MB = metric->value(nodeB);

    Eigen::Matrix3d MA_modified = metric->value(nodeA);
    Eigen::Matrix3d MB_modified = metric->value(nodeB);

    MA_modified(0,0) = 1.0 /sqrt(MA(0,0)); MA_modified(1,1) = 1.0 /sqrt(MA(1,1)); MA_modified(2,2) = 1.0 /sqrt(MA(2,2));
    MB_modified(0,0) = 1.0 /sqrt(MB(0,0)); MB_modified(1,1) = 1.0 /sqrt(MB(1,1)); MB_modified(2,2) = 1.0 /sqrt(MB(2,2));

    const Eigen::Vector3d vecA = MA_modified * dir;
    const Eigen::Vector3d vecB = MB_modified * dir;

    const double mA = vecA.norm();
    const double mB = vecB.norm();

    const double sizeInterval = edgeU[i + 1] - edgeU[i];

    //https://fr.wikipedia.org/wiki/Calcul_num%C3%A9rique_d%27une_int%C3%A9grale
    //basic discret square integral
    den += sizeInterval *  (0.5 * mA + 0.5 * mB);
  }
  if(den == 0.0)
  {
    throw gmds::GMDSException("den == 0.0");
  }

  const unsigned int n = std::ceil(static_cast<double>(sizeEdge) / den);
  std::vector<double>res{};
  metricSamplingEdge(n, res, edge, edgeU);

  for(unsigned int i = 0 ; i < res.size()  ; i++)
  {
    const double u = res[i];
    Point pt = computeTheEdgeNodeCoordinate(u, edge, edgeU);
    TInt newNodeId = m_nodesMesh.addNode(pt);
    (*BND_CURVE_COLOR_NODE)[newNodeId] = edgeId;
    m_nodesMesh.setAnalyticMetric(newNodeId);
    nodesAdded.push_back(newNodeId);
  }
}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::metricSamplingEdge(const unsigned int n, std::vector<double>& res, const std::vector<TInt>& edge, const std::vector<double>& edgeU) const
{
  if(edge.size() < 2)
  {
    throw gmds::GMDSException("edge.size < 2");
  }
  double epsilon = 0.01;
  //first we subdivide the edge uniformly (whitout using any random generator for the position)
  std::vector<double> params_u{};
  const double e = 1.0 / static_cast<double>(n);
  for(unsigned int i = 0 ; i < n + 1 ; i++)
  {
    params_u.push_back(e * static_cast<double>(i));
  }
  //Compute the right sample node position in u coordinate with Loyd methodes in CVT -> Fast Methods for Computing Centroidal Voronoi Tessellations
  //compute V(zk) & the corresponding zk -> [0, 1]
  double error_max = 0.0;
  do
  {
    error_max = 0.0;
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

      const Point pt_min = computeTheEdgeNodeCoordinate(u_min, edge, edgeU);
      const Point pt_max = computeTheEdgeNodeCoordinate(u_max, edge, edgeU);

      Eigen::Vector3d dir = Eigen::Vector3d(pt_max.X() - pt_min.X(), pt_max.Y() - pt_min.Y(), pt_max.Z() - pt_min.Z());
      dir.normalize();

      Eigen::Matrix3d M_min  = m_simplexMesh->getAnalyticMetric(pt_min);
      Eigen::Matrix3d M_max  = m_simplexMesh->getAnalyticMetric(pt_max);

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
      error_max = std::max(error_max, std::abs(zk - params_u[i + 1]));
    }
    p.push_back(1.0);
    params_u = p;
  }while(error_max > epsilon);

  res = params_u;
}
/*----------------------------------------------------------------------------*/
void MetricFFPointgeneration::subdivideEdgeUsingMetric_Dichotomie(std::vector<TInt>& nodesAdded, const std::vector<TInt>& edge, const std::vector<double>& edgeU, const double sizeEdge) const
{
  if(edge.size() < 2)
  {
    return;
  }

  double den = 0.0;
  Variable<Eigen::Matrix3d>* metric  = nullptr;

  try{
    metric = m_simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  for(unsigned int i = 0 ; i < edge.size() - 1 ; i++)
  {
    const TInt nodeA = edge[i];
    const TInt nodeB = edge[i + 1];

    const math::Point ptA = SimplicesNode(m_simplexMesh, nodeA).getCoords();
    const math::Point ptB = SimplicesNode(m_simplexMesh, nodeB).getCoords();

    Eigen::Vector3d dir = Eigen::Vector3d(ptB.X() - ptA.X(), ptB.Y() - ptA.Y(), ptB.Z() - ptA.Z());
    dir.normalize();

    const Eigen::Matrix3d MA = metric->value(nodeA);
    const Eigen::Matrix3d MB = metric->value(nodeB);

    Eigen::Matrix3d MA_modified = metric->value(nodeA);
    Eigen::Matrix3d MB_modified = metric->value(nodeB);

    MA_modified(0,0) = 1.0 /sqrt(MA(0,0)); MA_modified(1,1) = 1.0 /sqrt(MA(1,1)); MA_modified(2,2) = 1.0 /sqrt(MA(2,2));
    MB_modified(0,0) = 1.0 /sqrt(MB(0,0)); MB_modified(1,1) = 1.0 /sqrt(MB(1,1)); MB_modified(2,2) = 1.0 /sqrt(MB(2,2));

    const Eigen::Vector3d vecA = MA_modified * dir;
    const Eigen::Vector3d vecB = MB_modified * dir;

    const double mA = vecA.norm();
    const double mB = vecB.norm();

    const double sizeInterval = edgeU[i + 1] - edgeU[i];

    //https://fr.wikipedia.org/wiki/Calcul_num%C3%A9rique_d%27une_int%C3%A9grale
    //basic discret intergral square integral
    den += sizeInterval *  (0.5 * mA + 0.5 * mB);
  }

  if(den == 0.0)
  {
    throw gmds::GMDSException("den == 0.0");
  }
  else if(sizeEdge / den > 1.0)
  {
    //compute the center of the parametrized Edge usinig U
    for(unsigned int i = 0 ; i < edgeU.size() - 1 ; i++)
    {
      const double uA = edgeU[i];
      const double uB = edgeU[i + 1];
      if(uA <= 0.5 && uB >= 0.5)
      {
        //interpolation (linear) of the middle position using the interval ua, ub and u = 0.5
        if(uB != uA)
        {
          const TInt nodeA = edge[i];
          const TInt nodeB = edge[i + 1];

          const Point ptA = SimplicesNode(m_simplexMesh, nodeA).getCoords();
          const Point ptB = SimplicesNode(m_simplexMesh, nodeB).getCoords();

          const double t = (0.5-uA) / (uB-uA);

          std::vector<TInt> newEdge0{};
          std::vector<TInt> newEdge1{};

          std::vector<double> edgeU_0 {};
          std::vector<double> edgeU_1 {};

          if(t == 0.0)
          {
            nodesAdded.push_back(nodeA);
            std::copy(edge.begin(), edge.begin() + i + 1, std::back_inserter(newEdge0));
            std::copy(edge.begin() + i, edge.end(), std::back_inserter(newEdge1));
          }
          else if(t == 1.0)
          {
            nodesAdded.push_back(nodeB);
            std::copy(edge.begin(), edge.begin() + i + 2, std::back_inserter(newEdge0));
            std::copy(edge.begin() + i + 1, edge.end(), std::back_inserter(newEdge1));
          }
          else
          {
            const Point pt = ptA * (1.0 - t) + ptB * t;
            std::vector<TSimplexID> tetraContenaingPt{};
            bool alreadyAdd = false;
            TInt newNodeId = m_simplexMesh->addNodeAndcheck(pt, tetraContenaingPt, alreadyAdd);
            m_simplexMesh->setAnalyticMetric(newNodeId);
            std::copy(edge.begin(), edge.begin() + i + 1, std::back_inserter(newEdge0));
            newEdge0.push_back(newNodeId);
            newEdge1.push_back(newNodeId);
            std::copy(edge.begin() + i + 2, edge.end(), std::back_inserter(newEdge1));
            nodesAdded.push_back(newNodeId);
          }

          std::vector<double> edges_length0{};
          std::vector<double> edges_length1{};
          std::map<unsigned int, std::vector<TInt>> sortedEdges0;
          std::map<unsigned int, std::vector<TInt>> sortedEdges1;
          sortedEdges0[0] = newEdge0;
          sortedEdges1[0] = newEdge1;
          const std::vector<std::vector<double>> edgesU0 = buildParamEdgeU(sortedEdges0, edges_length0);
          const std::vector<std::vector<double>> edgesU1 = buildParamEdgeU(sortedEdges1, edges_length1);

          subdivideEdgeUsingMetric_Dichotomie(nodesAdded, newEdge0,  edgesU0.front(), edges_length0.front());
          subdivideEdgeUsingMetric_Dichotomie(nodesAdded, newEdge1,  edgesU1.front(), edges_length1.front());
          break;
        }
        else
        {
          throw gmds::GMDSException("uB == uA");
        }
      }
    }
  }
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

  std::vector<std::set<TInt>> sortedEdges{};
  for(unsigned int node = 0 ; node < meshNode.capacity() ; node++)
  {
    if((*BND_VERTEX_COLOR)[node] != 0)
    {
      const std::vector<TInt> neighborNodes = SimplicesNode(m_simplexMesh, node).neighborNodes();
      for(auto const neighborNode : neighborNodes)
      {
        if((*BND_CURVE_COLOR)[neighborNode] != 0 || (*BND_VERTEX_COLOR)[neighborNode] != 0)
        {
          std::set<TInt> sortedEdge{static_cast<int>(node)};
          sortedEdge.insert(neighborNode);
          sortedEdges.push_back(sortedEdge);
        }
      }
    }
  }


  std::vector<std::vector<TInt>> finalEdges{};
  for(auto const & sortedEdge : sortedEdges)
  {
    const TInt firstNode = *(sortedEdge.begin());
    TInt secondNode = *(--sortedEdge.end());
    if((*BND_CURVE_COLOR)[secondNode] != 0)
    {
      std::vector<TInt> finalEdge{firstNode, secondNode};
      do{
        const std::vector<TInt> neighborNodes = SimplicesNode(m_simplexMesh, secondNode).neighborNodes();
        for(auto const node : neighborNodes)
        {
          if((((*BND_VERTEX_COLOR)[node] != 0 && node != firstNode) || ((*BND_CURVE_COLOR)[node] == (*BND_CURVE_COLOR)[secondNode] && node != secondNode)) &&
            std::find(finalEdge.begin(), finalEdge.end(), node) == finalEdge.end())
          {
            finalEdge.push_back(node);
            secondNode = finalEdge.back();
            break;
          }
        }
      }while((*BND_VERTEX_COLOR)[finalEdge.back()] == 0);
      finalEdges.push_back(finalEdge);
    }
    else
    {
      throw gmds::GMDSException("(*BND_CURVE_COLOR)[secondNode] != 0");
    }
  }


  for(auto const & finalEdge : finalEdges)
  {
    res.insert(std::pair<unsigned int, std::vector<TInt>>((*BND_CURVE_COLOR)[finalEdge[1]], finalEdge));
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
