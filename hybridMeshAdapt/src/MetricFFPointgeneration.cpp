#include "gmds/hybridMeshAdapt/PointInsertion.h"
#include "gmds/hybridMeshAdapt/MetricFFPointgeneration.h"
#include "gmds/hybridMeshAdapt/SimplexMesh.h"
#include "gmds/hybridMeshAdapt/SimplicesCell.h"
#include "gmds/hybridMeshAdapt/SimplicesNode.h"
#include "gmds/hybridMeshAdapt/SimplicesTriangle.h"
#include "gmds/hybridMeshAdapt/Metric.h"
#include <gmds/hybridMeshAdapt/ICriterion.h>
#include <gmds/hybridMeshAdapt/ISimplexMeshIOService.h>
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
MetricFFPointgeneration::MetricFFPointgeneration(SimplexMesh* simplexMesh):m_simplexMesh(simplexMesh)
{

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

  const std::map<unsigned int, std::vector<TInt>> sortedEdges = buildSortedEdges();
  const std::vector<std::vector<double>> edgesU = buildParamEdgeU(sortedEdges);

  unsigned int i = 0;
  for(auto const & sortedEdge : sortedEdges)
  {
    std::vector<double> edgeU = edgesU[i];
    //dichotomie algo reduction edge to see if the edge is sample enough
    
    i++;
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
std::vector<std::vector<double>> MetricFFPointgeneration::buildParamEdgeU(const std::map<unsigned int, std::vector<TInt>>& sortedEdge) const
{
  std::vector<double> sizeEdges{};
  std::vector<std::vector<double>> sizeSubEdges{};

  std::cout << "sortedEdge.size() -> " << sortedEdge.size() << std::endl;
  for(auto const & edge : sortedEdge)
  {
    double AB_length = 0.0;
    std::vector<double> sizeSubEdge{};
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
