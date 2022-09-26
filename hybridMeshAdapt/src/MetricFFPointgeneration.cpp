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

  const gmds::BitVector& meshNode = m_simplexMesh->getBitVectorNodes();
  const std::map<unsigned int, std::vector<unsigned int>>& C2E = m_simplexMesh->getCornerEdgeConnexion();
  const std::multimap<TInt, std::pair<TInt,TInt>>& edgeStructure = m_simplexMesh->getEdgeStructure();
  for(auto const & corner : C2E)
  {
    std::vector<unsigned int> edges = corner.second;
    for(auto const edges)
    if()
    {

    }
  }
}
