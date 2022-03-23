#include "gmds/hybridMeshAdapt/PointInsertion.h"
#include "gmds/hybridMeshAdapt/SimplexMesh.h"
#include "gmds/hybridMeshAdapt/SimplicesNode.h"
#include "gmds/hybridMeshAdapt/ICriterion.h"
#include "gmds/hybridMeshAdapt/FrontalInsertion.h"
/******************************************************************************/
#include <math.h>       /* sqrt */
/******************************************************************************/
using namespace gmds;
using namespace hybrid;
using namespace simplicesNode;
using namespace operators;

void FrontalInsertion::execute()
{
  if(m_simplexMesh != nullptr)
  {
    //use the node metric and metric lenght between 2 ridge nodes in order
    //to increase/decrease the ridge's node sampling
    std::set<TInt> ridgesLabel{};

    gmds::Variable<Eigen::Matrix3d>* var = nullptr;
    Variable<int>* BND_CURVE_COLOR       = nullptr;

    try{
      BND_CURVE_COLOR =  m_simplexMesh->getVariable<TInt,SimplicesNode>("BND_CURVE_COLOR");
      var =  m_simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
    } catch(gmds::GMDSException e)
    {
      throw gmds::GMDSException(e);
    }
    const std::multimap<TInt, std::pair<TInt,TInt>>& edgeConstStructure = m_simplexMesh->getConstEdgeStructure();

    for(auto const & map : edgeConstStructure)
    {
        ridgesLabel.insert(map.first);
    }

    for(auto const ridgeLabel : ridgesLabel)
    {
      bool status = false;
      std::multimap<TInt, std::pair<TInt,TInt>>& edgeStructure = m_simplexMesh->getEdgeStructure();
      auto it = edgeStructure.equal_range(ridgeLabel);
      for(auto itr = it.first ; itr != it.second ; itr++)
      {
        const SimplicesNode nodeA(m_simplexMesh, itr->second.first);
        const SimplicesNode nodeB(m_simplexMesh, itr->second.second);
        std::cout << nodeA.getGlobalNode() << std::endl;
        std::cout << nodeB.getGlobalNode() << std::endl;
        //condition sur les metriques
        //lorsque la longueur est trop grande par rapport a la metrique local on insert un node aux millieux de cette arete
        //et on interpole la metrique au nouveau node inseré
        math::Vector3d _AB = nodeB.getCoords() - nodeA.getCoords();
        Eigen::Vector3d AB = Eigen::Vector3d(_AB.X(), _AB.Y(), _AB.Z());
        double l  = AB.norm();
        double lm = lenghtMetric(nodeA, nodeB);

        std::cout << "lm -> " << lm << std::endl;
        std::cout << "l -> " << l << std::endl;
        if(l >= lm)
        {
          std::vector<TSimplexID> tetraContenaingPt{};
          bool alreadyAdd = false;
          math::Point  point = (nodeA.getCoords() + nodeB.getCoords()) * 0.5;
          TInt node = m_simplexMesh->addNodeAndcheck(point, tetraContenaingPt, alreadyAdd);
          CriterionRAIS criterionRAIS(new VolumeCriterion());
          std::vector<TSimplexID> shell = nodeA.shell(nodeB);
          gmds::BitVector nodesAdded(m_simplexMesh->nodesCapacity());
          std::vector<TInt> deletedNodes{};
          const std::multimap<TInt, TInt> facesAlreadyBuilt{};

          BND_CURVE_COLOR->set(node, ridgeLabel);
          PointInsertion(m_simplexMesh, SimplicesNode(m_simplexMesh, node), criterionRAIS, status, shell,nodesAdded, deletedNodes, facesAlreadyBuilt);
          var->set(node, (*var)[itr->second.first]);
          std::cout << "OKOK" << std::endl;
          if(status){
            it = edgeStructure.equal_range(ridgeLabel);
          }
        }
        else
        {
          throw gmds::GMDSException("ne doit pas rentrer");
        }
      }
    }
  }
}

double FrontalInsertion::lenghtMetric(const SimplicesNode& nodeA, const SimplicesNode& nodeB)
{
  gmds::Variable<Eigen::Matrix3d>* var;
  try{
    var =  m_simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
  } catch(gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  double iterationMax = 2.0;
  double l = 0.0;
  math::Vector3d _AB = nodeB.getCoords() - nodeA.getCoords();
  Eigen::Vector3d AB = Eigen::Vector3d(_AB.X(), _AB.Y(), _AB.Z());
  //for(unsigned int iter = 0 ; iter < iterationMax ; iter++)
  //{
  //    l += math::sqrt()
  //}
  Eigen::Matrix3d mA = (*var)[nodeA.getGlobalNode()];
  Eigen::Matrix3d mB = (*var)[nodeB.getGlobalNode()];
  l = sqrt(AB.dot(mA* AB)) + sqrt(AB.dot(mB * AB));
  l / iterationMax;
  return l;
}
