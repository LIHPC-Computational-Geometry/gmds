#include "gmds/hybridMeshAdapt/PointInsertion.h"
#include "gmds/hybridMeshAdapt/MetricAdaptation.h"
#include "gmds/hybridMeshAdapt/SimplexMesh.h"
#include "gmds/hybridMeshAdapt/SimplicesCell.h"
#include "gmds/hybridMeshAdapt/SimplicesNode.h"
#include "gmds/hybridMeshAdapt/Metric.h"
#include <gmds/hybridMeshAdapt/ICriterion.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
using namespace math;
using namespace hybrid;
using namespace operators;
using namespace simplicesNode;
using namespace simplicesCell;
/*----------------------------------------------------------------------------*/
MetricAdaptation::MetricAdaptation(SimplexMesh* simplexMesh):m_simplexMesh(simplexMesh)
{
  buildEdgesMap();
  metricCorrection();
}
/*----------------------------------------------------------------------------*/
MetricAdaptation::~MetricAdaptation()
{

}
/*----------------------------------------------------------------------------*/
void MetricAdaptation::metricCorrection()
{

}
/*----------------------------------------------------------------------------*/
void MetricAdaptation::buildEdgesMap()
{
  //generate the edge structure in order to iterate over edge and not the mesh's node
  const gmds::BitVector& meshNode = m_simplexMesh->getBitVectorNodes();
  for(unsigned int nodeId = 0 ; nodeId < meshNode.capacity() ; nodeId++)
  {
    if(meshNode[nodeId] != 0)
    {
      SimplicesNode node(m_simplexMesh, nodeId);
      const std::vector<TInt> && directNodes = node.directNeighboorNodeId();
      for(auto const dNode : directNodes)
      {
        if(dNode < nodeId){m_edgesMap.insert(std::pair<TInt, TInt>(dNode, nodeId));}
        else{m_edgesMap.insert(std::pair<TInt, TInt>(dNode, nodeId));}
      }
    }
  }
}
/*----------------------------------------------------------------------------*/
void MetricAdaptation::execute()
{
  Variable<Eigen::Matrix3d>* metric = m_simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");

  if(metric == nullptr)
  {
    gmds::GMDSException("NO METRIC ATTACHED TO THE MESH");
  }
  else
  {
    Variable<int>* BND_VERTEX_COLOR  = m_simplexMesh->getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
    Variable<int>* BND_CURVE_COLOR   = m_simplexMesh->getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
    Variable<int>* BND_SURFACE_COLOR = m_simplexMesh->getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");



    CriterionRAIS criterionRAIS(new VolumeCriterion());
    const gmds::BitVector& meshNode = m_simplexMesh->getBitVectorNodes();
    for(auto const & edges : m_edgesMap)
    {
      const TInt node0 = edges.first;
      const TInt node1 = edges.second;
      if(meshNode[node0] != 0 && meshNode[node1])
      {
        SimplicesNode sNode0(m_simplexMesh, node0);
        SimplicesNode sNode1(m_simplexMesh, node1);

        if(sNode0.shell(sNode1).size() == 0){continue;}
        math::Point nodeCoord0 = sNode0.getCoords();
        math::Point nodeCoord1 = sNode1.getCoords();
        Metric<Eigen::Matrix3d> M0 = Metric<Eigen::Matrix3d>((*metric)[node0]);
        Metric<Eigen::Matrix3d> M1 = Metric<Eigen::Matrix3d>((*metric)[node1]);
        double metricLenght  = M0.metricDist(nodeCoord0, nodeCoord1, M1);

        const Vector3d vec =  nodeCoord0 - nodeCoord1;
        //config data for the point insertion algorithm
        double lenght_l2 = vec.norm();
        bool status = false;
        const gmds::BitVector markedNodes{};
        std::vector<TSimplexID> deletedSimplex{};
        std::vector<TInt> deletedNodes{};
        const std::multimap<TInt, std::pair<TInt, TInt>> facesAlreadyBuilt{};

        if(lenght_l2 < metricLenght)
        {
          //lenght_l2 < metricLenght  mean vec is too small so we build anew edge taller than vec
          const std::vector<TSimplexID>&& ball = sNode1.ballOf();
          //SimplicesNode node10(m_simplexMesh, 10);
          PointInsertion pi(m_simplexMesh, sNode0, criterionRAIS, status, ball, markedNodes, deletedNodes, facesAlreadyBuilt);
          //return;
        }
        else
        {
          //lenght_l2 < metricLenght means the lenght is too long for the mesh's metric so we move dNode in the middle of the current edge
          //(spliting by 2 the size of the edge)
          std::vector<TSimplexID> tetraContenaingPt{};
          bool alreadyAdd = false;
          const Point pt = 0.5 * (nodeCoord0 + nodeCoord1);
          TInt newNodeId = m_simplexMesh->addNodeAndcheck(pt, tetraContenaingPt, alreadyAdd);

          unsigned int nodeDim = ((*BND_VERTEX_COLOR)[node0] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[node0] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[node0] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
          unsigned int nodeLabel = ((*BND_VERTEX_COLOR)[node0] != 0)?(*BND_VERTEX_COLOR)[node0]:((*BND_CURVE_COLOR)[node0] != 0)?(*BND_CURVE_COLOR)[node0]:((*BND_SURFACE_COLOR)[node0] != 0)?(*BND_SURFACE_COLOR)[node0]:0;
          unsigned int directNodeDim = ((*BND_VERTEX_COLOR)[node1] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[node1] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[node1] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
          unsigned int directNodeLabel = ((*BND_VERTEX_COLOR)[node1] != 0)?(*BND_VERTEX_COLOR)[node1]:((*BND_CURVE_COLOR)[node1] != 0)?(*BND_CURVE_COLOR)[node1]:((*BND_SURFACE_COLOR)[node1] != 0)?(*BND_SURFACE_COLOR)[node1]:0;

          unsigned int newNodeDim    = std::max(nodeDim, directNodeDim);

          //labelization of the node being inserted
          if(newNodeDim == SimplexMesh::topo::CORNER)
          {
            //TODO
            //define a new curve color...
            //BND_CURVE_COLOR->set(newNodeId, newColorToDefine)
          }
          else if(newNodeDim == SimplexMesh::topo::RIDGE)
          {
            if(nodeDim == directNodeDim)
            {
              if(nodeLabel != directNodeLabel)
              {
                //node sur la surface...
                break;
              }
              else
              {
                BND_CURVE_COLOR->set(newNodeId, nodeLabel);
              }
            }
            else if(nodeDim != directNodeDim)
            {
              if(nodeDim > directNodeDim)
              {
                BND_CURVE_COLOR->set(newNodeId, nodeLabel);
              }
              else
              {
                BND_CURVE_COLOR->set(newNodeId, directNodeLabel);
              }
            }
          }
          else if(newNodeDim == SimplexMesh::topo::SURFACE)
          {
            if(nodeDim == directNodeDim)
            {
              if(nodeLabel == directNodeLabel)
              {
                BND_SURFACE_COLOR->set(newNodeId, nodeLabel);
              }
              //node dans le volume donc RAF sinon
            }
            else if(nodeDim != directNodeDim)
            {
              if(nodeDim > directNodeDim)
              {
                BND_SURFACE_COLOR->set(newNodeId, nodeLabel);
              }
              else
              {
                BND_SURFACE_COLOR->set(newNodeId, directNodeLabel);
              }
            }
          }

          const std::vector<TSimplexID>&& shell = sNode0.shell(sNode1);
          PointInsertion pi(m_simplexMesh, SimplicesNode(m_simplexMesh, newNodeId), criterionRAIS, status, shell, markedNodes, deletedNodes, facesAlreadyBuilt);
        }
      }
    }
  }
}
