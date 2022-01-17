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
  metricCorrection();
}


MetricAdaptation::~MetricAdaptation()
{

}

void MetricAdaptation::metricCorrection()
{

}

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
    for(unsigned int nodeId = 0 ; nodeId < meshNode.capacity() ; nodeId++)
    {
      if(meshNode[nodeId] != 0)
      {
        nodesAlreadyDone.resize(meshNode.capacity());
        SimplicesNode node(m_simplexMesh, nodeId);
        const std::vector<TInt> && directNodes = node.directNeighboorNodeId();
        math::Point nodeCoord = node.getCoords();
        Metric<Eigen::Matrix3d> M0 = Metric<Eigen::Matrix3d>((*metric)[nodeId]);
        for(auto const directNode : directNodes)
        {
          SimplicesNode dNode(m_simplexMesh, directNode);
          math::Point dNodeCoord = dNode.getCoords();
          Metric<Eigen::Matrix3d> M1 = Metric<Eigen::Matrix3d>((*metric)[directNode]);
          double metricLenght  = M0.metricDist(nodeCoord, dNodeCoord, M1);

          const Vector3d vec =  nodeCoord - dNodeCoord;
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
            const std::vector<TSimplexID>&& ball = dNode.ballOf();
            PointInsertion pi(m_simplexMesh, node, criterionRAIS, status, ball, markedNodes, deletedNodes, facesAlreadyBuilt);
            if(status)
            {
              nodesAlreadyDone.assign(node.getGlobalNode());
              //break;
            }
          }
          else
          {
            //lenght_l2 < metricLenght means the lenght is too long for the mesh's metric so we move dNode in the middle of the current edge
            //(spliting by 2 the size of the edge)
            std::vector<TSimplexID> tetraContenaingPt{};
            bool alreadyAdd = false;
            const Point pt = 0.5 * (dNodeCoord + nodeCoord);
            TInt newNodeId = m_simplexMesh->addNodeAndcheck(pt, tetraContenaingPt, alreadyAdd);

            unsigned int nodeDim = ((*BND_VERTEX_COLOR)[nodeId] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[nodeId] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[nodeId] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
            unsigned int nodeLabel = ((*BND_VERTEX_COLOR)[nodeId] != 0)?(*BND_VERTEX_COLOR)[nodeId]:((*BND_CURVE_COLOR)[nodeId] != 0)?(*BND_CURVE_COLOR)[nodeId]:((*BND_SURFACE_COLOR)[nodeId] != 0)?(*BND_SURFACE_COLOR)[nodeId]:0;
            unsigned int directNodeDim = ((*BND_VERTEX_COLOR)[directNode] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[directNode] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[directNode] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
            unsigned int directNodeLabel = ((*BND_VERTEX_COLOR)[directNode] != 0)?(*BND_VERTEX_COLOR)[directNode]:((*BND_CURVE_COLOR)[directNode] != 0)?(*BND_CURVE_COLOR)[directNode]:((*BND_SURFACE_COLOR)[directNode] != 0)?(*BND_SURFACE_COLOR)[directNode]:0;

            unsigned int newNodeDim    = std::max(nodeDim, directNodeDim);

            //labelization of the node being inserted
            if(newNodeDim == SimplexMesh::topo::CORNER)
            {
              //TODO
              //define a new curve color...
              //BND_CURVE_COLOR->set(newNodeId, newColorToDefine)
              break;
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

            //BND_CURVE_COLOR->set(node, (*BND_CURVE_COLOR_NODES)[idx])
            const std::vector<TSimplexID>&& shell = node.shell(dNode);
            PointInsertion pi(m_simplexMesh, SimplicesNode(m_simplexMesh, newNodeId), criterionRAIS, status, shell, markedNodes, deletedNodes, facesAlreadyBuilt);
            if(status)
            {
              break;
            }
          }
        }
      }
    }
  }
}
