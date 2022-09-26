#include "gmds/hybridMeshAdapt/PointInsertion.h"
#include "gmds/hybridMeshAdapt/MetricAdaptation.h"
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
MetricAdaptation::MetricAdaptation(SimplexMesh* simplexMesh):m_simplexMesh(simplexMesh)
{
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
  m_edgesMap.clear();
  //generate the edge structure in order to iterate over edge and not the mesh's node
  const gmds::BitVector& meshNode = m_simplexMesh->getBitVectorNodes();
  bool flag = false;
  for(TInt nodeId = 0 ; nodeId < meshNode.capacity() ; nodeId++)
  {
    if(meshNode[nodeId] != 0)
    {
      SimplicesNode node(m_simplexMesh, nodeId);
      const std::vector<TInt> && directNodes = node.directNeighboorNodeId();
      for(auto dNode : directNodes)
      {
        flag = false;
        TInt minNode = std::min(dNode, nodeId);
        std::pair<std::multimap<TInt,TInt>::iterator, std::multimap<TInt,TInt>::iterator> ret = m_edgesMap.equal_range(minNode);
        for(std::multimap<TInt,TInt>::iterator it = ret.first; it != ret.second ; it++)
        {
            if(it->second == std::max(dNode, nodeId))
            {
              flag = true;
              break;
            }
        }
        if(!flag)
        {
          m_edgesMap.insert(std::pair<TInt, TInt>(std::min(dNode, nodeId), std::max(dNode, nodeId)));
          //break;
        }
      }
    }
  }

  /*for(auto const m : m_edgesMap){
    std::cout << m.first << " | " << m.second << std::endl;
  }*/
}
/*----------------------------------------------------------------------------*/
bool MetricAdaptation::computeSlicing(const TInt nodeA, const TInt nodeB)
{
  static unsigned int cptTest = 0;
  Variable<int>* BND_VERTEX_COLOR   = nullptr;
  Variable<int>* BND_CURVE_COLOR    = nullptr;
  Variable<int>* BND_SURFACE_COLOR  = nullptr;
  Variable<Eigen::Matrix3d>* metric = nullptr;
  try{
    BND_VERTEX_COLOR  = m_simplexMesh->getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
    BND_CURVE_COLOR   = m_simplexMesh->getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
    BND_SURFACE_COLOR = m_simplexMesh->getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
    metric            = m_simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");

  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  gmds::ISimplexMeshIOService ioServiceMesh(m_simplexMesh);

  CriterionRAIS criterionRAIS(new VolumeCriterion());
  const gmds::BitVector& meshNode = m_simplexMesh->getBitVectorNodes();

  //Metric adaptation SLICING
  if(meshNode[nodeA] != 0 && meshNode[nodeB] != 0)
  {
    SimplicesNode node0(m_simplexMesh, nodeA);
    SimplicesNode node1(m_simplexMesh, nodeB);
    std::vector<TSimplexID> shell = node0.shell(node1);
    if(shell.size() == 0)
    {
      return false;
    }

    unsigned int nodeDim = ((*BND_VERTEX_COLOR)[nodeA] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[nodeA] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[nodeA] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
    unsigned int nodeLabel = ((*BND_VERTEX_COLOR)[nodeA] != 0)?(*BND_VERTEX_COLOR)[nodeA]:((*BND_CURVE_COLOR)[nodeA] != 0)?(*BND_CURVE_COLOR)[nodeA]:((*BND_SURFACE_COLOR)[nodeA] != 0)?(*BND_SURFACE_COLOR)[nodeA]:0;
    unsigned int directNodeDim = ((*BND_VERTEX_COLOR)[nodeB] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[nodeB] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[nodeB] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
    unsigned int directNodeLabel = ((*BND_VERTEX_COLOR)[nodeB] != 0)?(*BND_VERTEX_COLOR)[nodeB]:((*BND_CURVE_COLOR)[nodeB] != 0)?(*BND_CURVE_COLOR)[nodeB]:((*BND_SURFACE_COLOR)[nodeB] != 0)?(*BND_SURFACE_COLOR)[nodeB]:0;
    unsigned int newNodeDim = std::max(nodeDim, directNodeDim);

    const std::map<unsigned int, std::pair<unsigned int, unsigned int>>& E2S_indices = m_simplexMesh->getEdgeTianglesIndices();
    const std::map<unsigned int, std::vector<unsigned int>>& C2S_indices = m_simplexMesh->getCornerSurfaceConnexion() ;
    const std::map<unsigned int, std::vector<unsigned int>>& C2E_indices = m_simplexMesh->getCornerEdgeConnexion() ;

    bool alreadyAdd = false;
    std::vector<TSimplexID> tetraContenaingPt{};
    math::Point nodeCoord0 = node0.getCoords();
    math::Point nodeCoord1 = node1.getCoords();
    const Point pt = 0.5 * (nodeCoord0 + nodeCoord1);
    TInt newNodeId = m_simplexMesh->addNodeAndcheck(pt, tetraContenaingPt, alreadyAdd);

    if(alreadyAdd)
    {
      //TODO not sure about that we can reinsert the point
      return false;
    }

    m_simplexMesh->setAnalyticMetric(newNodeId);
    //labelization of the node being inserted
    if(newNodeDim == SimplexMesh::topo::CORNER)
    {
      //TODO
      //define a new curve color...
      //BND_CURVE_COLOR->set(newNodeId, newColorToDefine)
      return false;
      /*std::cout << "(*BND_VERTEX_COLOR)[nodeA] -> " << (*BND_VERTEX_COLOR)[nodeA] << std::endl;
      std::cout << "(*BND_VERTEX_COLOR)[nodeB] -> " << (*BND_VERTEX_COLOR)[nodeB] << std::endl;
      std::vector<unsigned int> edgeIdx_A = C2E_indices.at((*BND_VERTEX_COLOR)[nodeA]);
      std::vector<unsigned int> edgeIdx_B = C2E_indices.at((*BND_VERTEX_COLOR)[nodeB]);
      std::find_if(edgeIdx_A.begin(), edgeIdx_A.end(), [=](unsigned int e){
        if(std::find(edgeIdx_B.begin(), edgeIdx_B.end(), e) != edgeIdx_B.end())
        {
          BND_CURVE_COLOR->set(newNodeId, e);
          return true;
        }
      });

      std::vector<unsigned int> surfaceIdx_A = C2S_indices.at((*BND_VERTEX_COLOR)[nodeA]);
      std::vector<unsigned int> surfaceIdx_B = C2S_indices.at((*BND_VERTEX_COLOR)[nodeB]);
      std::find_if(surfaceIdx_A.begin(), surfaceIdx_A.end(), [=](unsigned int s){
        if(std::find(surfaceIdx_B.begin(), surfaceIdx_B.end(), s) != surfaceIdx_B.end())
        {
          BND_SURFACE_COLOR->set(newNodeId, s);
          return true;
        }
      });*/
    }
    else if(newNodeDim == SimplexMesh::topo::RIDGE)
    {
      if(nodeDim == directNodeDim)
      {
        if(nodeLabel != directNodeLabel)
        {
          //node sur la surface or the volume...
          std::pair<unsigned int, unsigned int> p0 = E2S_indices.at(nodeLabel);
          std::pair<unsigned int, unsigned int> p1 = E2S_indices.at(directNodeLabel);
          bool flag = false;
          std::vector<unsigned int> p0Vec{p0.first, p0.second};
          std::vector<unsigned int> p1Vec{p1.first, p1.second};

          for(auto const surface0 : p0Vec)
          {
            for(auto const surface1 : p1Vec)
            {
              if(surface0 == surface1)
              {
                BND_SURFACE_COLOR->set(newNodeId, surface0);
                flag = !flag;
                break;
              }
            }
            if(flag){
              break;
            } // else it's on the volume
          }
        }
        else
        {
          BND_CURVE_COLOR->set(newNodeId, nodeLabel);
        }
      }
      else if(nodeDim > directNodeDim)
      {
        std::vector<unsigned int> edgeIdx = C2E_indices.at((*BND_VERTEX_COLOR)[nodeB]);
        if(std::find(edgeIdx.begin(), edgeIdx.end(), (*BND_CURVE_COLOR)[nodeA]) != edgeIdx.end())
        {
          BND_CURVE_COLOR->set(newNodeId, nodeLabel);
        }
        else
        {
          std::vector<unsigned int> surfaceIdx0 = C2S_indices.at((*BND_VERTEX_COLOR)[nodeB]);
          std::pair<unsigned int, unsigned int> surfaceIdx1 = E2S_indices.at((*BND_CURVE_COLOR)[nodeA]);
          for(auto const s : surfaceIdx0)
          {
            if(s == surfaceIdx1.first || s == surfaceIdx1.second)
            {
              BND_SURFACE_COLOR->set(newNodeId, s);
              break;
            }
          }
        }
      }
      else
      {
        std::vector<unsigned int> edgeIdx = C2E_indices.at((*BND_VERTEX_COLOR)[nodeA]);
        if(std::find(edgeIdx.begin(), edgeIdx.end(), (*BND_CURVE_COLOR)[nodeB]) != edgeIdx.end())
        {
          BND_CURVE_COLOR->set(newNodeId, directNodeLabel);
        }
        else
        {
          std::vector<unsigned int> surfaceIdx0 = C2S_indices.at((*BND_VERTEX_COLOR)[nodeA]);
          std::pair<unsigned int, unsigned int> surfaceIdx1 = E2S_indices.at((*BND_CURVE_COLOR)[nodeB]);
          for(auto const s : surfaceIdx0)
          {
            if(s == surfaceIdx1.first || s == surfaceIdx1.second)
            {
              BND_SURFACE_COLOR->set(newNodeId, s);
              break;
            }
          }
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
      else if(nodeDim > directNodeDim)
      {
        //on verifie si la surface est bien relié au ridge
        if(directNodeDim == SimplexMesh::topo::RIDGE)
        {
          const std::map<unsigned int, std::pair<unsigned int, unsigned int>>& ridgeEdge = m_simplexMesh->getEdgeTianglesIndices();
          std::pair<unsigned int, unsigned int> p1 = ridgeEdge.at(directNodeLabel);
          if(p1.first == nodeLabel || p1.second == nodeLabel)
          {
            BND_SURFACE_COLOR->set(newNodeId, nodeLabel);
          }
        }
        else if(directNodeDim == SimplexMesh::topo::CORNER)
        {
          std::vector<unsigned int> surfaceIdx = C2S_indices.at((*BND_VERTEX_COLOR)[nodeB]);
          if(std::find(surfaceIdx.begin(), surfaceIdx.end(), (*BND_SURFACE_COLOR)[nodeA]) != surfaceIdx.end())
          {
            BND_SURFACE_COLOR->set(newNodeId, nodeLabel);
          } //else it's on the volume
        }
      }
      else // nodeDim < directNodeDim
      {
        if(nodeDim == SimplexMesh::topo::RIDGE)
        {
          const std::map<unsigned int, std::pair<unsigned int, unsigned int>>& ridgeEdge = m_simplexMesh->getEdgeTianglesIndices();
          std::pair<unsigned int, unsigned int> p1 = ridgeEdge.at(nodeLabel);
          if(p1.first == directNodeLabel || p1.second == directNodeLabel)
          {
            BND_SURFACE_COLOR->set(newNodeId, directNodeLabel);
          }
        }
        else if(nodeDim == SimplexMesh::topo::CORNER)
        {
          std::vector<unsigned int> surfaceIdx = C2S_indices.at((*BND_VERTEX_COLOR)[nodeA]);
          if(std::find(surfaceIdx.begin(), surfaceIdx.end(), (*BND_SURFACE_COLOR)[nodeB]) != surfaceIdx.end())
          {
            BND_SURFACE_COLOR->set(newNodeId, directNodeLabel);
          }
        }
      }
    }
    //config data for the point insertion algorithm
    bool status = false;
    const gmds::BitVector markedNodes{};
    std::vector<TSimplexID> deletedSimplex{};
    std::vector<TInt> deletedNodes{};
    const std::multimap<TInt, TInt>& facesAlreadyBuilt{};
    std::vector<TSimplexID> cellsCreated{};
    PointInsertion pi(m_simplexMesh, SimplicesNode(m_simplexMesh, newNodeId), criterionRAIS, status, shell, markedNodes, deletedNodes, facesAlreadyBuilt, cellsCreated);

    return status;
  }
}
/*----------------------------------------------------------------------------*/
bool MetricAdaptation::computeEdgeRemove(const TInt nodeA, const TInt nodeB) const
{
  const gmds::BitVector& meshNode = m_simplexMesh->getBitVectorNodes();
  if(meshNode[nodeA] != 0 && meshNode[nodeB] != 0)
  {
    SimplicesNode sNodeA(m_simplexMesh, nodeA);
    SimplicesNode sNodeB(m_simplexMesh, nodeB);
    std::vector<TSimplexID> shell = sNodeA.shell(sNodeB);

    if(shell.size() == 0){return false;}
    bool flag = false;
    flag = m_simplexMesh->edgeRemove(nodeA, nodeB);
    if((nodeA == 1 && nodeB == 21) || (nodeB == 1 && nodeA == 21))
    {
      std::cout << sNodeA << std::endl;
      std::cout << sNodeB << std::endl;
      std::cout << "flag -> " << flag << std::endl;
    }
    return flag;
  }
}
/*----------------------------------------------------------------------------*/
unsigned int MetricAdaptation::computeSurfaceEdgeSwap()
{
  Variable<int>* BND_VERTEX_COLOR   = nullptr;
  Variable<int>* BND_CURVE_COLOR    = nullptr;
  Variable<int>* BND_SURFACE_COLOR  = nullptr;
  Variable<Eigen::Matrix3d>* metric = nullptr;
  try{
    BND_VERTEX_COLOR  = m_simplexMesh->getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
    BND_CURVE_COLOR   = m_simplexMesh->getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
    BND_SURFACE_COLOR = m_simplexMesh->getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
    metric            = m_simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");

  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }
  const gmds::BitVector& tetIDS = m_simplexMesh->getBitVectorTet();

  CriterionRAIS criterionRAIS(new VolumeCriterion());
  const gmds::BitVector& meshNode = m_simplexMesh->getBitVectorNodes();
  unsigned int computeSurfaceEdgeSwap = 0;

  const TInt border = std::numeric_limits<TInt>::min();
  unsigned int sizeFace = 4;

  buildEdgesMap();
  double halfInter = 3.0 * sqrt(2.0) * 0.25; // it represents the half valueof the interval [sqrt(2)*0.5 ; sqrt(2)]
  for(auto const edge : m_edgesMap)
  {
    TInt node0 = edge.first;
    TInt node1 = edge.second;
    double metricLenght0_1 = 0.0;
    //std::cout << "nodes -> " << node0 << " | " << node1 << " --> " << meshNode[node0] << " | " << meshNode[node1] <<  std::endl;
    if(meshNode[node0] != 0 && meshNode[node1] != 0)
    {
      unsigned int node0Dim = ((*BND_VERTEX_COLOR)[node0] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[node0] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[node0] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
      unsigned int node1Dim = ((*BND_VERTEX_COLOR)[node1] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[node1] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[node1] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
      unsigned int label0   = (node0Dim == SimplexMesh::topo::CORNER)?(*BND_VERTEX_COLOR)[node0]:(node0Dim == SimplexMesh::topo::RIDGE)?(*BND_CURVE_COLOR)[node0]:(node0Dim == SimplexMesh::topo::SURFACE)?(*BND_SURFACE_COLOR)[node0]:0;
      unsigned int label1   = (node1Dim == SimplexMesh::topo::CORNER)?(*BND_VERTEX_COLOR)[node1]:(node1Dim == SimplexMesh::topo::RIDGE)?(*BND_CURVE_COLOR)[node1]:(node1Dim == SimplexMesh::topo::SURFACE)?(*BND_SURFACE_COLOR)[node1]:0;

      if(node0Dim == SimplexMesh::topo::SURFACE && node1Dim == SimplexMesh::topo::SURFACE)
      {
        if(label0 == label1)
        {
          SimplicesNode sNode0(m_simplexMesh, node0);
          SimplicesNode sNode1(m_simplexMesh, node1);
          std::vector<TSimplexID> shell = sNode0.shell(sNode1);
          if(shell.size() == 0){continue;}

          math::Point nodeCoord0 = sNode0.getCoords();
          math::Point nodeCoord1 = sNode1.getCoords();
          Metric<Eigen::Matrix3d> M0 = Metric<Eigen::Matrix3d>((*metric)[node0]);
          Metric<Eigen::Matrix3d> M1 = Metric<Eigen::Matrix3d>((*metric)[node1]);
          metricLenght0_1        = M0.metricDist(nodeCoord0, nodeCoord1, M1);

          if(metricLenght0_1 < sqrt(2.0)*0.5 || metricLenght0_1 > sqrt(2.0))
          {
            std::vector<TInt> nodesSurface{};
            std::vector<TSimplexID> cavity = sNode0.shell(sNode1);
            for(auto const simplex : cavity)
            {
              if(simplex < 0)
              {
                const SimplicesTriangle triangle(m_simplexMesh, -simplex);
                std::vector<TInt> nodesTriangle = triangle.getNodes();
                std::copy(nodesTriangle.begin(), nodesTriangle.end(), std::back_inserter(nodesSurface));
              }
            }

            std::set<TInt> sameLabelNode{};
            for(auto const node : nodesSurface)
            {
              if(node != node0 && node != node1 && (*BND_SURFACE_COLOR)[node] == label0)
              {
                sameLabelNode.insert(node);
              }
            }

            if(sameLabelNode.size() == 2)
            {
              TInt n00 = *(sameLabelNode.begin());
              TInt n11 = *(--sameLabelNode.end());
              SimplicesNode sNode00(m_simplexMesh, n00);
              SimplicesNode sNode11(m_simplexMesh, n11);

              math::Point nodeCoord00 = sNode00.getCoords();
              math::Point nodeCoord11 = sNode11.getCoords();
              Metric<Eigen::Matrix3d> M0 = Metric<Eigen::Matrix3d>((*metric)[n00]);
              Metric<Eigen::Matrix3d> M1 = Metric<Eigen::Matrix3d>((*metric)[n11]);
              double metricLenght00_11   = M0.metricDist(nodeCoord00, nodeCoord11, M1);

              if(std::abs(metricLenght00_11 - halfInter) < std::abs(metricLenght0_1 - halfInter))
              {
                const gmds::BitVector markedNodes{};
                std::vector<TSimplexID> deletedSimplex{};
                std::vector<TInt> deletedNodes{};
                const std::multimap<TInt, TInt>& facesAlreadyBuilt{};
                std::vector<TSimplexID> cellsCreated{};
                bool status = false;
                PointInsertion pi(m_simplexMesh, sNode00, criterionRAIS, status, cavity, markedNodes, deletedNodes, facesAlreadyBuilt, cellsCreated);
                if(status)
                {
                  computeSurfaceEdgeSwap++;
                  /*gmds::ISimplexMeshIOService ioServiceMesh(m_simplexMesh);
                  gmds::VTKWriter vtkWriterCS(&ioServiceMesh);
                  vtkWriterCS.setCellOptions(gmds::N|gmds::R|gmds::F);
                  vtkWriterCS.setDataOptions(gmds::N|gmds::R|gmds::F);
                  vtkWriterCS.write("ES_" + std::to_string(computeSurfaceEdgeSwap) + ".vtk");*/
                }
              }
            }
          }
        }
      }
    }
  }
  return computeSurfaceEdgeSwap;
}
/*----------------------------------------------------------------------------*/
unsigned int MetricAdaptation::computeFaceSwap()
{
  static unsigned int cptEXE = 0;
  std::cout << "cptEXE -> " << cptEXE << std::endl;
  gmds::ISimplexMeshIOService ioServiceMesh(m_simplexMesh);
  gmds::VTKWriter vtkWriterCS(&ioServiceMesh);
  vtkWriterCS.setCellOptions(gmds::N|gmds::R|gmds::F);
  vtkWriterCS.setDataOptions(gmds::N|gmds::R|gmds::F);
  vtkWriterCS.write("cptEXE_" + std::to_string(cptEXE) + ".vtk");

  Variable<int>* BND_VERTEX_COLOR   = nullptr;
  Variable<int>* BND_CURVE_COLOR    = nullptr;
  Variable<int>* BND_SURFACE_COLOR  = nullptr;
  Variable<Eigen::Matrix3d>* metric = nullptr;
  try{
    BND_VERTEX_COLOR  = m_simplexMesh->getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
    BND_CURVE_COLOR   = m_simplexMesh->getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
    BND_SURFACE_COLOR = m_simplexMesh->getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
    metric            = m_simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");

  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }
  const gmds::BitVector& tetIDS = m_simplexMesh->getBitVectorTet();
  //gmds::ISimplexMeshIOService ioServiceMesh(m_simplexMesh);

  CriterionRAIS criterionRAIS(new VolumeCriterion());
  const gmds::BitVector& meshNode = m_simplexMesh->getBitVectorNodes();
  unsigned int computeFaceSwap = 0;

  const TInt border = std::numeric_limits<TInt>::min();
  unsigned int sizeFace = 4;

  for(TSimplexID tetId = 0 ; tetId < tetIDS.capacity() ; tetId++)
  {
    if(tetIDS[tetId] != 0)
    {
      bool status = false;
      for(unsigned int face = 0 ; face < sizeFace ; face++)
      {
        SimplicesCell cell0(m_simplexMesh, tetId);
        const TInt node_0 = cell0.getNodes()[face];
        TSimplexID oppositeTet = cell0.oppositeTetraIdx(face);
        if((*BND_VERTEX_COLOR)[node_0] != 0)
        {
          continue;
        }
        if(oppositeTet < 0)//it is a triangle
        {
          continue;
        }
        else
        {
          SimplicesCell cell1(m_simplexMesh, oppositeTet);
          //check if the swapping is possible
          std::vector<TInt> visibleFaces = cell1.visibleFaces(SimplicesNode(m_simplexMesh, node_0).getCoords());
          std::vector<TInt> nodesCell    = cell1.getNodes();
          if(visibleFaces.size() == 3)
          {
            std::vector<TInt> nodesInSurface{};
            std::vector<TSimplexID> cavity{};
            if((*BND_SURFACE_COLOR)[node_0] != 0)
            {
              for(auto const node : visibleFaces)
              {
                if((*BND_SURFACE_COLOR)[nodesCell[node]] != 0 && (*BND_SURFACE_COLOR)[nodesCell[node]] == (*BND_SURFACE_COLOR)[node_0])
                {
                  nodesInSurface.push_back(nodesCell[node]);
                }
              }

              if(nodesInSurface.size() == 2)
              {
                if((*BND_SURFACE_COLOR)[nodesInSurface.front()] == (*BND_SURFACE_COLOR)[nodesInSurface.back()])
                {
                  cavity = SimplicesNode(m_simplexMesh, nodesInSurface.front()).shell(SimplicesNode(m_simplexMesh, nodesInSurface.back()));
                }
                else
                {
                  continue;
                }
              }
            }
            else
            {
              continue;
            }
            if(cavity.size() == 0)
            {
              ///////
              continue;
              cavity.push_back(tetId);
              cavity.push_back(oppositeTet);
            }
            //std::cout << "cavity.size() -> " << cavity.size() << std::endl;
            double qualityCell0 = m_simplexMesh->computeQualityElement(tetId);
            double qualityCell1 = m_simplexMesh->computeQualityElement(oppositeTet);
            double worstQuality = std::min(qualityCell0, qualityCell1);

            std::vector<TInt> nodesFace = cell0.getOrderedFace(face);


            if(nodesFace.size() == 3)
            {
              std::vector<TInt> node_1 = cell1.getOtherNodeInSimplex(nodesFace);
              if(node_1.size() == 1)
              {
                std::vector<TInt> newTet_0{nodesFace[0], nodesFace[1], node_1.front(), node_0};
                std::vector<TInt> newTet_1{nodesFace[1], nodesFace[2], node_1.front(), node_0};
                std::vector<TInt> newTet_2{nodesFace[2], nodesFace[0], node_1.front(), node_0};

                double qualityElementtet_0 = m_simplexMesh->computeQualityElement(newTet_0[0], newTet_0[1], newTet_0[2], newTet_0[3]);
                double qualityElementtet_1 = m_simplexMesh->computeQualityElement(newTet_1[0], newTet_1[1], newTet_1[2], newTet_1[3]);
                double qualityElementtet_2 = m_simplexMesh->computeQualityElement(newTet_2[0], newTet_2[1], newTet_2[2], newTet_2[3]);

                double worstNewElement = std::min(qualityElementtet_0, std::min(qualityElementtet_1, qualityElementtet_2));

                if(worstNewElement > worstQuality)
                {
                  const gmds::BitVector markedNodes{};
                  std::vector<TSimplexID> deletedSimplex{};
                  std::vector<TInt> deletedNodes{};
                  const std::multimap<TInt, TInt>& facesAlreadyBuilt{};
                  std::vector<TSimplexID> cellsCreated{};
                  PointInsertion pi(m_simplexMesh, SimplicesNode(m_simplexMesh, node_0), criterionRAIS, status, cavity, markedNodes, deletedNodes, facesAlreadyBuilt, cellsCreated);
                  if(status)
                  {
                    computeFaceSwap++;
                    gmds::ISimplexMeshIOService ioServiceMesh(m_simplexMesh);
                    gmds::VTKWriter vtkWriterCS(&ioServiceMesh);
                    vtkWriterCS.setCellOptions(gmds::N|gmds::R|gmds::F);
                    vtkWriterCS.setDataOptions(gmds::N|gmds::R|gmds::F);
                    vtkWriterCS.write("FS_" + std::to_string(cptEXE) + "_" + std::to_string(computeFaceSwap) + ".vtk");
                    break;
                  }
                }
              }
              else
              {
                throw gmds::GMDSException("node_1.size() != 1");
              }
            }
            else
            {
              throw gmds::GMDSException("nodesFace.size() != 3");
            }
          }
        }
      }
    }
  }
  cptEXE++;
  return computeFaceSwap;
}
/*----------------------------------------------------------------------------*/
unsigned int MetricAdaptation::computePointSmoothing()
{
  static unsigned int cptTT = 0;
  unsigned int i = 0;
  Variable<int>* BND_VERTEX_COLOR   = nullptr;
  Variable<int>* BND_CURVE_COLOR    = nullptr;
  Variable<int>* BND_SURFACE_COLOR  = nullptr;
  Variable<Eigen::Matrix3d>* metric = nullptr;
  try{
    BND_VERTEX_COLOR  = m_simplexMesh->getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
    BND_CURVE_COLOR   = m_simplexMesh->getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
    BND_SURFACE_COLOR = m_simplexMesh->getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
    metric            = m_simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");

  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  gmds::ISimplexMeshIOService ioServiceMesh(m_simplexMesh);
  CriterionRAIS criterionRAIS(new VolumeCriterion());
  const gmds::BitVector& meshNode = m_simplexMesh->getBitVectorNodes();
  unsigned int computePointSmoothing = 0;
  double epsilon = 0.01;

  for(TInt nodeId = 0 ; nodeId < meshNode.capacity() ; nodeId++)
  {
    if(meshNode[nodeId] != 0)
    {
      unsigned int nodeIdDim   = ((*BND_VERTEX_COLOR)[nodeId] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[nodeId] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[nodeId] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
      unsigned int nodeIdLabel = (nodeIdDim == SimplexMesh::topo::CORNER)?(*BND_VERTEX_COLOR)[nodeId]:(nodeIdDim == SimplexMesh::topo::RIDGE)?(*BND_CURVE_COLOR)[nodeId]:(nodeIdDim == SimplexMesh::topo::SURFACE)?(*BND_SURFACE_COLOR)[nodeId]:0;

      if(!(nodeIdDim == SimplexMesh::topo::SURFACE))
      {
        continue;
      }

      Metric<Eigen::Matrix3d> M0 = Metric<Eigen::Matrix3d>((*metric)[nodeId]);
      const SimplicesNode   node = SimplicesNode(m_simplexMesh, nodeId);
      const math::Point        p = node.getCoords();
      const std::vector<TInt> nodes = node.neighborSubSurfaceNodes();
      std::vector<TSimplexID> ball = node.ballOf();
      if(ball.size() == 0){continue;}
      math::Point p_opti;

      for(auto const n : nodes)
      {
        if(meshNode[n] != 0)
        {
          const SimplicesNode sNode   = SimplicesNode(m_simplexMesh, n);
          const math::Point pj        = sNode.getCoords();
          Metric<Eigen::Matrix3d> M1  = Metric<Eigen::Matrix3d>((*metric)[n]);
          double metricLenght         = M0.metricDist(p, pj, M1);

          if(metricLenght != 0.0)
          {
            math::Vector3d vec = pj - p;
            p_opti = p_opti + (p + (vec / metricLenght));
          }
          else
          {
            throw gmds::GMDSException("metricLenght != 0.0");
          }
        }
      }

      math::Point new_P = p_opti;
      new_P = math::Point(new_P.X() / (double)nodes.size(), new_P.Y() / (double)nodes.size(), new_P.Z() / (double)nodes.size());

      //test pour calculer si new_P est étoilable par rapport a la boule volumique ou surfacique de nodeId
      bool flag = false;
      for(auto const simplex : ball)
      {
        if(simplex >= 0/* && nodeIdDim == SimplexMesh::topo::VOLUME*/)
        {
          unsigned int localNode = SimplicesCell(m_simplexMesh, simplex).getLocalNode(nodeId);
          if(criterionRAIS.execute(m_simplexMesh, simplex, localNode, new_P))
          {
            flag = true;
            break;
          }
        }
      }
      if(flag)
      {
        continue;
      }

      m_simplexMesh->moveNodeCoord(nodeId, new_P);
      m_simplexMesh->setAnalyticMetric(nodeId);
      computePointSmoothing++;
      /*gmds::ISimplexMeshIOService ioServiceMesh(m_simplexMesh);
      gmds::VTKWriter vtkWriterCS(&ioServiceMesh);
      vtkWriterCS.setCellOptions(gmds::N|gmds::R|gmds::F);
      vtkWriterCS.setDataOptions(gmds::N|gmds::R|gmds::F);
      vtkWriterCS.write("PS_" + std::to_string(cptTT) + "_" + std::to_string(i) + ".vtk");
      i++;*/
    }
  }
  //cptTT++;
  return computePointSmoothing;
}
/*----------------------------------------------------------------------------*/
void MetricAdaptation::execute()
{
  Variable<int>* BND_VERTEX_COLOR   = nullptr;
  Variable<int>* BND_SURFACE_COLOR   = nullptr;
  Variable<int>* BND_CURVE_COLOR   = nullptr;

  Variable<Eigen::Matrix3d>* metric = nullptr;
  try{
    metric = m_simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
    BND_VERTEX_COLOR  = m_simplexMesh->getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
    BND_SURFACE_COLOR  = m_simplexMesh->getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
    BND_CURVE_COLOR  = m_simplexMesh->getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  unsigned int cpt = 0;
  unsigned int maxIterationAdaptation = 50;
  CriterionRAIS criterionRAIS(new VolumeCriterion());
  gmds::ISimplexMeshIOService ioServiceMesh(m_simplexMesh);
  const gmds::BitVector& meshNode = m_simplexMesh->getBitVectorNodes();

  //ADAPTATION
  std::cout << "ADAPTATION  START" << std::endl;
  static unsigned int cptTest = 0;
  for(unsigned int iter = 0 ; iter < maxIterationAdaptation ; iter++)
  {
    std::cout << std::endl;
    std::cout << " cpt -> " << cpt << std::endl;
    unsigned int cptSmoothing = 0;
    unsigned int cptFaceSwap = 0;
    unsigned int cptEdgeSwap = 0;
    double meanEdge = 0.0; double minEdge = 0.0; double maxEdge = 0.0;

    //Metric adaptation BY adding point or removing them
    std::cout << "SLICING/EDGE REMOVE START" << std::endl;
    buildEdgesMap();
    unsigned int cptt = 0;
    for(auto const edge : m_edgesMap)
    {
      TInt nodeA = edge.first;
      TInt nodeB = edge.second;
      if(meshNode[nodeA] != 0 && meshNode[nodeB] != 0)
      {
        SimplicesNode sNodeA(m_simplexMesh, nodeA);
        SimplicesNode sNodeB(m_simplexMesh, nodeB);

        std::vector<TSimplexID> shell = sNodeA.shell(sNodeB);
        if(shell.size() == 0){continue;}

        math::Point nodeCoord0 = sNodeA.getCoords();
        math::Point nodeCoord1 = sNodeB.getCoords();
        Metric<Eigen::Matrix3d> M0 = Metric<Eigen::Matrix3d>((*metric)[nodeA]);
        Metric<Eigen::Matrix3d> M1 = Metric<Eigen::Matrix3d>((*metric)[nodeB]);
        double metricLenght        = M0.metricDist(nodeCoord0, nodeCoord1, M1);

        if(metricLenght > sqrt(2))
        {
          computeSlicing(nodeA, nodeB);
        }
        else if(metricLenght < 0.5 * sqrt(2))
        {

          computeEdgeRemove(nodeA, nodeB);
        }
        /*gmds::ISimplexMeshIOService ioServiceMesh(m_simplexMesh);
        gmds::VTKWriter vtkWriterCS(&ioServiceMesh);
        vtkWriterCS.setCellOptions(gmds::N|gmds::R|gmds::F);
        vtkWriterCS.setDataOptions(gmds::N|gmds::R|gmds::F);
        vtkWriterCS.write("TEST_" + std::to_string(cptTest) + "_" + std::to_string(cptt) + ".vtk");
        cptt++;*/
      }
    }
    //cptTest++;

    gmds::VTKWriter vtkWriterCS(&ioServiceMesh);
    vtkWriterCS.setCellOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterCS.setDataOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterCS.write("ComputeSlicing_EdgeRemove_" + std::to_string(cpt) + ".vtk");

    //BOUGE DE POINTS
    std::cout << "POINT SMOOTHING START" << std::endl;
    cptSmoothing = computePointSmoothing();
    std::cout << "    cptSmoothing -> " << cptSmoothing << std::endl;
    gmds::VTKWriter vtkWriterPS(&ioServiceMesh);
    vtkWriterPS.setCellOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterPS.setDataOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterPS.write("ComputePoint_Smoothing_" + std::to_string(cpt) + ".vtk");
    std::cout << "MESH VALIDITY CHECK" << std::endl;


    /*std::cout << "  SWAPPING EDGE START" << std::endl;
    cptEdgeSwap = computeSurfaceEdgeSwap();
    std::cout << "    cptEdgeSwap -> " << cptEdgeSwap << std::endl;
    gmds::VTKWriter vtkWriterMeshSwapping(&ioServiceMesh);
    vtkWriterMeshSwapping.setCellOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterMeshSwapping.setDataOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterMeshSwapping.write("ComputeSurfaceSwap_" + std::to_string(cpt) + ".vtk");*/

    /*//SWAPPING FACE
    std::cout << "  SWAPPING FACE START" << std::endl;
    //do{
      cptSwcptFaceSwapap = computeFaceSwap();
      std::cout << "    cptFaceSwap -> " << cptFaceSwap << std::endl;
    //}
    //while(cptSwap != 0);
    gmds::VTKWriter vtkWriterMeshSwapping(&ioServiceMesh);
    vtkWriterMeshSwapping.setCellOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterMeshSwapping.setDataOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterMeshSwapping.write("ComputeFaceSwap_" + std::to_string(cpt) + ".vtk");*/
    cpt++;
    m_simplexMesh->getEdgeSizeInfowithMetric(meanEdge, minEdge, maxEdge);
  }
}
/*----------------------------------------------------------------------------*/
void MetricAdaptation::executeCustomMethod()
{
  CriterionRAIS criterionRAIS(new VolumeCriterion());

  Variable<int>* BND_VERTEX_COLOR   = nullptr;
  Variable<int>* BND_SURFACE_COLOR   = nullptr;
  Variable<int>* BND_CURVE_COLOR   = nullptr;
  Variable<Eigen::Matrix3d>* metric = nullptr;

  try{
    metric = m_simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");
    BND_VERTEX_COLOR  = m_simplexMesh->getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
    BND_SURFACE_COLOR  = m_simplexMesh->getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
    BND_CURVE_COLOR  = m_simplexMesh->getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  const math::Vector3d v0 = math::Vector3d(-sqrt(2.0) / 2.0 , std::sqrt(2.0) / 2.0 , 0.0);
  const gmds::BitVector& meshNode = m_simplexMesh->getBitVectorNodes();


  static unsigned int cptTEST = 0;
  unsigned int maxIterationAdaptation = 50;
  for(unsigned int iter = 0 ; iter < maxIterationAdaptation ; iter++)
  {
    std::cout << "iteration -> " << cptTEST << std::endl;
    unsigned int cpt = 0;
    buildEdgesMap();
    for(auto const edge : m_edgesMap)
    {
      TInt nodeA = edge.first;
      TInt nodeB = edge.second;
      if(meshNode[nodeA] != 0 && meshNode[nodeB] != 0)
      {
        const SimplicesNode sNodeA(m_simplexMesh, nodeA);
        const SimplicesNode sNodeB(m_simplexMesh, nodeB);
        const math::Point ptA = sNodeA.getCoords();

        /*std::cout << "nodeA -> " << nodeA << std::endl;
        std::cout << "nodeB -> " << nodeB << std::endl;*/
        if((*BND_SURFACE_COLOR)[nodeA] == 5 && (*BND_SURFACE_COLOR)[nodeB] == 5)
        {
          const double m = metric->value(nodeB)(0, 0);

          unsigned int nodeDim = ((*BND_VERTEX_COLOR)[nodeB] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[nodeB] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[nodeB] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
          unsigned int nodeLabel = ((*BND_VERTEX_COLOR)[nodeB] != 0)?(*BND_VERTEX_COLOR)[nodeA]:((*BND_CURVE_COLOR)[nodeB] != 0)?(*BND_CURVE_COLOR)[nodeA]:((*BND_SURFACE_COLOR)[nodeB] != 0)?(*BND_SURFACE_COLOR)[nodeB]:0;

          double borneMin = sqrt(1.0 / (2.0 * m * (std::pow(v0.Y() / v0.X(), 2) + 1))) + ptA.Y();
          double borneMax = sqrt(2.0 / ( m * (std::pow(v0.Y() / v0.X(), 2) + 1))) + ptA.Y();

          double newPosition_yB = (borneMin + borneMax) / 2.0;
          double newPosition_xB = ptA.X() + (ptA.Y() - newPosition_yB)*(v0.Y() / v0.X());
          const Point pt = Point(newPosition_xB, newPosition_yB, 0.0);
          const Vector3d p = Vector3d(pt.X() - ptA.X(), pt.Y() - ptA.Y(), pt.Z() - ptA.Z()) ;
          //std::cout << "p length -> " << p.norm() << std::endl;
          bool alreadyAdd = false;
          std::vector<TSimplexID> tetraContenaingPt{};
          TInt newNodeId = m_simplexMesh->addNodeAndcheck(pt, tetraContenaingPt, alreadyAdd);
          if(!alreadyAdd)
          {
            m_simplexMesh->setAnalyticMetric(newNodeId);
            if(nodeDim == SimplexMesh::topo::CORNER)
            {
              BND_VERTEX_COLOR->set(newNodeId, nodeLabel);
            }
            else if(nodeDim == SimplexMesh::topo::RIDGE)
            {
              BND_CURVE_COLOR->set(newNodeId, nodeLabel);
            }
            else if(nodeDim == SimplexMesh::topo::SURFACE)
            {
              BND_SURFACE_COLOR->set(newNodeId, nodeLabel);
            }
          }

          bool status = false;
          const gmds::BitVector markedNodes{};
          std::vector<TSimplexID> deletedSimplex{};
          std::vector<TInt> deletedNodes{};
          const std::multimap<TInt, TInt>& facesAlreadyBuilt{};
          std::vector<TSimplexID> cellsCreated{};
          //std::vector<TSimplexID> cavity = sNodeB.ballOf();
          std::vector<TSimplexID> cavity = sNodeB.shell(sNodeA);
          gmds::BitVector nodesAdded(m_simplexMesh->nodesCapacity());
          //DelaunayPointInsertion DI(m_simplexMesh, SimplicesNode(m_simplexMesh, newNodeId), criterionRAIS, cavity, status, nodesAdded, deletedSimplex, facesAlreadyBuilt);
          PointInsertion pi(m_simplexMesh, SimplicesNode(m_simplexMesh, newNodeId), criterionRAIS, status, cavity, markedNodes, deletedNodes, facesAlreadyBuilt, cellsCreated);
          /*std::cout << "status -> " << status << std::endl;
          std::cout << std::endl;*/
          if(status)
          {
            gmds::ISimplexMeshIOService ioServiceMesh(m_simplexMesh);
            gmds::VTKWriter vtkWriterCS(&ioServiceMesh);
            vtkWriterCS.setCellOptions(gmds::N|gmds::R|gmds::F);
            vtkWriterCS.setDataOptions(gmds::N|gmds::R|gmds::F);
            vtkWriterCS.write("TEST2_Metric_Vector_"+ std::to_string(cptTEST)+ "_" + std::to_string(cpt) + ".vtk");
            cpt++;
          }
        }
      }
    }
    cptTEST++;
  }
}
