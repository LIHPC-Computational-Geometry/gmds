#include "gmds/hybridMeshAdapt/PointInsertion.h"
#include "gmds/hybridMeshAdapt/MetricAdaptation.h"
#include "gmds/hybridMeshAdapt/SimplexMesh.h"
#include "gmds/hybridMeshAdapt/SimplicesCell.h"
#include "gmds/hybridMeshAdapt/SimplicesNode.h"
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
  m_edgesMap.clear();
  //generate the edge structure in order to iterate over edge and not the mesh's node
  const gmds::BitVector& meshNode = m_simplexMesh->getBitVectorNodes();
  bool flag = false;
  for(TInt nodeId = 0 ; nodeId < meshNode.capacity() ; nodeId++)
  {
    flag = false;
    if(meshNode[nodeId] != 0)
    {
      SimplicesNode node(m_simplexMesh, nodeId);
      const std::vector<TInt> && directNodes = node.directNeighboorNodeId();
      for(auto dNode : directNodes)
      {
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
        if(flag)
        {
          break;
        }
        m_edgesMap.insert(std::pair<TInt, TInt>(std::min(dNode, nodeId), std::max(dNode, nodeId)));
      }
    }
  }

  /*for(auto const m : m_edgesMap){
    std::cout << m.first << " | " << m.second << std::endl;
  }*/
}
/*----------------------------------------------------------------------------*/
void MetricAdaptation::execute()
{
  Variable<int>* BND_VERTEX_COLOR ;
  Variable<int>* BND_CURVE_COLOR  ;
  Variable<int>* BND_SURFACE_COLOR;
  Variable<Eigen::Matrix3d>* metric;
  try{
    BND_VERTEX_COLOR  = m_simplexMesh->getVariable<int,SimplicesNode>("BND_VERTEX_COLOR");
    BND_CURVE_COLOR   = m_simplexMesh->getVariable<int,SimplicesNode>("BND_CURVE_COLOR");
    BND_SURFACE_COLOR = m_simplexMesh->getVariable<int,SimplicesNode>("BND_SURFACE_COLOR");
    metric            = m_simplexMesh->getVariable<Eigen::Matrix3d, SimplicesNode>("NODE_METRIC");

  }catch (gmds::GMDSException e)
  {
    throw gmds::GMDSException(e);
  }

  unsigned int maxIterationAdaptation = 1;
  CriterionRAIS criterionRAIS(new VolumeCriterion());
  const gmds::BitVector& meshNode = m_simplexMesh->getBitVectorNodes();
  unsigned int cpt = 0;
  for(unsigned int iter = 0 ; iter < maxIterationAdaptation ; iter++)
  {
    buildEdgesMap();
    //Metric adaptation
    for(auto const edge : m_edgesMap)
    {
      TInt nodeId = edge.first;
      TInt dNode  = edge.second;

      if(meshNode[nodeId] != 0 && meshNode[dNode])
      {
        SimplicesNode node0(m_simplexMesh, nodeId);
        SimplicesNode node1(m_simplexMesh, dNode);
        if(node0.shell(node1).size() == 0){continue;}

        math::Point nodeCoord0 = node0.getCoords();
        math::Point nodeCoord1 = node1.getCoords();
        Metric<Eigen::Matrix3d> M0 = Metric<Eigen::Matrix3d>((*metric)[nodeId]);
        Metric<Eigen::Matrix3d> M1 = Metric<Eigen::Matrix3d>((*metric)[dNode]);
        double metricLenght        = M0.metricDist(nodeCoord0, nodeCoord1, M1);

        //config data for the point insertion algorithm
        bool status = false;
        const gmds::BitVector markedNodes{};
        std::vector<TSimplexID> deletedSimplex{};
        std::vector<TInt> deletedNodes{};
        const std::multimap<TInt, TInt>& facesAlreadyBuilt{};
        std::vector<TSimplexID> cellsCreated{};

        unsigned int nodeDim = ((*BND_VERTEX_COLOR)[nodeId] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[nodeId] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[nodeId] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
        unsigned int nodeLabel = ((*BND_VERTEX_COLOR)[nodeId] != 0)?(*BND_VERTEX_COLOR)[nodeId]:((*BND_CURVE_COLOR)[nodeId] != 0)?(*BND_CURVE_COLOR)[nodeId]:((*BND_SURFACE_COLOR)[nodeId] != 0)?(*BND_SURFACE_COLOR)[nodeId]:0;
        unsigned int directNodeDim = ((*BND_VERTEX_COLOR)[dNode] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[dNode] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[dNode] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
        unsigned int directNodeLabel = ((*BND_VERTEX_COLOR)[dNode] != 0)?(*BND_VERTEX_COLOR)[dNode]:((*BND_CURVE_COLOR)[dNode] != 0)?(*BND_CURVE_COLOR)[dNode]:((*BND_SURFACE_COLOR)[dNode] != 0)?(*BND_SURFACE_COLOR)[dNode]:0;
        unsigned int newNodeDim = std::max(nodeDim, directNodeDim);


        if(metricLenght < 1.0 / sqrt(2.0))
        {
          gmds::BitVector nodeToRemove(meshNode.capacity());
          for(unsigned int node = 0 ; node < meshNode.capacity() ; node++)
          {
            if(meshNode[node] != 0)
            {
              nodeToRemove.assign(node);
            }
          }
          nodeToRemove.unselect(nodeId);
          unsigned int edgesRemoved = m_simplexMesh->edgesRemove(nodeToRemove, deletedNodes);
          if(deletedNodes.size() != 0)
          {
            m_simplexMesh->deleteNode(deletedNodes.front(), true);
          }
        }
        else if(metricLenght > sqrt(2.0))
        {
          bool alreadyAdd = false;
          std::vector<TSimplexID> tetraContenaingPt{};
          const Point pt = 0.5 * (nodeCoord0 + nodeCoord1);
          TInt newNodeId = m_simplexMesh->addNodeAndcheck(pt, tetraContenaingPt, alreadyAdd);

          metric->set(newNodeId, (*metric)[nodeId]);
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
                const std::map<unsigned int, std::pair<unsigned int, unsigned int>>& ridgeEdge = m_simplexMesh->getEdgeTianglesIndices();
                std::pair<unsigned int, unsigned int> p0 = ridgeEdge.at(nodeId);
                std::pair<unsigned int, unsigned int> p1 = ridgeEdge.at(dNode);
                bool flag = false;
                std::vector<unsigned int> p0Vec{p0.first, p0.second};
                std::vector<unsigned int> p1Vec{p1.first, p1.second};

                for(auto const surface0 : p0Vec)
                {
                  for(auto const surface1 : p1Vec)
                  {
                    if(surface0 == surface1)
                    {
                      BND_CURVE_COLOR->set(newNodeId, surface0);
                      flag = !flag;
                      break;
                    }
                  }
                  if(flag){
                    break;
                  }
                }

                if(!flag)
                {
                  throw gmds::GMDSException("no common surface");
                }
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
                //on verifie si la surface est bien relié au ridge
                if(directNodeDim == SimplexMesh::topo::RIDGE)
                {
                  const std::map<unsigned int, std::pair<unsigned int, unsigned int>>& ridgeEdge = m_simplexMesh->getEdgeTianglesIndices();
                  std::pair<unsigned int, unsigned int> p1 = ridgeEdge.at(directNodeLabel);
                  std::vector<unsigned int> p1Vec{p1.first, p1.second};
                  if(p1Vec.front() == nodeLabel || p1Vec.back() == nodeLabel)
                  {
                    BND_SURFACE_COLOR->set(newNodeId, nodeLabel);
                  }
                }
                else if(directNodeDim == SimplexMesh::topo::CORNER)
                {
                  //TODO when merging with the master (main) because on this branch there is no information about Corner label -> Surface label
                }
              }
              else
              {
                BND_SURFACE_COLOR->set(newNodeId, directNodeLabel);
              }
            }
          }
          const std::vector<TSimplexID>&& shell = node0.shell(node1);
          PointInsertion pi(m_simplexMesh, SimplicesNode(m_simplexMesh, newNodeId), criterionRAIS, status, shell, markedNodes, deletedNodes, facesAlreadyBuilt, cellsCreated);
        }
      }
    }

    gmds::ISimplexMeshIOService ioServiceMesh(m_simplexMesh);
    gmds::VTKWriter vtkWriterMesh(&ioServiceMesh);
    vtkWriterMesh.setCellOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterMesh.setDataOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterMesh.write("ADAPTAION_FIRST_STEP.vtk");

    //continue;
    //Swapping face
    std::cout << "SWAPPING FACE" << std::endl;
    const gmds::BitVector& tetBitVec = m_simplexMesh->getBitVectorTet();
    unsigned int sizeFace = 4;
    unsigned int cpt = 0;
    for(TSimplexID tetId = 0 ; tetId < meshNode.capacity() ; tetId++)
    {
      bool flag = false;
      if(tetBitVec[tetId] != 0)
      {
        for(unsigned int face = 0 ; face < sizeFace ; face++)
        {
          if(flag)
          {
            break;
          }
          SimplicesCell cell0(m_simplexMesh, tetId);
          std::vector<TInt> nodesFace = cell0.getOrderedFace(face);
          unsigned int nodeDim0 = ((*BND_VERTEX_COLOR)[nodesFace.front()] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[nodesFace.front()] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[nodesFace.front()] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
          unsigned int nodeDim1 = ((*BND_VERTEX_COLOR)[nodesFace[1]] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[nodesFace[1]] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[nodesFace[1]] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
          unsigned int nodeDim2 = ((*BND_VERTEX_COLOR)[nodesFace.back()] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[nodesFace.back()] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[nodesFace.back()] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;

          TSimplexID oppositeTet = cell0.oppositeTetraIdx(face);
          if(oppositeTet < 0)//it is a triangle
          {
            continue;
          }
          SimplicesCell cell1(m_simplexMesh, oppositeTet);
          std::vector<TInt> node_0 = cell0.getOtherNodeInSimplex(nodesFace);
          std::vector<TInt> node_1 = cell1.getOtherNodeInSimplex(nodesFace);
          if(node_0.size() != 1 || node_1.size() != 1)
          {
            throw gmds::GMDSException("node_0.size() != 1 || node_1.size() != 1");
          }
          else
          {
            std::vector<TSimplexID> cavity{tetId, oppositeTet};
            double qualityCell0 = m_simplexMesh->computeQualityElement(tetId);
            double qualityCell1 = m_simplexMesh->computeQualityElement(oppositeTet);
            double worstQuality = std::min(qualityCell0, qualityCell1);

            std::vector<TInt> nodesToInsert{node_0.front(), node_1.front()};
            for(auto const node : nodesToInsert)
            {
              //test to check if one of the inserted node formed a better quality mesh
              bool status = false;
              const gmds::BitVector markedNodes{};
              std::vector<TSimplexID> deletedSimplex{};
              std::vector<TInt> deletedNodes{};
              const std::multimap<TInt, TInt>& facesAlreadyBuilt{};
              std::vector<TSimplexID> cellsCreated{};
              SimplexMesh simplexMeshCopy = *m_simplexMesh;
              PointInsertion pi(&simplexMeshCopy, SimplicesNode(&simplexMeshCopy, node), criterionRAIS, status, cavity, markedNodes, deletedNodes, facesAlreadyBuilt, cellsCreated);

              if(status && cellsCreated.size() != 0)
              {
                double qual = std::numeric_limits<double>::max();
                for(auto const simplex : cellsCreated)
                {
                  double quality = simplexMeshCopy.computeQualityElement(simplex);
                  if(quality < qual)
                  {
                    qual = quality;
                  }
                }

                if(qual > worstQuality)
                {
                  flag = true;
                  //do the reinsertion on the real mesh
                  PointInsertion pi(m_simplexMesh, SimplicesNode(m_simplexMesh, node), criterionRAIS, status, cavity, markedNodes, deletedNodes, facesAlreadyBuilt, cellsCreated);
                  gmds::ISimplexMeshIOService ioServiceMesh(m_simplexMesh);
                  gmds::VTKWriter vtkWriterMesh(&ioServiceMesh);
                  vtkWriterMesh.setCellOptions(gmds::N|gmds::R);
                  vtkWriterMesh.setDataOptions(gmds::N|gmds::R);
                  vtkWriterMesh.write("FACE_ADAPTAION_TEST_" + std::to_string(cpt) + ".vtk");
                  cpt++;
                }
                break;
                /*
                std::cout << "worstQuality -> " << worstQuality << std::endl;
                std::cout << "qual -> " << qual << std::endl;
                m_simplexMesh->deleteAllSimplicesBut(cavity);
                simplexMeshCopy.deleteAllSimplicesBut(cellsCreated);
                gmds::ISimplexMeshIOService ioServiceMesh(m_simplexMesh);
                gmds::VTKWriter vtkWriterMesh(&ioServiceMesh);
                vtkWriterMesh.setCellOptions(gmds::N|gmds::R);
                vtkWriterMesh.setDataOptions(gmds::N|gmds::R);
                vtkWriterMesh.write("TEST_tet_Before.vtk");

                gmds::ISimplexMeshIOService ioServiceMeshCopy(&simplexMeshCopy);
                gmds::VTKWriter vtkWriterMeshCopy(&ioServiceMeshCopy);
                vtkWriterMeshCopy.setCellOptions(gmds::N|gmds::R);
                vtkWriterMeshCopy.setDataOptions(gmds::N|gmds::R);
                vtkWriterMeshCopy.write("TEST_tet_After.vtk");
                return;*/
              }
            }
          }
        }
      }
    }
  }
}
