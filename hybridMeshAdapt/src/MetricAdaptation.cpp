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
unsigned int MetricAdaptation::computeSlicing()
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
  unsigned int cptSlice = 0;

  //Metric adaptation SLICING
  buildEdgesMap();
  for(auto const edge : m_edgesMap)
  {
    TInt nodeId = edge.first;
    TInt dNode  = edge.second;
    if(meshNode[nodeId] != 0 && meshNode[dNode] != 0)
    {
      SimplicesNode node0(m_simplexMesh, nodeId);
      SimplicesNode node1(m_simplexMesh, dNode);
      std::vector<TSimplexID> shell = node0.shell(node1);
      if(shell.size() == 0){continue;}

      math::Point nodeCoord0 = node0.getCoords();
      math::Point nodeCoord1 = node1.getCoords();
      Metric<Eigen::Matrix3d> M0 = Metric<Eigen::Matrix3d>((*metric)[nodeId]);
      Metric<Eigen::Matrix3d> M1 = Metric<Eigen::Matrix3d>((*metric)[dNode]);
      double metricLenght        = M0.metricDist(nodeCoord0, nodeCoord1, M1);


      unsigned int nodeDim = ((*BND_VERTEX_COLOR)[nodeId] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[nodeId] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[nodeId] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
      unsigned int nodeLabel = ((*BND_VERTEX_COLOR)[nodeId] != 0)?(*BND_VERTEX_COLOR)[nodeId]:((*BND_CURVE_COLOR)[nodeId] != 0)?(*BND_CURVE_COLOR)[nodeId]:((*BND_SURFACE_COLOR)[nodeId] != 0)?(*BND_SURFACE_COLOR)[nodeId]:0;
      unsigned int directNodeDim = ((*BND_VERTEX_COLOR)[dNode] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[dNode] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[dNode] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
      unsigned int directNodeLabel = ((*BND_VERTEX_COLOR)[dNode] != 0)?(*BND_VERTEX_COLOR)[dNode]:((*BND_CURVE_COLOR)[dNode] != 0)?(*BND_CURVE_COLOR)[dNode]:((*BND_SURFACE_COLOR)[dNode] != 0)?(*BND_SURFACE_COLOR)[dNode]:0;
      unsigned int newNodeDim = std::max(nodeDim, directNodeDim);

      if(metricLenght > sqrt(2.0))
      {
        //std::cout << " node0 -> " << nodeId << " | dimNode -> " << nodeDim << std::endl;
        //std::cout << " node1 -> " << dNode << " | dimNode -> " << directNodeDim << std::endl;
        //std::cout << " metricLenght -> " << metricLenght << std::endl;
        bool alreadyAdd = false;
        std::vector<TSimplexID> tetraContenaingPt{};
        const Point pt = 0.5 * (nodeCoord0 + nodeCoord1);
        TInt newNodeId = m_simplexMesh->addNodeAndcheck(pt, tetraContenaingPt, alreadyAdd);

        if(alreadyAdd){continue;}

        metric->set(newNodeId, (*metric)[nodeId]);
        //labelization of the node being inserted
        if(newNodeDim == SimplexMesh::topo::CORNER)
        {
          //TODO
          //define a new curve color...
          //BND_CURVE_COLOR->set(newNodeId, newColorToDefine)
          continue;
        }
        else if(newNodeDim == SimplexMesh::topo::RIDGE)
        {
          if(nodeDim == directNodeDim)
          {
            if(nodeLabel != directNodeLabel)
            {
              //node sur la surface...
              const std::map<unsigned int, std::pair<unsigned int, unsigned int>>& ridgeEdge = m_simplexMesh->getEdgeTianglesIndices();
              std::pair<unsigned int, unsigned int> p0 = ridgeEdge.at(nodeLabel);
              std::pair<unsigned int, unsigned int> p1 = ridgeEdge.at(directNodeLabel);
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
                }
              }

              if(!flag)
              {
                std::cout << "edge To Insert -> " << nodeId << " | " << dNode << std::endl;
                std::cout << "p0.first | p0.second -> " << p0.first << " | " << p0.second << std::endl;
                std::cout << "p1.first | p1.second -> " << p1.first << " | " << p1.second << std::endl;

                SimplexMesh nodeMesh = SimplexMesh();
                nodeMesh.addNode(SimplicesNode(m_simplexMesh, nodeId).getCoords());
                nodeMesh.addNode(SimplicesNode(m_simplexMesh, dNode).getCoords());
                nodeMesh.addTetraedre(0, 0, 0, 0);
                gmds::ISimplexMeshIOService ioServiceNode(&nodeMesh);
                gmds::VTKWriter vtkWriterNode(&ioServiceNode);
                vtkWriterNode.setCellOptions(gmds::N|gmds::R);
                vtkWriterNode.setDataOptions(gmds::N|gmds::R);
                vtkWriterNode.write("NODE_" + std::to_string(nodeId) + "_" + std::to_string(dNode) +  ".vtk");

                gmds::VTKWriter vtkWriterDEBUG(&ioServiceMesh);
                vtkWriterDEBUG.setCellOptions(gmds::N|gmds::R|gmds::F);
                vtkWriterDEBUG.setDataOptions(gmds::N|gmds::R|gmds::F);
                vtkWriterDEBUG.write("No_COMMON_SURFACE.vtk");
                throw gmds::GMDSException("no common surface");
              }
            }
            else
            {
              BND_CURVE_COLOR->set(newNodeId, nodeLabel);
            }
          }
          else
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
          else
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
                else
                {
                  std::cout << "nodeLabel -> " << nodeLabel << std::endl;
                  std::cout << "directNodeLabel -> " << directNodeLabel << std::endl;
                  std::cout << "continue -> " << p1.first << " | " << p1.second << std::endl;
                  for(auto const connexion : ridgeEdge)
                  {
                    std::cout << "CONNEXION -> " << connexion.first << " | " << connexion.second.first << " & " << connexion.second.second << std::endl;
                  }
                  throw gmds::GMDSException("surface no linked to the ridge");
                  //continue;
                }
              }
              else if(directNodeDim == SimplexMesh::topo::CORNER)
              {
                BND_SURFACE_COLOR->set(newNodeId, nodeLabel);
              }
            }
            else // nodeDim < directNodeDim
            {
              if(nodeDim == SimplexMesh::topo::RIDGE)
              {
                const std::map<unsigned int, std::pair<unsigned int, unsigned int>>& ridgeEdge = m_simplexMesh->getEdgeTianglesIndices();
                std::pair<unsigned int, unsigned int> p1 = ridgeEdge.at(nodeLabel);
                std::vector<unsigned int> p1Vec{p1.first, p1.second};
                if(p1Vec.front() == directNodeLabel || p1Vec.back() == directNodeLabel)
                {
                  BND_SURFACE_COLOR->set(newNodeId, directNodeLabel);
                }
                else
                {
                  continue;
                }
              }
              else if(nodeDim == SimplexMesh::topo::CORNER)
              {
                BND_SURFACE_COLOR->set(newNodeId, directNodeLabel);
              }
              //BND_SURFACE_COLOR->set(newNodeId, directNodeLabel);
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

        if(status)
        {
          cptSlice++;
          //if(cptSlice == 1)
          {
            //std::cout << "    INSERTION " << std::endl;
            //gmds::VTKWriter vtkWriterMesh(&ioServiceMesh);
            //vtkWriterMesh.setCellOptions(gmds::N|gmds::R|gmds::F);
            //vtkWriterMesh.setDataOptions(gmds::N|gmds::R|gmds::F);
            //vtkWriterMesh.write("SLICING_TEST_" + std::to_string(cptSlice) + ".vtk");
          }
        }
      }
    }
  }
  cptTest++;
  return cptSlice;
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

    //std::cout << " node0 -> " << nodeId << " | dimNode -> " << nodeDim << std::endl;
    //std::cout << " node1 -> " << dNode << " | dimNode -> " << directNodeDim << std::endl;
    //std::cout << " metricLenght -> " << metricLenght << std::endl;
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
    }
    else if(newNodeDim == SimplexMesh::topo::RIDGE)
    {
      if(nodeDim == directNodeDim)
      {
        if(nodeLabel != directNodeLabel)
        {
          //node sur la surface or the volume...
          const std::map<unsigned int, std::pair<unsigned int, unsigned int>>& ridgeEdge = m_simplexMesh->getEdgeTianglesIndices();
          std::pair<unsigned int, unsigned int> p0 = ridgeEdge.at(nodeLabel);
          std::pair<unsigned int, unsigned int> p1 = ridgeEdge.at(directNodeLabel);
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
            }
          }

          //node à inserer dans le volume donc rien a faire
          if(!flag)
          {
            /*std::cout << "edge To Insert -> " << nodeA << " | " << nodeB << std::endl;
            std::cout << "p0.first | p0.second -> " << p0.first << " | " << p0.second << std::endl;
            std::cout << "p1.first | p1.second -> " << p1.first << " | " << p1.second << std::endl;

            SimplexMesh nodeMesh = SimplexMesh();
            nodeMesh.addNode(SimplicesNode(m_simplexMesh, nodeA).getCoords());
            nodeMesh.addNode(SimplicesNode(m_simplexMesh, nodeB).getCoords());
            nodeMesh.addTetraedre(0, 0, 0, 0);
            gmds::ISimplexMeshIOService ioServiceNode(&nodeMesh);
            gmds::VTKWriter vtkWriterNode(&ioServiceNode);
            vtkWriterNode.setCellOptions(gmds::N|gmds::R);
            vtkWriterNode.setDataOptions(gmds::N|gmds::R);
            vtkWriterNode.write("NODE_" + std::to_string(nodeA) + "_" + std::to_string(nodeB) +  ".vtk");

            gmds::VTKWriter vtkWriterDEBUG(&ioServiceMesh);
            vtkWriterDEBUG.setCellOptions(gmds::N|gmds::R|gmds::F);
            vtkWriterDEBUG.setDataOptions(gmds::N|gmds::R|gmds::F);
            vtkWriterDEBUG.write("No_COMMON_SURFACE.vtk");
            throw gmds::GMDSException("no common surface");*/
          }
        }
        else
        {
          BND_CURVE_COLOR->set(newNodeId, nodeLabel);
        }
      }
      else
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
      else
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
            else
            {
              //The node is on the volume

              /*std::cout << "nodeA -> " << nodeA << " | dimNode -> " << nodeDim << " | label -> " << nodeLabel << std::endl;
              std::cout << "nodeB -> " << nodeB << " | dimNode -> " << directNodeDim << " | label -> " << directNodeLabel << std::endl;
              std::cout << "nodeLabel -> " << nodeLabel << std::endl;
              std::cout << "directNodeLabel -> " << directNodeLabel << std::endl;
              std::cout << "continue -> " << p1.first << " | " << p1.second << std::endl;
              for(auto const connexion : ridgeEdge)
              {
                std::cout << "CONNEXION -> " << connexion.first << " | " << connexion.second.first << " & " << connexion.second.second << std::endl;
              }
              throw gmds::GMDSException("surface no linked to the ridge");*/
              //continue;
            }
          }
          else if(directNodeDim == SimplexMesh::topo::CORNER)
          {
            BND_SURFACE_COLOR->set(newNodeId, nodeLabel);
          }
        }
        else // nodeDim < directNodeDim
        {
          if(nodeDim == SimplexMesh::topo::RIDGE)
          {
            const std::map<unsigned int, std::pair<unsigned int, unsigned int>>& ridgeEdge = m_simplexMesh->getEdgeTianglesIndices();
            std::pair<unsigned int, unsigned int> p1 = ridgeEdge.at(nodeLabel);
            std::vector<unsigned int> p1Vec{p1.first, p1.second};
            if(p1Vec.front() == directNodeLabel || p1Vec.back() == directNodeLabel)
            {
              BND_SURFACE_COLOR->set(newNodeId, directNodeLabel);
            }
            else
            {
              return false;
            }
          }
          else if(nodeDim == SimplexMesh::topo::CORNER)
          {
            BND_SURFACE_COLOR->set(newNodeId, directNodeLabel);
          }
          //BND_SURFACE_COLOR->set(newNodeId, directNodeLabel);
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
unsigned int MetricAdaptation::computeEdgeRemove()
{
  std::vector<TSimplexID> v{};
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
  //gmds::ISimplexMeshIOService ioServiceMesh(m_simplexMesh);

  CriterionRAIS criterionRAIS(new VolumeCriterion());
  const gmds::BitVector& meshNode = m_simplexMesh->getBitVectorNodes();
  unsigned int cptEdgeRemove = 0;

  buildEdgesMap();
  for(auto const edge : m_edgesMap)
  {
    TInt node0 = edge.first;
    TInt node1 = edge.second;
    //std::cout << "nodes -> " << node0 << " | " << node1 << " --> " << meshNode[node0] << " | " << meshNode[node1] <<  std::endl;
    if(meshNode[node0] != 0 && meshNode[node1] != 0)
    {
      unsigned int nodeDim = ((*BND_VERTEX_COLOR)[node0] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[node0] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[node0] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
      unsigned int directNodeDim = ((*BND_VERTEX_COLOR)[node1] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[node1] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[node1] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
      //we do not remove and edge with a node on the surface and node on a curve
      //std::cout << "node0 label -> " << nodeDim << " | node1 label ->" << directNodeDim <<   std::endl;

      if((nodeDim == SimplexMesh::topo::RIDGE))
      {
        if((directNodeDim == SimplexMesh::topo::SURFACE))
        {
          continue;
        }
      }
      if((nodeDim == SimplexMesh::topo::SURFACE))
      {
        if((directNodeDim == SimplexMesh::topo::RIDGE))
        {
          continue;
        }
      }
      //
      if((nodeDim == SimplexMesh::topo::SURFACE))
      {
        if((directNodeDim == SimplexMesh::topo::VOLUME))
        {
          continue;
        }
      }
      if((nodeDim == SimplexMesh::topo::VOLUME))
      {
        if((directNodeDim == SimplexMesh::topo::SURFACE))
        {
          continue;
        }
      }
      //
      if((nodeDim == SimplexMesh::topo::RIDGE))
      {
        if((directNodeDim == SimplexMesh::topo::VOLUME))
        {
          continue;
        }
      }
      if((nodeDim == SimplexMesh::topo::VOLUME))
      {
        if((directNodeDim == SimplexMesh::topo::RIDGE))
        {
          continue;
        }
      }
      //
      unsigned int newNodeDim = std::max(nodeDim, directNodeDim);
      //if(newNodeDim > 2){continue;}
      //
      SimplicesNode sNode0(m_simplexMesh, node0);
      SimplicesNode sNode1(m_simplexMesh, node1);
      std::vector<TSimplexID> shell = sNode0.shell(sNode1);

      //std::cout << "shell.size() -> " << shell.size() << std::endl;
      if(shell.size() == 0){continue;}

      math::Point nodeCoord0 = sNode0.getCoords();
      math::Point nodeCoord1 = sNode1.getCoords();
      Metric<Eigen::Matrix3d> M0 = Metric<Eigen::Matrix3d>((*metric)[node0]);
      Metric<Eigen::Matrix3d> M1 = Metric<Eigen::Matrix3d>((*metric)[node1]);
      double metricLenght        = M0.metricDist(nodeCoord0, nodeCoord1, M1);
      //std::cout << "metricLenght -> " << metricLenght << std::endl;
      if(metricLenght < (sqrt(2.0) * 0.5))
      {
        //v.push_back(1);
        const std::vector<TSimplexID>& base = m_simplexMesh->getBase();
        const gmds::BitVector& bitVectorTet = m_simplexMesh->getBitVectorTet();
        const gmds::BitVector& bitVectorNode = m_simplexMesh->getBitVectorNodes();
        //std::cout << "  edgeRemove( " <<  node0 << " | " << node1  << " ) dimension -> " << nodeDim << " | " << directNodeDim << "  | label -> " << << std::endl;
        if(m_simplexMesh->edgeRemove(node0, node1))
        {
          //std::cout << "edge remove " << std::endl;
          /*gmds::VTKWriter vtkWriterMeshER(&ioServiceMesh);
          vtkWriterMeshER.setCellOptions(gmds::N|gmds::R|gmds::F);
          vtkWriterMeshER.setDataOptions(gmds::N|gmds::R|gmds::F);
          vtkWriterMeshER.write("TEST_BETA_" + std::to_string(cptt) + ".vtk");
          std::cout << "cpt -> " << cptt << std::endl;
          cptt++;*/
          cptEdgeRemove++;
          //std::cout << std::endl;
          //std::cout << std::endl;
        }
      }
    }
  }

  return cptEdgeRemove;
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
    return m_simplexMesh->edgeRemove(nodeA, nodeB);
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

  for(TInt nodeId = 0 ; nodeId < meshNode.capacity() ; nodeId++)
  {
    if(meshNode[nodeId] != 0)
    {
      std::cout << std::endl;
      unsigned int nodeIdDim   = ((*BND_VERTEX_COLOR)[nodeId] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[nodeId] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[nodeId] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
      unsigned int nodeIdLabel = (nodeIdDim == SimplexMesh::topo::CORNER)?(*BND_VERTEX_COLOR)[nodeId]:(nodeIdDim == SimplexMesh::topo::RIDGE)?(*BND_CURVE_COLOR)[nodeId]:(nodeIdDim == SimplexMesh::topo::SURFACE)?(*BND_SURFACE_COLOR)[nodeId]:0;
      //std::cout << std::endl;
      //std::cout << "nodeId -> " << nodeId << std::endl;
      //std::cout << "nodeIdDim -> " << nodeIdDim << std::endl;
      //std::cout << "nodeIdLabel -> " << nodeIdLabel << std::endl;

      if(nodeIdDim == SimplexMesh::topo::SURFACE || nodeIdDim == SimplexMesh::topo::RIDGE)
      {
        continue;
      }

      //computing the optimal position of P in order to have Lm(PPj) ~ 1
      //https://tel.archives-ouvertes.fr/tel-00120327/document
      Metric<Eigen::Matrix3d> M0 = Metric<Eigen::Matrix3d>((*metric)[nodeId]);
      const SimplicesNode   node = SimplicesNode(m_simplexMesh, nodeId);
      const math::Point        p = node.getCoords();
      const std::vector<TInt> nodes = node.neighborSubSurfaceNodes();
      std::vector<math::Point> P_opti{};
      std::vector<TSimplexID> ball = node.ballOf();
      if(ball.size() == 0){continue;}

      for(auto const simplex : ball)
      {
        if(simplex >= 0)
        {
          unsigned int sizeTet = 4;
          unsigned int localNode = SimplicesCell(m_simplexMesh, simplex).getLocalNode(nodeId);
          math::Point p_opti(0.0, 0.0, 0.0);

          for(unsigned int n = 1 ; n < sizeTet ; n++)
          {
            unsigned int nodeFaceLocal = (localNode + n) % sizeTet;
            const SimplicesNode sNode = SimplicesCell(m_simplexMesh, simplex).getNode(nodeFaceLocal);
            const math::Point pj        = sNode.getCoords();
            Metric<Eigen::Matrix3d> M1  = Metric<Eigen::Matrix3d>((*metric)[sNode.getGlobalNode()]);
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
          P_opti.push_back(p_opti * (1.0 / 3.0));
        }
      }

      math::Point new_P;
      if(P_opti.size() != 0)
      {
        for(auto const p : P_opti)
        {
          new_P = new_P + p;
        }
        new_P = math::Point(new_P.X() / (double)P_opti.size(), new_P.Y() / (double)P_opti.size(), new_P.Z() / (double)P_opti.size());

        //test pour calculer si new_P est étoilable par rapport a la boule de nodeId
        std::set<double> tetQuality{};
        for(auto const simplex : ball)
        {
          if(simplex >= 0)
          {
            //compute quality of the tetrahedron for the next step
            tetQuality.insert(m_simplexMesh->computeQualityElement(simplex));
            unsigned int localNode = SimplicesCell(m_simplexMesh, simplex).getLocalNode(nodeId);
            if(criterionRAIS.execute(m_simplexMesh, simplex, localNode, new_P))
            {
              return computePointSmoothing;
            }
          }
        }

        //test pour calculer le plus pire tetra créer est meilleur que le plus pire de l'ancienne boule liée a nodeId
        std::set<double> newQuality{};
        for(auto const simplex : ball)
        {
          if(simplex >= 0)
          {
            unsigned int sizeTet = 4;
            unsigned int localNode = SimplicesCell(m_simplexMesh, simplex).getLocalNode(nodeId);
            std::vector<TInt> faceNodes = SimplicesCell(m_simplexMesh, simplex).getOrderedFace(localNode);

            const math::Point pA = SimplicesNode(m_simplexMesh, faceNodes.front()).getCoords();
            const math::Point pB = SimplicesNode(m_simplexMesh, faceNodes[1]).getCoords();
            const math::Point pC = SimplicesNode(m_simplexMesh, faceNodes.back()).getCoords();

            Eigen::Matrix3d MA = (*metric)[faceNodes.front()];
            Eigen::Matrix3d MB = (*metric)[faceNodes[1]];
            Eigen::Matrix3d MC = (*metric)[faceNodes.back()];

            Eigen::Matrix3d M_newpt =  Eigen::MatrixXd::Identity(3, 3);
            double metricX = 0.1;
            double metricY = 0.1;
            double metricZ = 0.1;

            M_newpt(0,0) = 1.0 / (metricX*metricX);
            M_newpt(1,1) = 1.0 / (metricY*metricY);
            M_newpt(2,2) = 1.0 / (metricZ*metricZ);

            double qualityNewTet = m_simplexMesh->computeQualityElement(pA, pB, pC, new_P, MA, MB, MC, M_newpt);
            //if a new tet vuild is worst than the worst previous tet we do not rebuild the cav
            if(qualityNewTet < *(tetQuality.begin()))
            {
              return computePointSmoothing;
            }
          }
        }
      }
      m_simplexMesh->moveNodeCoord(nodeId, new_P);
    }
  }

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
  unsigned int maxIterationAdaptation = 100;
  CriterionRAIS criterionRAIS(new VolumeCriterion());
  gmds::ISimplexMeshIOService ioServiceMesh(m_simplexMesh);
  const gmds::BitVector& meshNode = m_simplexMesh->getBitVectorNodes();

  //ADAPTATION
  std::cout << "ADAPTATION  START" << std::endl;
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
      //std::cout << "nodes -> " << nodeA << " | " << nodeB << std::endl;
      //if((*BND_VERTEX_COLOR)[nodeA] == 0 && (*BND_SURFACE_COLOR)[nodeA] == 0 && (*BND_CURVE_COLOR)[nodeA] == 0){continue;}
      //if((*BND_VERTEX_COLOR)[nodeB] == 0 && (*BND_SURFACE_COLOR)[nodeB] == 0 && (*BND_CURVE_COLOR)[nodeB] == 0){continue;}
      if((*BND_SURFACE_COLOR)[nodeA] == (*BND_SURFACE_COLOR)[nodeB]  && (*BND_SURFACE_COLOR)[nodeB] == 1){
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
        }
      }
    }

    gmds::VTKWriter vtkWriterCS(&ioServiceMesh);
    vtkWriterCS.setCellOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterCS.setDataOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterCS.write("ComputeSlicing_EdgeRemove_" + std::to_string(cpt) + ".vtk");
    /*std::multimap<TInt, std::pair<TInt,TInt>>& edgeStructure = m_simplexMesh->getEdgeStructure();
    for(auto const p : edgeStructure)
    {
      std::cout << "edge idx -> " << p.first << " | nodes -> " << p.second.first << " : " << p.second.second << std::endl;
    }*/
    //throw gmds::GMDSException("END");

    //BOUGE DE POINTS
    /*std::cout << "POINT SMOOTHING START" << std::endl;
    cptSmoothing = computePointSmoothing();
    std::cout << "    cptSmoothing -> " << cptSmoothing << std::endl;
    gmds::VTKWriter vtkWriterPS(&ioServiceMesh);
    vtkWriterPS.setCellOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterPS.setDataOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterPS.write("ComputePoint_Smoothing_" + std::to_string(cpt) + ".vtk");*/

    /*const gmds::BitVector& tetQuality = m_simplexMesh->getBitVectorTet();
    std::set<double> s{};
    for(unsigned int tet = 0 ; tet < tetQuality.capacity() ; tet++)
    {
      if(tetQuality[tet] != 0)
      {
        double qualityelmt = m_simplexMesh->computeQualityElement(tet);
        s.insert(qualityelmt);
        std::cout << "tet -> " << tet << std::endl;
        std::cout << "qualityelmt -> " << qualityelmt << std::endl;
      }
    }
    std::cout << "BEST quality Element -> " << *(--s.end()) << std::endl;
    std::vector<TSimplexID> tet{6178};
    m_simplexMesh->deleteAllSimplicesBut(tet);
    gmds::VTKWriter vtkWriterPS(&ioServiceMesh);
    vtkWriterPS.setCellOptions(gmds::N|gmds::R);
    vtkWriterPS.setDataOptions(gmds::N|gmds::R);
    vtkWriterPS.write("TEST.vtk");
    throw gmds::GMDSException("END");*/
    /*for(auto const qual : s)
    {
      std::cout << "qualityelmt -> " << qual << std::endl;
    }*/

    /*std::cout << "  SWAPPING EDGE START" << std::endl;
    cptEdgeSwap = computeSurfaceEdgeSwap();
    std::cout << "    cptEdgeSwap -> " << cptEdgeSwap << std::endl;
    gmds::VTKWriter vtkWriterMeshSwapping(&ioServiceMesh);
    vtkWriterMeshSwapping.setCellOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterMeshSwapping.setDataOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterMeshSwapping.write("ComputeSurfaceSwap_" + std::to_string(cpt) + ".vtk");

    std::multimap<TInt, std::pair<TInt,TInt>>& edgeStructure1 = m_simplexMesh->getEdgeStructure();
    for(auto const p : edgeStructure1)
    {
      std::cout << "edge idx -> " << p.first << " | nodes -> " << p.second.first << " : " << p.second.second << std::endl;
    }
    throw gmds::GMDSException("END");*/
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
    if(cpt == 10)
    {
      //throw gmds::GMDSException("END");
    }
    /*m_simplexMesh->getEdgeSizeInfowithMetric(meanEdge, minEdge, maxEdge);
    std::cout << "      meanEdge END -> " << meanEdge << std::endl;
    std::cout << "      minEdge  END -> " << minEdge << std::endl;
    std::cout << "      maxEdge  END -> " << maxEdge << std::endl;*/
  }
/*----------------------------------------------------------------------------*/
/*void MetricAdaptation::execute()
{
  unsigned int cpt = 0;
  unsigned int maxIterationAdaptation = 100;
  CriterionRAIS criterionRAIS(new VolumeCriterion());
  gmds::ISimplexMeshIOService ioServiceMesh(m_simplexMesh);

  //ADAPTATION
  std::cout << "ADAPTATION  START" << std::endl;
  for(unsigned int iter = 0 ; iter < maxIterationAdaptation ; iter++)
  {
    std::cout << std::endl;
    std::cout << " cpt -> " << cpt << std::endl;

    unsigned int cptSlice = 0;
    unsigned int cptEdgeRemove = 0;
    unsigned int cptSwap = 0;
    double meanEdge = 0.0; double minEdge = 0.0; double maxEdge = 0.0;

    //Metric adaptation SLICING
    std::cout << "  SLICING START" << std::endl;
    cptEdgeRemove = computeEdgeRemove();

    unsigned int previousCptSlice = cptSlice;
    do{
      previousCptSlice = cptSlice;
      cptSlice = computeSlicing();
      std::cout << "    cptSlice -> " << cptSlice << std::endl;
      //cptEdgeRemove = computeEdgeRemove();
      if(cptSlice == 1 & cptSlice == previousCptSlice)
      {
        //break;
      }
    }while(cptSlice != 0);
    gmds::VTKWriter vtkWriterMesh(&ioServiceMesh);
    vtkWriterMesh.setCellOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterMesh.setDataOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterMesh.write("ADAPTAION_LOOP_SLICING_" + std::to_string(cpt) + ".vtk");
    m_simplexMesh->getEdgeSizeInfo(meanEdge, minEdge, maxEdge);
    std::cout << "      meanEdge END -> " << meanEdge << std::endl;
    std::cout << "      minEdge  END -> " << minEdge << std::endl;
    std::cout << "      maxEdge  END -> " << maxEdge << std::endl;

    //Metric adaptation EDGE REMOVE
    std::cout << "  EDGE REMOVE START" << std::endl;
    do{
      cptEdgeRemove = computeEdgeRemove();
      std::cout << "    cptEdgeRemove -> " << cptEdgeRemove << std::endl;
    }while(cptEdgeRemove != 0);
    gmds::VTKWriter vtkWriterMeshER(&ioServiceMesh);
    vtkWriterMeshER.setCellOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterMeshER.setDataOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterMeshER.write("ADAPTAION_LOOP_EDGEREMOVE_" + std::to_string(cpt) + ".vtk");
    m_simplexMesh->getEdgeSizeInfo(meanEdge, minEdge, maxEdge);
    std::cout << "      meanEdge END -> " << meanEdge << std::endl;
    std::cout << "      minEdge  END -> " << minEdge << std::endl;
    std::cout << "      maxEdge  END -> " << maxEdge << std::endl;


    //Swapping face
    std::cout << "  SWAPPING FACE START" << std::endl;
    do{
      cptSwap = computeFaceSwap();
      std::cout << "    cptSwap -> " << cptSwap << std::endl;
    }while(cptSwap != 0);
    gmds::VTKWriter vtkWriterMeshSwapping(&ioServiceMesh);
    vtkWriterMeshSwapping.setCellOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterMeshSwapping.setDataOptions(gmds::N|gmds::R|gmds::F);
    vtkWriterMeshSwapping.write("ADAPTAION_LOOP_SWAPPING_" + std::to_string(cpt) + ".vtk");
    m_simplexMesh->getEdgeSizeInfo(meanEdge, minEdge, maxEdge);
    std::cout << "      meanEdge END -> " << meanEdge << std::endl;
    std::cout << "      minEdge  END -> " << minEdge << std::endl;
    std::cout << "      maxEdge  END -> " << maxEdge << std::endl;

    cpt++;

    //BOUGE DE POINTS
    //cpt++;
  }*/
}
