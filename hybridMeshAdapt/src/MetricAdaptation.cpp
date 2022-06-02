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
  //ADAPTATION face
  std::cout << "ADAPTATION FACE START" << std::endl;
  gmds::ISimplexMeshIOService ioServiceMesh(m_simplexMesh);

  for(unsigned int iter = 0 ; iter < maxIterationAdaptation ; iter++)
  {
    std::cout << "MESHNODE.CAPACITY() -> " << meshNode.capacity() << std::endl;

    unsigned int cptEdgeRemove = 0;
    unsigned int cptSlice = 0;
    unsigned int cptSwap = 0;

    //Metric adaptation SLICING
    std::cout << "  SLICING START" << std::endl;
    /*buildEdgesMap();
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
                    continue;
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

          //gmds::VTKWriter vtkWriterMesh(&ioServiceMesh);
          //vtkWriterMesh.setCellOptions(gmds::N|gmds::R|gmds::F);
          //vtkWriterMesh.setDataOptions(gmds::N|gmds::R|gmds::F);
          //vtkWriterMesh.write("LOOP_" + std::to_string(cpt) + ".vtk");
          //cpt++;
          if(status)
          {
            cptSlice++;
          }
        }
      }
    }*/
    /*gmds::VTKWriter vtkWriterMesh(&ioServiceMesh);
    vtkWriterMesh.setCellOptions(gmds::N|gmds::R);
    vtkWriterMesh.setDataOptions(gmds::N|gmds::R);
    vtkWriterMesh.write("ADAPTAION_LOOP_SLICING_" + std::to_string(cpt) + ".vtk");
    cpt++;*/

    std::cout << "  SLINCING END" << std::endl;
    //continue;

    //Metric adaptation EDGE REMOVE
    std::cout << "  EDGE REMOVE START" << std::endl;
    buildEdgesMap();
    const TInt border = std::numeric_limits<TInt>::min();
    for(auto const edge : m_edgesMap)
    {
      TInt node0 = edge.first;
      TInt node1 = edge.second;
      if(meshNode[node0] != 0 && meshNode[node1] != 0)
      {
        unsigned int nodeDim = ((*BND_VERTEX_COLOR)[node0] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[node0] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[node0] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
        unsigned int directNodeDim = ((*BND_VERTEX_COLOR)[node1] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[node1] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[node1] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
        //
        unsigned int newNodeDim = std::max(nodeDim, directNodeDim);
        if(newNodeDim > 2){continue;}
        //
        SimplicesNode sNode0(m_simplexMesh, node0);
        SimplicesNode sNode1(m_simplexMesh, node1);
        std::vector<TSimplexID> shell = sNode0.shell(sNode1);

        if(shell.size() == 0){continue;}

        math::Point nodeCoord0 = sNode0.getCoords();
        math::Point nodeCoord1 = sNode1.getCoords();
        Metric<Eigen::Matrix3d> M0 = Metric<Eigen::Matrix3d>((*metric)[node0]);
        Metric<Eigen::Matrix3d> M1 = Metric<Eigen::Matrix3d>((*metric)[node1]);
        double metricLenght        = M0.metricDist(nodeCoord0, nodeCoord1, M1);
        if(metricLenght < (1.0 / sqrt(2.0)))
        {
          const std::vector<TSimplexID>& base = m_simplexMesh->getBase();
          const gmds::BitVector& bitVectorTet = m_simplexMesh->getBitVectorTet();
          const gmds::BitVector& bitVectorNode = m_simplexMesh->getBitVectorNodes();
          std::cout << "base[439] -> " << base[439] << std::endl;
          std::cout << "bitVectorTet[base[439]] -> " << bitVectorTet[base[439]] << std::endl;
          if(bitVectorNode[439] == 1 && bitVectorTet[base[439]] == 0)
          {
            return;
          }
          std::cout << "node0 node1 -> " << node0 << " | " << node1 << std::endl;

          if(m_simplexMesh->edgeRemove(node0, node1))
          {
            cptEdgeRemove++;
            gmds::VTKWriter vtkWriterMeshER(&ioServiceMesh);
            vtkWriterMeshER.setCellOptions(gmds::N|gmds::R|gmds::F);
            vtkWriterMeshER.setDataOptions(gmds::N|gmds::R|gmds::F);
            vtkWriterMeshER.write("ADAPTAION_LOOP_EDGEREMOVE_" + std::to_string(cptEdgeRemove) + ".vtk");
            std::cout << "cptEdgeRemove -> " << cptEdgeRemove << std::endl;
          }
        }
      }
    }
    std::cout << "  EDGE REMOVE END" << std::endl;
    continue;

    /*gmds::VTKWriter vtkWriterMeshER(&ioServiceMesh);
    vtkWriterMeshER.setCellOptions(gmds::N|gmds::R);
    vtkWriterMeshER.setDataOptions(gmds::N|gmds::R);
    vtkWriterMeshER.write("ADAPTAION_LOOP_EDGEREMOVE_" + std::to_string(cpt) + ".vtk");*/


    //Swapping face
    /*std::cout << "  SWAPPING FACE START" << std::endl;
    const gmds::BitVector& tetBitVec = m_simplexMesh->getBitVectorTet();
    unsigned int sizeFace = 4;
    for(TSimplexID tetId = 0 ; tetId < tetBitVec.capacity() ; tetId++)
    {
      if(tetBitVec[tetId] != 0)
      {
        bool status = false;
        for(unsigned int face = 0 ; face < sizeFace ; face++)
        {
          if(status)
          {
            break;
          }
          SimplicesCell cell0(m_simplexMesh, tetId);
          const TInt node_0 = cell0.getNodes()[face];
          TSimplexID oppositeTet = cell0.oppositeTetraIdx(face);
          if(oppositeTet < 0)//it is a triangle
          {
            continue;
          }
          else
          {
            SimplicesCell cell1(m_simplexMesh, oppositeTet);
            //check if the swapping is possible
            std::vector<TInt> visibleFaces = cell1.visibleFaces(SimplicesNode(m_simplexMesh, node_0).getCoords());
            if(visibleFaces.size() == 3)
            {
              std::vector<TSimplexID> cavity{tetId, oppositeTet};
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

                  //std::cout << "element quality factor -> " << worstQuality << std::endl;
                  /std::cout << "new element quality factor -> " << worstNewElement << std::endl;
                  //std::cout << std::endl;

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
                      cptSwap++;
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
    }*/
    std::cout << "  SWAPPING FACE END" << std::endl;
    std::cout << "cptEdgeRemove -> " << cptEdgeRemove << " | " << "cptSlice -> " << cptSlice << " | " << "cptSwap -> " << cptSwap << std::endl;

    std::cout << std::endl;
    gmds::VTKWriter vtkWriterMeshSwapping(&ioServiceMesh);
    vtkWriterMeshSwapping.setCellOptions(gmds::N|gmds::R);
    vtkWriterMeshSwapping.setDataOptions(gmds::N|gmds::R);
    vtkWriterMeshSwapping.write("ADAPTAION_LOOP_SWAPPING_" + std::to_string(cpt) + ".vtk");
    cpt++;
    //BOUGE DE POINTS
    /*std::cout << "BOUGE DE POINTS" << std::endl;
    const gmds::BitVector& meshNode = m_simplexMesh->getBitVectorNodes();
    const gmds::BitVector& meshTet  = m_simplexMesh->getBitVectorTet();
    gmds::BitVector nodesAdded(meshNode.capacity());
    unsigned int cptBouge = 0;
    for(TInt nodeId = 0 ; nodeId < meshNode.capacity() ; nodeId++)
    {
      if(meshNode[nodeId] != 0 && nodesAdded[nodeId] == 0)
      {
        //std::cout << "nodeId -> " << nodeId << std::endl;
        unsigned int nodeIdDim   = ((*BND_VERTEX_COLOR)[nodeId] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[nodeId] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[nodeId] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
        unsigned int nodeIdLabel = (nodeIdDim == SimplexMesh::topo::CORNER)?(*BND_VERTEX_COLOR)[nodeId]:(nodeIdDim == SimplexMesh::topo::RIDGE)?(*BND_CURVE_COLOR)[nodeId]:(nodeIdDim == SimplexMesh::topo::SURFACE)?(*BND_SURFACE_COLOR)[nodeId]:0;

        if(nodeIdDim == SimplexMesh::topo::VOLUME || nodeIdDim == SimplexMesh::topo::CORNER || nodeIdDim == SimplexMesh::topo::RIDGE)
        {
          continue;
        }
        //std::cout << "nodeIdDim -> " << nodeIdDim << std::endl;

        //computing the optimal position of P in order to have Lm(PPj) ~ 1
        //https://tel.archives-ouvertes.fr/tel-00120327/document
        Metric<Eigen::Matrix3d> M0 = Metric<Eigen::Matrix3d>((*metric)[nodeId]);
        const SimplicesNode   node = SimplicesNode(m_simplexMesh, nodeId);
        const math::Point        p = node.getCoords();
        const std::vector<TSimplexID> ball = node.ballOf();
        std::vector<math::Point> P_opti{};
        if(ball.size() == 0){continue;}
        for(auto const simplex : ball)
        {
          if(simplex >= 0)
          {
            if(meshTet[simplex] != 0)
            {
              const SimplicesCell cell = SimplicesCell(m_simplexMesh, simplex);
              std::vector<TInt> nodes = cell.getNodes();
              math::Point p_opti(0.0, 0.0, 0.0);
              double factor = 0.0;
              for(auto const node : nodes)
              {
                if(node != nodeId)
                {
                  if(meshNode[node] != 0)
                  {
                    unsigned int nodeDim = ((*BND_VERTEX_COLOR)[node] != 0)?SimplexMesh::topo::CORNER:((*BND_CURVE_COLOR)[node] != 0)?SimplexMesh::topo::RIDGE:((*BND_SURFACE_COLOR)[node] != 0)?SimplexMesh::topo::SURFACE:SimplexMesh::topo::VOLUME;
                    unsigned int nodeLabel = (nodeDim == SimplexMesh::topo::CORNER)?(*BND_VERTEX_COLOR)[node]:(nodeDim == SimplexMesh::topo::RIDGE)?(*BND_CURVE_COLOR)[node]:(nodeDim == SimplexMesh::topo::SURFACE)?(*BND_SURFACE_COLOR)[node]:0;
                    if(nodeDim <= nodeIdDim)
                    {
                      if(nodeDim == nodeIdDim)
                      {
                        if(nodeLabel != nodeIdLabel)
                        {
                          continue;
                        }
                      }
                      if(nodeDim <=  nodeIdDim)
                      {
                        //std::cout << "node -> " << node << std::endl;
                        //std::cout << "nodeDim -> " << nodeDim << std::endl;
                        //std::cout <<std::endl;
                        const SimplicesNode sNode   = SimplicesNode(m_simplexMesh, node);
                        const math::Point pj        = sNode.getCoords();
                        Metric<Eigen::Matrix3d> M1  = Metric<Eigen::Matrix3d>((*metric)[node]);
                        double metricLenght         = M0.metricDist(p, pj, M1);
                        if(metricLenght != 0.0)
                        {
                          math::Vector3d vec = sNode.getCoords() - p;
                          p_opti = p_opti + (p + (vec / metricLenght));
                          factor += 1.0;
                        }
                        else
                        {
                          //DO SOMETHING ?
                        }
                      }
                    }
                  }
                  else
                  {
                    throw gmds::GMDSException("meshNode[node] == 0");
                  }
                }
              }
              if(factor != 0 )
              {
                P_opti.push_back(1.0 / factor * p_opti);
              }
            }
            else
            {
              throw gmds::GMDSException("meshTet[simplex] == 0");
            }
          }
        }


        math::Point new_P;
        if(P_opti.size() != 0)
        {
          for(auto const p : P_opti)
          {
            //std::cout << "p_opti -> " << p << std::endl;
            new_P = new_P + p;
          }
          new_P = math::Point(new_P.X() / (double)P_opti.size(), new_P.Y() / (double)P_opti.size(), new_P.Z() / (double)P_opti.size());
        }
        else
        {
          throw gmds::GMDSException("P_opti.size() == 0");
        }

        bool status = false;
        const gmds::BitVector markedNodes{};
        std::vector<TSimplexID> deletedSimplex{};
        std::vector<TInt> deletedNodes{};
        const std::multimap<TInt, TInt>& facesAlreadyBuilt{};
        std::vector<TSimplexID> cellsCreated{};
        bool alreadyAdd = false;
        std::vector<TSimplexID> tetraContenaingPt{};
        TInt new_Node = m_simplexMesh->addNodeAndcheck(new_P, tetraContenaingPt, alreadyAdd);
        if(nodeIdDim == SimplexMesh::topo::CORNER){BND_VERTEX_COLOR->set(new_Node, (*BND_VERTEX_COLOR)[nodeId]);}
        else if(nodeIdDim == SimplexMesh::topo::RIDGE){BND_CURVE_COLOR->set(new_Node, (*BND_CURVE_COLOR)[nodeId]);}
        else if(nodeIdDim == SimplexMesh::topo::SURFACE){BND_SURFACE_COLOR->set(new_Node, (*BND_SURFACE_COLOR)[nodeId]);}
        //std::cout << "new_Node -> " << new_Node << std::endl;
        //std::cout << "new_P -> " << new_P << std::endl;


        if(nodesAdded.capacity() != meshNode.capacity())
        {
          nodesAdded.resize(meshNode.capacity());
        }
        nodesAdded.assign(new_Node);
        if(!alreadyAdd)
        {
          PointInsertion pi(m_simplexMesh, SimplicesNode(m_simplexMesh, new_Node), criterionRAIS, status, ball, markedNodes, deletedNodes, facesAlreadyBuilt, cellsCreated);
          if(status)
          {
            //constant metric for now so no need to use interpolation metric
            metric->set(new_Node, (*metric)[nodeId]);
            m_simplexMesh->deleteNode(nodeId, true);
          }
        }
      }
    }*/
    std::cout << "LOOP ending -> " << cpt << std::endl;
    //cpt++;

  }
  std::cout << "meshNode.capacity() END -> " << meshNode.capacity() << std::endl;

}
