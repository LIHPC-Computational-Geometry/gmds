/******************************************************************************/
#include <gmds/hybridMeshAdapt/EdgeInsertion.h>
/******************************************************************************/
using namespace gmds;
using namespace hybrid;
using namespace operators;
using namespace simplicesNode;
using namespace simplicesCell;
/******************************************************************************/
EdgeInsertion::EdgeInsertion(SimplexMesh* simplexMesh, const simplicesNode::SimplicesNode& simpliceNodeA, const simplicesNode::SimplicesNode& simpliceNodeB, const CriterionRAIS& criterionRAIS, std::vector<TInt> markedNodes)
{
  if(simplexMesh != nullptr || simpliceNodeA.shell(simpliceNodeB).size() != 0)
  {
    int border = std::numeric_limits<int>::min();
    math::Vector3d vec = simpliceNodeB.getCoords() - simpliceNodeA.getCoords();
    std::vector<TSimplexID> tetraContainingPt = simpliceNodeA.directSimplices(vec);
    std::cout << "tetraContainingPt.size() = " << tetraContainingPt.size() << std::endl;
    //PointInsertion pi(simplexMesh, simpliceNodeB, criterionRAIS, tetraContainingPt, {}, markedNodes);
  }
}
/******************************************************************************/
/*EdgeInsertion::EdgeInsertion(SimplexMesh* simplexMesh, const simplicesNode::SimplicesNode& simpliceNodeA, const simplicesNode::SimplicesNode& simpliceNodeB, const CriterionRAIS& criterionRAIS, std::vector<TSimplexID> markedSimplex)
{
  if(simplexMesh != nullptr)
  {
    int border = std::numeric_limits<int>::min();
    std::vector<TSimplexID> cavity{};
    //biilding the teporary node and storing them in the nodesIntermediaire vector..
    if(simpliceNodeA.shell(simpliceNodeB).size() == 0)
    {
      std::vector<TInt> nodes{};
      bool cavityRebuild = false;
      bool EdgeInsert = insertionEdge(simplexMesh, simpliceNodeA, simpliceNodeB, cavityRebuild);

      fillNodeAdj(simplexMesh, simpliceNodeA, simpliceNodeB, nodes);
      if(!EdgeInsert)
      {
        std::vector<TInt> nodesTmp = nodes;

        for(;;)
        {

          std::cout << "EdgeInsertion" << std::endl;
          std::cout << simpliceNodeA << std::endl;
          std::cout << simpliceNodeB << std::endl;
          std::cout << "nodes.size() = " << nodes.size() << std::endl;
          for(auto const node : nodes)
          {
            std::cout << node << std::endl;
          }
          std::cout << std::endl;

          std::vector<TInt> nodesToDelete{};
          for(auto const node : nodes)
          {
            std::cout << node << std::endl;
            cavityRebuild = false;
            if(SimplicesNode(simplexMesh, node).shell(simpliceNodeB).size() == 0)
            {
              insertionEdge(simplexMesh, SimplicesNode(simplexMesh, node), simpliceNodeB, cavityRebuild);
            }
            std::cout << "cavityRebuild : " << cavityRebuild << std::endl;
            nodesToDelete.push_back(node);
            if(cavityRebuild)
            {
              break;
            }
          }
          std::cout << "nodesToDelete.size() : " << nodesToDelete.size() << std::endl;
          for(auto const& nodeToDelete : nodesToDelete)
          {
            std::cout << "nodeToDelete : " << nodeToDelete << std::endl;
            nodes.erase(std::remove(nodes.begin(), nodes.end(), nodeToDelete), nodes.end());
          }
          cavityRebuild = false;
          EdgeInsert = insertionEdge(simplexMesh, simpliceNodeA, simpliceNodeB, cavityRebuild);
          std::cout << "cavityRebuild : " << cavityRebuild << std::endl;
          if(EdgeInsert)
          {
            return;
          }
          else
          {
            if(!cavityRebuild && nodes.size() == 0)
            {
              return;
            }
            if(cavityRebuild)
            {
              nodes.clear();
              std::cout << "fillNodeAdj" << std::endl;
              fillNodeAdj(simplexMesh, simpliceNodeA, simpliceNodeB, nodes);
            }

            if(nodes == nodesTmp)
            {
              std::cout << "return" << std::endl;
              return;
            }
            else
            {
              nodesTmp = nodes;
            }
          }
        }
      }
    }
  }
}*/
/******************************************************************************/
bool EdgeInsertion::insertionEdge(SimplexMesh* simplexMesh, const simplicesNode::SimplicesNode& simpliceNodeA, const simplicesNode::SimplicesNode& simpliceNodeB, bool& cavityRebuild)
{
  const TInt border = std::numeric_limits<int>::min();
  bool flag = false;
  const operators::CriterionRAIS criterion;
  if(simpliceNodeA.shell(simpliceNodeB).size() == 0)
  {
    const math::Vector3d dir =  simpliceNodeB.getCoords() - simpliceNodeA.getCoords();

      for(;;)
      {
        TSimplexID simplex    = simpliceNodeA.directSimplex(dir);
        if(simplex != border)
        {
          TSimplexID simplexAdj = SimplicesCell(simplexMesh,simplex).oppositeTetraIdx(simpliceNodeA);
          if(simplexAdj == border || simplex == border)
          {
            break;
          }

          std::vector<TSimplexID> cavity{simplex, simplexAdj};
          if(simpliceNodeA.checkStar(cavity, criterion))
          {
            cavityRebuild = true;
            //PointInsertion pi(simplexMesh, simpliceNodeA, criterion, cavity);
            if(simpliceNodeA.shell(simpliceNodeB).size() != 0)
            {
              flag = true;
            }
            break;

          }
          else
          {
            break;
          }
        }
        else
        {
          flag = false;
          break;
        }
      }
  }
  else
  {
    flag = true;
  }
  return flag;
}
/******************************************************************************/
 void EdgeInsertion::fillNodeAdj(SimplexMesh* simplexMesh, const simplicesNode::SimplicesNode& simpliceNodeA, const simplicesNode::SimplicesNode& simpliceNodeB, std::vector<TInt>& nodes)
{

  const math::Vector3d dir =  simpliceNodeB.getCoords() - simpliceNodeA.getCoords();
  const TInt border = std::numeric_limits<int>::min();
  TSimplexID simplex = simpliceNodeA.directSimplex(dir);

  if(simplex != border)
  {
    std::vector<TInt> currentNode{simpliceNodeA.getGlobalNode()};
    nodes = SimplicesCell(simplexMesh, simplex).getOtherNodeInSimplex(currentNode);
    if(std::find(nodes.begin(), nodes.end(), simpliceNodeB.getGlobalNode()) != nodes.end())
    {
      nodes.clear();
    }
  }

}
/******************************************************************************/
/*EdgeInsertion::EdgeInsertion(SimplexMesh* simplexMesh, const simplicesNode::SimplicesNode& simpliceNodeA, const simplicesNode::SimplicesNode& simpliceNodeB, const CriterionRAIS& criterionRAIS, std::vector<TSimplexID> markedSimplex, unsigned int iter)
{
  //std::cout << "EdgeInsertion" << std::endl;
  //std::cout << simpliceNodeA << std::endl;
  //std::cout << simpliceNodeB << std::endl;

  if(simplexMesh != nullptr)
  {
    SimplicesNode A_tmp = simpliceNodeA;
    std::vector<SimplicesNode> nodesIntermediaire{};
    bool simplexAdjContainingNodeB = false;

    //biilding the teporary node and storing them in the nodesIntermediaire vector..
    if(simpliceNodeA.shell(simpliceNodeB).size() == 0)
    {
      while(!simplexAdjContainingNodeB)
      {
        const math::Point& nodeACoord  = A_tmp.getCoords();
        const math::Point& nodeBCoord  = simpliceNodeB.getCoords();
        math::Vector3d normal =  nodeBCoord - nodeACoord;
        TSimplexID tetraContainingDirectNeighbor = A_tmp.directSimplex(normal);
        TSimplexID tetraAdj = SimplicesCell(simplexMesh, tetraContainingDirectNeighbor).oppositeTetraIdx(A_tmp);

        if(SimplicesCell(simplexMesh, tetraAdj).containNode(simpliceNodeB))
        {
          simplexAdjContainingNodeB = true;
        }

        std::vector<TSimplexID>&& faceIntersection = SimplicesCell(simplexMesh, tetraContainingDirectNeighbor).intersectionNodes(SimplicesCell(simplexMesh, tetraAdj));

        if(faceIntersection.size() == 3)
        {

          math::Vector3d tuv;
          bool intersectionFaceFlag = simplexMesh->MollerTriangleIntersection(A_tmp, normal, faceIntersection, tuv);

          if(intersectionFaceFlag)
          {

            math::Point nodeToInsert = A_tmp.getCoords() + tuv[0] * normal;
            std::vector<TSimplexID> initCavity{tetraContainingDirectNeighbor, tetraAdj};
            TInt nodeId = simplexMesh->addNode(nodeToInsert);

            PointInsertion pi(simplexMesh, SimplicesNode(simplexMesh, nodeId), criterionRAIS, initCavity);

            nodesIntermediaire.push_back(SimplicesNode(simplexMesh, nodeId));
            A_tmp = SimplicesNode(simplexMesh,nodeId);

          }
          else
          {
            //to_do
            //std::vector<TSimplexID> initCavity{tetraContainingDirectNeighbor, tetraAdj};
            //simplexMesh->deleteAllSimplicesBut(initCavity);
            return;
          }
        }
        else
        {
          std::cout << "faceIntersection.size() != 3" << std::endl;
        }

      }
    }


    //Now we do a Point Smoothing of the intermediate node
    SimplicesNode nodeFromWhereRebuild = simpliceNodeA;
    std::vector<SimplicesNode> persistanteNodes{};
    for(std::vector<SimplicesNode>::iterator it = nodesIntermediaire.begin() ; it != nodesIntermediaire.end() ; it++)
    {
      std::cout << "it  : " << it->getGlobalNode() << std::endl;
      std::vector<TSimplexID>&& initCavity = it->ballOf();
      std::map<TSimplexID, std::vector<std::vector<TInt>>> map;
      bool cavityCanBeRebuildFromA = nodeFromWhereRebuild.checkStar(initCavity, map);

      if(cavityCanBeRebuildFromA)
      {
        PointInsertion pi;//(simplexMesh, nodeFromWhereRebuild, criterionRAIS, initCavity);
        std::vector<TInt> nodeToDelete = pi.PointInsertionWithNodeDelete(simplexMesh, nodeFromWhereRebuild, criterionRAIS, initCavity);

        nodesIntermediaire.erase(std::remove_if(nodesIntermediaire.begin(), nodesIntermediaire.end(),[&](SimplicesNode node){
          bool flag = false;
          if(std::find(nodeToDelete.begin(), nodeToDelete.end(), node.getGlobalNode()) != nodeToDelete.end()){
            //flag = true;
          }
          return flag;
        }), nodesIntermediaire.end());
      }
      else
      {
        nodeFromWhereRebuild = *it;
        persistanteNodes.push_back(nodeFromWhereRebuild);
      }
      //break;
    }


    //()petite optimisation check de si la cavité est etoilable par B et non plus par A
    for(std::vector<SimplicesNode>::reverse_iterator rit = persistanteNodes.rbegin(); rit < persistanteNodes.rend(); rit++)
    {
      SimplicesNode persistanteNode = *rit;
      std::vector<TSimplexID>&& persistanteNodeBall = persistanteNode.ballOf();
      std::map<TSimplexID, std::vector<std::vector<TInt>>> map;
      bool cavityCanBeRebuildFromB = nodeFromWhereRebuild.checkStar(persistanteNodeBall, map);
      if(cavityCanBeRebuildFromB)
      {
        PointInsertion pi(simplexMesh, simpliceNodeB, criterionRAIS, persistanteNodeBall);
        persistanteNodes.erase(std::remove(persistanteNodes.begin(),persistanteNodes.end(),*rit), persistanteNodes.end());
      }
      else
      {
        break;
      }
    }


    //découpage des maillage pour rendre la cavité convexe..
    int errorId = std::numeric_limits<int>::min();
    for(auto const & persistanteNode : persistanteNodes)
    {
      std::map<TSimplexID, std::vector<std::vector<TInt>>> map;
      std::vector<TSimplexID>&& persistanteNodeBall = persistanteNode.ballOf();
      while(!simpliceNodeA.checkStar(persistanteNodeBall, map))
      {
        for(auto & it : map)
        {
          TSimplexID simplex0 = it.first;
          for(auto const & nodes : it.second)
          {
            for(auto const & node : nodes)
            {
              bool flagNodeToMove = true;
              SimplicesNode simpliceNode = SimplicesNode(simplexMesh, node);
              std::vector<TSimplexID>&& ballNode = simpliceNode.ballOf();
              for(auto const & tet: ballNode)
              {
                if( std::find(persistanteNodeBall.begin(), persistanteNodeBall.end(), tet) != persistanteNodeBall.end()
                   && SimplicesCell(simplexMesh, tet).containNode(simpliceNodeA))
                {
                  flagNodeToMove = false;
                  break;
                }
              }

              if(flagNodeToMove)
              {
                ballNode.erase(std::remove_if(ballNode.begin(), ballNode.end(), [&](TSimplexID simplex){
                  bool flag = false;
                  if(std::find(persistanteNodeBall.begin(), persistanteNodeBall.end(), simplex) == persistanteNodeBall.end())
                  {
                    flag = true;
                  }
                  return flag;
                }),ballNode.end());


                double t = 0.5;
                math::Vector3d dir = persistanteNode.getCoords() - simpliceNode.getCoords();
                math::Point nodeToInsert = simpliceNode.getCoords() + t * dir;
                TInt nodeToInsertIdx = simplexMesh->addNode(nodeToInsert);
                PointInsertion pi(simplexMesh, SimplicesNode(simplexMesh, nodeToInsertIdx), criterionRAIS, ballNode);
              }
            }
          }
        }

        persistanteNodeBall = persistanteNode.ballOf();
        map.clear();
      }
      std::vector<TSimplexID> PersitanteNodeBall = persistanteNode.ballOf();
      PointInsertion pi(simplexMesh, simpliceNodeA, criterionRAIS, PersitanteNodeBall);
    }
  }
  //std::cout << "Ending EdgeInsertion" << std::endl;
}*/
