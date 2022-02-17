/*----------------------------------------------------------------------------*/
#ifndef SIMPLEX_MESH_H_
#define SIMPLEX_MESH_H_
/*----------------------------------------------------------------------------*/
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
#include <gmds/math/Matrix.h>
#include <gmds/utils/BitVector.h>
#include <gmds/utils/VariableManager.h>
#include <gmds/utils/Variable.h>
#include <gmds/utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
#include <gmds/hybridMeshAdapt/CommonInfo.h>
#include <gmds/hybridMeshAdapt/SimplicesNode.h>
#include <gmds/hybridMeshAdapt/SimplicesCell.h>
#include <gmds/hybridMeshAdapt/SimplicesTriangle.h>
#include <gmds/hybridMeshAdapt/StructuredGrid.h>
#include <gmds/hybridMeshAdapt/DelaunayPointInsertion.h>
#include <gmds/hybridMeshAdapt/PointInsertion.h>
#include <gmds/hybridMeshAdapt/CavityOperator.h>
#include <gmds/hybridMeshAdapt/Octree.h>
/*----------------------------------------------------------------------------*/
#include <gmds/math/Orientation.h>
#include <gmds/math/Plane.h>
#include <gmds/math/Ray.h>
/*----------------------------------------------------------------------------*/
#include <type_traits>
#include <map>
#include <set>
#include <bitset>
#include <iostream>
#include <fstream>
#include <math.h>
/*----------------------------------------------------------------------------*/
/** \class  SimplexMesh
 *  \brief  ???
 */
 namespace gmds{
   namespace hybrid{

     class Octree;
class SimplexMesh
{
 public:

   ///////////////
   ///////////////
   ///////////////
   enum topo
   {
     CORNER  = 0,
     RIDGE   = 1,
     SURFACE = 2,
     VOLUME  = 3
   };

   //To move on class somewhere else
   struct dataNode{
     TInt node;
     double lenght;
     int dim_Nj = 0;
     int index_Nj = 0;
     std::vector<TSimplexID> cavity{};
     std::vector<TInt> nodesInCavity{};
     TInt nodeToReinsert = std::numeric_limits<int>::min();
   };

   struct closestSimplex{
     TSimplexID closeSimplex;
     double minLenghtPtToNodeOfSimplex = std::numeric_limits<double>::max();
     double u = 0;
     double v = 0;
     double w = 0;
     double t = 0;
   };
   ///////////////
   ///////////////
   ///////////////
   friend class simplicesNode::SimplicesNode;
   friend class simplicesCell::SimplicesCell;
   friend class simplicesTriangle::SimplicesTriangle;
  /*------------------------------------------------------------------------*/
  /** \brief Constructor.
   */
  SimplexMesh();

  /*
   * Destructor
   */
  ~SimplexMesh();
  /*------------------------------------------------------------------------*/
  SimplexMesh& operator=(const SimplexMesh& simplexMesh);
  /*------------------------------------------------------------------------*/
  SimplexMesh(const SimplexMesh& simplexMesh);
  /*------------------------------------------------------------------------*/
  SimplexMesh(SimplexMesh&& simplexMesh);
  /*------------------------------------------------------------------------*/
  /** \brief get the number of ADim simplex
    */
  TInt getNbSimplex               ();

  /*return how many nodes the simplex mesh have*/
  TInt getNbNodes                 ();

  /*return how many tetra the simplex mesh have*/
  TInt getNbTetra                 ();

  /*return how many triangle the simplex mesh have*/
  TInt getNbTriangle              ();

  /*return the simplex containing the pt, if not return -1*/
  bool simplexContaining(const gmds::math::Point& pt, TSimplexID& tetraContainingPt);

  /*return the simplex containing the pt, if not return -1*/
  bool simplexContaining(const simplicesNode::SimplicesNode& node, TSimplexID& tetraContainingNode);

  bool checkSimplexContenaing(const gmds::math::Point& pt, TSimplexID& tetraContainingPt);

  bool checkSimplicesContenaing(const gmds::math::Point& pt, std::vector<TSimplexID>& tetraContainingPt, TSimplexID simplexToCheckFirst = std::numeric_limits<TSimplexID>::min());

  TSimplexID nextSimplexToCheck(const TSimplexID currentSimplex, const math::Point& pt, double& u, double& v, double& w, double& t, BitVector& cyclingCheck, closestSimplex& closerSimplexInfo);

  TSimplexID nextSimplexToCheckOrientation(const TSimplexID currentSimplex, const math::Point& pt, std::vector<math::Orientation::Sign>& uvwt, BitVector& cyclingCheck);

  /*return true if pt is in the tetra.. if not false*/
  bool isInSimplex(TSimplexID& tetra, const gmds::math::Point& pt);

  /*call reorientTet() for every Tetra*/
  void reorientAllTetra           ();

  /*return the opposite face of a triangle by a node*/
  TSimplexID buildOppositeFacesVector(const TSimplexID currentTriangle, const TInt currentNode);


  void reorientTetra(const TSimplexID & tetIndx); /*voir si on peut reorienté que le tetra crée*/

  void buildTetBaseAndAdjLocal(const TSimplexID & tetIndx);

  void buildTriBaseAndAdjLocal(const TSimplexID & triIndx);

  void buildAdjInfoGlobal();


  void buildBaseLocal(const TSimplexID& tetIndx);


  void buildOppFaces(const TSimplexID triIdx);

  /*return the idx of the bitset contenaing the first bit to 1*/
  TSimplexID firstBitTo1(const std::vector<bool>& vec );
  /*return true if at least 3 Node of array1 are in array2*/
  template <size_t N, size_t M>
  bool containNodes(const TSimplexID* array1, const TSimplexID* array2);

  /*compare 2 std::set and erase the componenents of set1 that are already in set2 (but not the last)*/
  template<typename T>
  void setComparator(std::vector<T>& set1, const std::vector<T>& set2);

  /* Return the opposite cell of  ATetID by the noed ANodeID --> if the return is <0 the opposite cell is a triangle index = -index*/
  TSimplexID getOppositeCell      (const TInt ANodeID, const TSimplexID ATetID);

  /* Return the triangle adjacente face of ATriID by the noed ANodeID*/
  TSimplexID getOppositeFace      (const TInt ANodeID, const TSimplexID ATriID);

  TSimplexID getSimplexFromBase   (const TInt ANodeID);

  /*adding some point to the mesh*/
  TInt addNode(const math::Point&& pt);

  /*adding some point to the mesh*/
  TInt addNode(const math::Point& pt);

  /*adding some point to the mesh*/
  TInt addNode(const TCoord X, const TCoord Y, const TCoord Z);


  /*adding some point to the mesh and check if  this point already exist*/
  TInt addNodeAndcheck(const math::Point& pt, std::vector<TSimplexID>& tetraContainingPt, bool& alreadyAdd, TSimplexID checkSimplicesContenaing = std::numeric_limits<TSimplexID>::min());

  /*return true if the node deleting passed well*/
  bool deleteNode(const TInt indexNode, bool eraseNode = true);

  /*return true if the node deleting passed well*/
  bool deleteNode(const simplicesNode::SimplicesNode& simpliceNode, bool eraseNode = true);


  /*adding some triangle to the mesh with existing SimplicesNode in the mesh*/
  TSimplexID addTriangle(const simplicesNode::SimplicesNode& ANode0,
                         const simplicesNode::SimplicesNode& ANode1,
                         const simplicesNode::SimplicesNode& ANode2,
                         bool flag = true);

  /*adding some triangle to the mesh with existing SimplicesNode in the mesh*/
  TSimplexID addTriangle(const TInt AIndexPoint0,
                         const TInt AIndexPoint1,
                         const TInt AIndexPoint2,
                         bool flag = true);

  /*return true if the triangle wad deleting with sucess*/
  std::vector<TInt> deleteTriangle(const TInt ATriangleIndex);

  /*adding some tetrahedre to the mesh with existing SimplicesNode in the mesh*/
  TSimplexID addTetraedre(const simplicesNode::SimplicesNode&& ANode0,
                         const simplicesNode::SimplicesNode&& ANode1,
                         const simplicesNode::SimplicesNode&& ANode2,
                         const simplicesNode::SimplicesNode&& ANode3,
                         const bool rebuildAdjinfo = true);

                         /*adding some tetrahedre to the mesh with existing SimplicesNode in the mesh*/
  TSimplexID addTetraedre(const TInt AIndexPoint0,
                          const TInt AIndexPoint1,
                          const TInt AIndexPoint2,
                          const TInt AIndexPoint3,
                          const bool rebuildAdjinfo = true);
  void fillHexahedron(const TInt ANode0, const TInt ANode1, const TInt ANode2, const TInt ANode3,
                      const TInt ANode4, const TInt ANode5, const TInt ANode6, const TInt ANode7);

  const std::vector<std::vector<TInt>>& getHexadronData(){return m_hexahedronData;}

  void setHexadronData(std::vector<std::vector<TInt>>& hexahedronData){m_hexahedronData = hexahedronData;}


  void setMarkedTet(const gmds::BitVector& markedTet){m_markedTet = markedTet;}

  bool checkMesh();

  bool checkMeshLocal(const simplicesNode::SimplicesNode node);

  bool doCellExist(const TSimplexID simplex) const ;

  bool doNodeExist(const TInt node) const ;

  void buildSimplexHull() ;

  TInt findRemainTriangleIdx(const TInt tri, gmds::BitVector& cyclingCheck);

  SimplexMesh& buildRobustLayerMesh(const unsigned int nbrLayer);

  SimplexMesh& buildRobustLayerMeshOrderedNode(const unsigned int nbrLayer);

  SimplexMesh& buildRobustLayerMeshOrderedNode01(const unsigned int nbrLayer);

  std::vector<TInt> buildQuadFaceFromCoordPt0(const std::vector<math::Point>& nodes, std::vector<TSimplexID>& markedSimplex);

  void buildQuadFaceFromNode(const std::vector<TInt>& nodesFace, std::vector<TSimplexID>& markedSimplex);

  std::vector<TInt> buildQuadFaceFromCoordPt(const std::vector<math::Point>& nodes, std::vector<TSimplexID>& markedSimplex);

  SimplexMesh& buildPrismFromFace0(const std::vector<TInt> & nodesDown, const std::vector<math::Point> & nodesUpCoord, std::vector<TSimplexID>& markedSimplex, const int debug = 0);

  SimplexMesh& buildPrismFromFace(const std::vector<TInt> & nodesDown, const std::vector<math::Point> & nodesUpCoord, std::vector<TSimplexID>& markedSimplex, const int debug = 0);

  SimplexMesh& buildQuadFromNodes0(const std::vector<math::Point>& nodes, std::vector<TSimplexID>& marquedSimplex);

  SimplexMesh& buildQuadFromNodes(const std::vector<math::Point>& nodes, std::vector<TSimplexID>& marquedSimplex);

  SimplexMesh& builHexaFromNodesAndSpline(const std::vector<math::Point>& nodes, std::vector<TSimplexID>& marquedSimplex);

  /*build intersection of the vector of tetra in parameter*/
  std::vector<TSimplexID> intersectionSimplex(std::vector<std::vector<TSimplexID>> & vectorOfBalls);

  gmds::BitVector FindTetInHexa(const std::vector<std::vector<TInt>>& nodesQuads) const;

  /*return true if the tetra has been deleted*/
  std::vector<TInt> deleteTetra(const TInt ATetraIndex);

  /*return the index of the next node*/
  TInt nextNode();

  /*initialize the nodeIndx*/
  TInt getFirstNode();

  /*return the m_node_ids capacity*/
  TInt nodesCapacity() {return m_node_ids.capacity();};

  /*return the index of the next Tet*/
  TInt nextTet();

  /*initialize the tetIndx*/
  TInt getFirstTet();

  /*return the m_node_ids capacity*/
  TInt tetCapacity() {return m_tet_ids.capacity();};

  /*return the index of the node*/
  TInt nextTri();

  /*initialize the triIndx*/
  TInt getFirstTri();

  /*return the m_node_ids capacity*/
  TInt triCapacity() {return m_tri_ids.capacity();};

  const gmds::BitVector& getBitVectorNodes() const {return m_node_ids;}

  const gmds::BitVector& getBitVectorTet() const {return m_tet_ids;}

  const gmds::BitVector& getMarkedTet() const {return m_markedTet;}

  const gmds::BitVector& getBitVectorTri() const {return m_tri_ids;}

  const std::map<unsigned int, std::pair<unsigned int, unsigned int>>& getEdgeTianglesIndices() const {return edgeTianglesIndices;}

  /*return a vector of simplex in border*/
  std::set<TSimplexID> simplexInBorder();

  /*delete All the tetra but tetra in arg*/
  void deleteAllSimplicesBut(const std::vector<TSimplexID> & simplices);

  void deleteAllTrianglesBut(const std::vector<TSimplexID> & triangles);

  void deleteAllTriangle();

  void clear();

  void buildHexaWithPattern(const math::Vector3i& nodeIdx, const std::unordered_map<math::Vector3i, TInt, StructuredGrid::hash_function>& uMapNodeIdx);

  template<typename T>
  void intersection(std::vector<std::vector<T>>& buffer0, const std::vector<T>& buffer1) const ;

  template<typename T>
  void intersection(std::vector<T>& buffer0, const std::vector<T>& buffer1) const ;

  void buildQuadEdgeFromNodes(const std::vector<TInt>& nodes, const operators::CriterionRAIS& criterion);

  void buildHexaFaceFromNodes(const std::vector<TInt>& nodes, const operators::CriterionRAIS& criterion);

  void buildHexaEdgeFromNodes(const std::vector<TInt>& nodes, const operators::CriterionRAIS& criterion);

  void buildFacesFromEdges(const std::vector<TInt>& nodes, const operators::CriterionRAIS& criterion);

  void buildHexaedre(const std::vector<TInt>& nodes, const operators::CriterionRAIS& criterion);

  bool MollerTriangleIntersection(const simplicesNode::SimplicesNode& nodeA, const math::Vector3d& dir, const std::vector<TInt>& triangleNodeId, math::Vector3d& tuv) ;

  std::vector<TSimplexID> cavityIntersectedByEdge(simplicesNode::SimplicesNode& simpliceNodeA, simplicesNode::SimplicesNode& simpliceNodeB);

  void hexaBuildPerformance(const std::vector<std::vector<TInt>>& orderedNodesQuads, double& edgePerformance, double& hexaPerformance) ;

  bool isFaceBuild(const std::vector<TInt>& nodes);

  void rebuildCavity(operators::CavityOperator::CavityIO& cavityIO, const TInt nodeToConnect);

  void pointsLabeling(const std::vector<math::Point> &points, std::vector<int>& pointsLabeling, std::vector<int>& topoIndex, std::vector<TInt>& nearNodes);

  void fillBNDVariable();

  unsigned int edgesRemove(const gmds::BitVector& nodeBitVector, std::vector<TSimplexID>& deletedNodes);



  unsigned int buildEdges(const std::multimap<TInt, TInt>& AEdges, const gmds::BitVector& nodeBitVector);

  bool isHexaEdgeBuild(const std::vector<std::vector<TInt>>& ANodesFaces);

  void whatFaceIsBuilt(const std::vector<std::vector<TInt>>& ANodesFaces, std::multimap<TInt, TInt>& facesAlreadyBuilt);

  std::vector<TSimplexID> hex2tet(const std::vector<TInt>& ANodesHex);

  std::vector<TSimplexID> initializeCavityWith(const TInt nodeA, const TInt nodeB);

  bool buildFace(const std::vector<TInt>& nodes, const gmds::BitVector& nodeAdded, const std::multimap<TInt, TInt>& facesAlreadyBuilt);

  bool pointInTriangle(const math::Point& query_point,
                       const math::Point& triangle_vertex_0,
                       const math::Point& triangle_vertex_1,
                       const math::Point& triangle_vertex_2,
                       double& distance,
                       math::Point& projectedPoint);

  /*adding some properties to the diffferent data of the mesh*/
  template<typename T, class C>
  gmds::Variable<T>* newVariable(const std::string& AName);

  /*getting  properties of the diffferent data of the mesh*/
  template<typename T, class C>
  gmds::Variable<T>* getVariable(const std::string&& AName);

  std::vector<VariableItf*> getAllVariables(ECellType AType) const;

  friend std::ostream&  operator<<(std::ostream& os, SimplexMesh& simplexMesh)
  {

    for(unsigned int idx = 0; idx < simplexMesh.m_tet_ids.capacity() ; idx++)
    {
      if(simplexMesh.m_tet_ids[idx] != 0)
      {
        const TInt Node0 = simplicesCell::SimplicesCell(&simplexMesh, idx).getNode(0).getGlobalNode();
        const TInt Node1 = simplicesCell::SimplicesCell(&simplexMesh, idx).getNode(1).getGlobalNode();
        const TInt Node2 = simplicesCell::SimplicesCell(&simplexMesh, idx).getNode(2).getGlobalNode();
        const TInt Node3 = simplicesCell::SimplicesCell(&simplexMesh, idx).getNode(3).getGlobalNode();

        os << "cellId = " << idx << "| Node : " << Node0 << " " << Node1 << " " << Node2 << " " << Node3 << std::endl;

      }
    }
    return os;
  }

  void setOctree(Octree* octree){m_octree = octree;}

  Octree* getOctree(){return m_octree;}

  double subSurfaceFactor(const std::vector<std::vector<TInt>>& faces);

private:

  //node coordinates with bitvector
  TInt nodeIndx = 0;
  gmds::BitVector m_node_ids; //a utiliser avec selectnewBit et assign //pour avoir l'index de la premiere case libre on utiliser la free_stach de bit vector
  std::vector<gmds::math::Point> m_coords;
  std::vector<TSimplexID> m_base;

  // Ordered node indices for each tet
  TInt tetIndx   = 0;
  TInt tetMarker  = 0; // for fast tetraedre building (peut mettre en atomic pour le parallelisme après)
  gmds::BitVector m_tet_ids;
  std::vector<std::vector<TInt>> m_tet_nodes;
  std::vector<std::vector<TInt>> m_tet_adj;

  TInt triIndx = 0;
  gmds::BitVector m_tri_ids;
  std::vector<std::vector<TInt>> m_tri_nodes; //3 sommets et 1 index du tetra voisin positif
  std::vector<std::vector<TInt>> m_tri_adj; //3 faces et 1 index du tetra voisin negatif

  Octree* m_octree = nullptr;
  std::map<unsigned int, std::pair<unsigned int, unsigned int>> edgeTianglesIndices{};

  //hexahedron data for tet extraction..
  std::vector<std::vector<TInt>> m_hexahedronData;

  gmds::BitVector m_markedTet;

  //VariableManager Node
  VariableManager* m_node_variable_manager     = nullptr;

  //VariableManager Cell (tetra for now)
  VariableManager* m_tet_variable_manager      = nullptr;

  //Variablemanager triangle
  VariableManager* m_tri_variable_manager = nullptr;
};

/******************************************************************************/
template<typename T>
void SimplexMesh::intersection(std::vector<std::vector<T>>& buffer0, const std::vector<T>& buffer1) const
{
  for(auto & simplices : buffer0)
  {
      simplices.erase(std::remove_if(simplices.begin(), simplices.end(), [&](T node)
      {
        bool flag = false;
        if(std::find(buffer1.begin(), buffer1.end(), node) == buffer1.end())
        {
          flag = true;
        };
        return flag;
      }), simplices.end());
  }
}
/******************************************************************************/
template<typename T>
void SimplexMesh::intersection(std::vector<T>& buffer0, const std::vector<T>& buffer1) const
{
  buffer0.erase(std::remove_if(buffer0.begin(), buffer0.end(), [&](T node)
                {
                  bool flag = false;
                  if(std::find(buffer1.begin(), buffer1.end(), node) == buffer1.end())
                  {
                    flag = true;
                  }
                  return flag;
                }
              ), buffer0.end());
}
/******************************************************************************/
template<typename T, class C>
gmds::Variable<T>* SimplexMesh::newVariable(const std::string& AName)
{
  gmds::Variable<T>* var = nullptr;

  if(std::is_same<C, simplicesNode::SimplicesNode>::value == true)
  {
    std::vector<TInt> ids;
    for(const auto & id : m_node_ids)
    {
      ids.push_back(id);
    }
    var = m_node_variable_manager->newVariable<T>(AName, m_node_ids.capacity(), &ids);
  }
  else if(std::is_same<C, simplicesCell::SimplicesCell>::value == true)
  {
    std::vector<TInt> ids;
    for(const auto & id : m_tet_ids)
    {
      ids.push_back(id);
    }
    var = m_tet_variable_manager->newVariable<T>(AName, m_tet_ids.capacity(), &ids);
  }
  else if(std::is_same<C,simplicesTriangle::SimplicesTriangle>::value == true)
  {
    std::vector<TInt> ids;
    for(const auto & id : m_tri_ids)
    {
      ids.push_back(id);
    }
    var = m_tri_variable_manager->newVariable<T>(AName, m_tri_ids.capacity(), &ids);
  }
  else
  {
    throw GMDSException("Unmanaged type of value -> impossible to create a variable");
  }

  return var;
}
/******************************************************************************/
template<typename T, class C>
gmds::Variable<T>* SimplexMesh::getVariable(const std::string&& AName)
{
  gmds::Variable<T>* var;

  if(std::is_same<C, simplicesNode::SimplicesNode>::value == true)
  {
    var = m_node_variable_manager->getVariable<T>(AName);
  }
  else if(std::is_same<C, simplicesCell::SimplicesCell>::value == true)
  {
    var = m_tet_variable_manager->getVariable<T>(AName);
  }
  else if(std::is_same<C,simplicesTriangle::SimplicesTriangle>::value == true)
  {
    var = m_tri_variable_manager->getVariable<T>(AName);
  }
  else
  {
    throw GMDSException("Unmanaged type of value -> impossible to access to a variable");
  }

  return var;
}


  }
}


/*----------------------------------------------------------------------------*/
#endif //A_H_
/*----------------------------------------------------------------------------*/
