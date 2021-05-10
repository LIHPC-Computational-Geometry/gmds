/*----------------------------------------------------------------------------*/
#ifndef SIMPLEX_MESH_H_
#define SIMPLEX_MESH_H_
/*----------------------------------------------------------------------------*/
#include <gmds/math/Point.h>
#include <gmds/math/Vector.h>
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
#include <gmds/hybridMeshAdapt/PointInsertion.h>
#include <gmds/hybridMeshAdapt/PointSmoothing.h>
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

class EXPORT_GMDS SimplexMesh
{
 public:


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
  /** \brief get the number of ADim simplex
   */
  SimplexMesh(gmds::BitVector m_node_ids,
              std::vector<gmds::math::Point> m_coords,
              std::vector<TSimplexID> m_base,
              gmds::BitVector m_tet_ids,
              std::vector<TSimplexID* > m_tet_nodes,
              std::vector<TSimplexID* > m_tet_adj,
              gmds::BitVector m_tri_ids,
              std::vector<TSimplexID* > m_tri_nodes,
              std::vector<TSimplexID* > m_tri_adj
            );
  /*------------------------------------------------------------------------*/
  SimplexMesh(
              std::vector<gmds::math::Point> coords,
              std::vector<TSimplexID* > tet_nodes,
              std::vector<TSimplexID* > tet_adj,
              std::vector<TSimplexID* > tri_nodes,
              std::vector<TSimplexID* > tri_adj);
  /*------------------------------------------------------------------------*/
  SimplexMesh(std::vector<gmds::math::Point> coords,
              std::vector<TSimplexID* > tet_nodes,
              std::vector<TSimplexID* > tri_nodes,
              const bool flag = false);
  /*------------------------------------------------------------------------*/
  SimplexMesh(const StructuredGrid& structuredGrid);
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
  bool simplexContaining(const gmds::math::Point& pt, /*std::vector<*/TSimplexID/*>*/& tetraContainingPt);

  /*return the simplex containing the pt, if not return -1*/
  bool simplexContaining(const simplicesNode::SimplicesNode& node, /*std::vector<*/TSimplexID/*>*/& tetraContainingNode);

  /*return true if pt is in the tetra.. if not false*/
  bool isInSimplex(TSimplexID& tetra, const gmds::math::Point& pt);

  /*build the base vector*/
  void buildBaseVector            ();

  /*build the Adj vector (m_tet_adj) for tetraedres*/
  void buildAdjTetVector          ();

  /*call reorientTet() for every Tetra*/
  void reorientAllTetra           ();

  /*rebuild all opposote face of triangles..*/
  void buildOppFacesVector        ();


  /*rebuild all adjacent triangle vector (and)..*/
  std::vector<TSimplexID> buildAdjTriVector          (const TSimplexID currentTriangle);

  /*return the opposite face of a triangle by a node*/
  TSimplexID buildOppositeFacesVector(const TSimplexID currentTriangle, const TInt currentNode);


  void reorientTetra(const TSimplexID & tetIndx); /*voir si on peut reorienté que le tetra crée*/

  void buildBaseAndAdjLocal(const TSimplexID & tetIndx, const simplicesNode::SimplicesNode& ANode0,
                                                        const simplicesNode::SimplicesNode& ANode1,
                                                        const simplicesNode::SimplicesNode& ANode2,
                                                        const simplicesNode::SimplicesNode& ANode3);

  void buildBaseLocal(const TSimplexID& tetIndx);


  void buildAdjTet(const TSimplexID & idx, const simplicesNode::SimplicesNode& ANode0,
                                           const simplicesNode::SimplicesNode& ANode1,
                                           const simplicesNode::SimplicesNode& ANode2,
                                           const simplicesNode::SimplicesNode& ANode3);

  void buildOppFaces(const simplicesNode::SimplicesNode& ANode0,
                     const simplicesNode::SimplicesNode& ANode1,
                     const simplicesNode::SimplicesNode& ANode2);

  /*search Adjacent Tetraedre of currentTet*/
  TSimplexID adjacentTet_or_tri(const TInt currentTet, const TInt ANode);

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

  /*adding some point to the mesh*/
  TInt addNode(const math::Point&& pt);

  /*adding some point to the mesh*/
  TInt addNode(const math::Point& pt);

  /*adding some point to the mesh*/
  TInt addNode(const TCoord X, const TCoord Y, const TCoord Z);

  /*return true if the node deleting passed well*/
  bool deleteNode(const TInt indexNode);

  /*return true if the node deleting passed well*/
  bool deleteNode(const simplicesNode::SimplicesNode& simpliceNode);


  /*adding some triangle to the mesh with existing SimplicesNode in the mesh*/
  TSimplexID addTriangle(const simplicesNode::SimplicesNode& ANode0,
                         const simplicesNode::SimplicesNode& ANode1,
                         const simplicesNode::SimplicesNode& ANode2);

  /*adding some triangle to the mesh with existing SimplicesNode in the mesh*/
  TSimplexID addTriangle(const TInt AIndexPoint0,
                         const TInt AIndexPoint1,
                         const TInt AIndexPoint2);

  /*return true if the triangle wad deleting with sucess*/
  std::vector<TInt> deleteTriangle(const TInt ATriangleIndex);

  /*adding some tetrahedre to the mesh with existing SimplicesNode in the mesh*/
  TSimplexID addTetraedre(const simplicesNode::SimplicesNode&& ANode0,
                         const simplicesNode::SimplicesNode&& ANode1,
                         const simplicesNode::SimplicesNode&& ANode2,
                         const simplicesNode::SimplicesNode&& ANode3);

                         /*adding some tetrahedre to the mesh with existing SimplicesNode in the mesh*/
  TSimplexID addTetraedre(const TInt AIndexPoint0,
                          const TInt AIndexPoint1,
                          const TInt AIndexPoint2,
                          const TInt AIndexPoint3);

  SimplexMesh& buildRobustLayerMesh(const unsigned int nbrLayer);

  /*build intersection of the vector of tetra in parameter*/
  std::vector<TSimplexID> intersectionSimplex(std::vector<std::vector<TSimplexID>> & vectorOfBalls);

  /*return true if the tetra has been deleted*/
  std::vector<TInt> deleteTetra(const TInt ATetraIndex);

  /*return the index of the next node*/
  TInt nextNode();

  /*initialize the nodeIndx*/
  TInt getFirstNode();

  /*return the m_node_ids capacity*/
  TInt nodeCapacity() {return m_node_ids.capacity();};

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

  /*return a vector of simplex in border*/
  std::set<TSimplexID> simplexInBorder();

  void clear();

  void buildHexaWithPattern(const math::Vector3i& nodeIdx, const std::unordered_map<math::Vector3i, TInt, StructuredGrid::hash_function>& uMapNodeIdx);

  void saveAtSimplexMeshFormat(const std::string && destination, const std::string&& name);

  void loadSimplexMeshFormat(const std::string && file);

  template<typename T>
  void intersection(std::vector<std::vector<T>>& buffer0, const std::vector<T>& buffer1);

  template<typename T>
  void intersection(std::vector<T>& buffer0, const std::vector<T>& buffer1);


  /*adding some properties to the diffferent data of the mesh*/
  template<typename T, class C>
  gmds::Variable<T>* newVariable(const std::string&& AName);

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
  std::vector<TSimplexID* > m_tet_nodes;
  // Ordered adj cell for each tet
  std::vector<TSimplexID* > m_tet_adj;

  TInt triIndx = 0;
  gmds::BitVector m_tri_ids;
  std::vector<TSimplexID* > m_tri_nodes; //3 sommets et 1 index du tetra voisin positif
  // Ordered adj triangle for each triangle
  std::vector<TSimplexID* > m_tri_adj; //3 faces et 1 index du tetra voisin negatif


  //VariableManager Node
  VariableManager* m_node_variable_manager     = nullptr;

  //VariableManager Cell (tetra for now)
  VariableManager* m_tet_variable_manager      = nullptr;

  //Variablemanager triangle
  VariableManager* m_triangle_variable_manager = nullptr;
};


template<typename T, class C>
gmds::Variable<T>* SimplexMesh::newVariable(const std::string&& AName)
{
  gmds::Variable<T>* var = nullptr;

  if(std::is_same<C, simplicesNode::SimplicesNode>::value == true)
  {
    /*std::vector<TInt> ids;
    for(auto const & id: m_node_ids)
    {
      ids.push_back(id);
    }*/
    var = m_node_variable_manager->newVariable<T>(AName, m_node_ids.capacity()/*, &ids*/);
  }
  else if(std::is_same<C, simplicesCell::SimplicesCell>::value == true)
  {
    var = m_tet_variable_manager->newVariable<T>(AName, m_tet_ids.capacity());
  }
  else if(std::is_same<C,simplicesTriangle::SimplicesTriangle>::value == true)
  {
    var = m_triangle_variable_manager->newVariable<T>(AName, m_tri_ids.capacity());
  }
  else
  {
    throw GMDSException("Unmanaged type of value -> impossible to create a variable");
  }

  return var;
}


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
    var = m_triangle_variable_manager->getVariable<T>(AName);
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
