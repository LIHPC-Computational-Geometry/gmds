/*----------------------------------------------------------------------------*/
/** \file    Cell.h
 *  \author  F. LEDOUX
 *  \date    January 6, 2014
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_CELL_H_
#define GMDS_CELL_H_
/*----------------------------------------------------------------------------*/
// gmds file headers
#include <gmds/utils/CommonTypes.h>
#include <gmds/utils/Exception.h>
#include <gmds/math/Point.h>
#include "GMDSIg_export.h"
/*----------------------------------------------------------------------------*/
// STL file headers
#include <vector>
#include <iostream>
/*----------------------------------------------------------------------------*/
namespace gmds{

  class Mesh;
  class Node;
  class Edge;
  class Face;
  class Region;
  /*----------------------------------------------------------------------------*/
  /** \class Cell
   *  \brief Defines the functions that are common to any type of cells (nodes,
   *  	   edges, faces, regions) in an Incidence Graph Representation.
   */
  /*----------------------------------------------------------------------------*/
  class GMDSIg_API Cell
  {
  public:
      /*------------------------------------------------------------------------*/
      /** \struct Data
       *  \brief nested structure to represent a cell by its dim and it id
       */
      /*------------------------------------------------------------------------*/
      struct Data{
          int dim;
          TCellID id;
          
          Data(const int ADim=-1, const TCellID AID=0):dim(ADim),id(AID){;}
      };
      
    /*------------------------------------------------------------------------*/
    /** \brief  Accessor to the global id. There exists a unique id for every
     * 			cell of a specific dimension but a face and an edge can have the
     * 			same id.
     *
     *  \return an id
     */
    TCellID id() const;
    /*------------------------------------------------------------------------*/
    /** \brief  Accesor to the cell dim.
     */
    virtual int dim() const = 0;
    /*------------------------------------------------------------------------*/
    /** \brief  Accesor to the cell type.
     */
    ECellType type() const;

    /*------------------------------------------------------------------------*/
    /** \brief Accessor th the number of incident nodes, edges, faces and
     * 		   adjacent regions
     */
    virtual TInt nbNodes()   const = 0;
    virtual TInt nbEdges()   const = 0;
    virtual TInt nbFaces()   const = 0;
    virtual TInt nbRegions() const = 0;

    /*------------------------------------------------------------------------*/
    /** \brief  Compute the center of the cell
     *
     * \return the center of the region
     */
    virtual math::Point center() const = 0;

    /*------------------------------------------------------------------------*/
    /** \brief  Compute the bounding box
     *
     * \param minXYZ the minimum corner of the bounding box
     * \param maxXYZ the maximum corner of the bounding box
     */
    void computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const;

    /*------------------------------------------------------------------------*/
    /** \brief  Accessor to the incident cells. Only the non-null cells are
     * 			provided.
     *
     * 			T can be Node, Edge, Face or Region
     */
    template<class T> std::vector<T> get() const;
    template<class T> void get(std::vector<T>&) const;

    /*------------------------------------------------------------------------*/
    /** \brief  Accessor to the ids of the incident cells.
     * 			Only the non-null cells are provided.
     *
     * 			T can be Node, Edge, Face or Region
     */
    template<class T> std::vector<TCellID> getIDs() const;
    template<class T> void getIDs(std::vector<TCellID>&) const;

    /*------------------------------------------------------------------------*/
    /** \brief  Accessor to the incident nodes including null cells. In this case
     * 			cells are provided in the storage order
     */
    template<class T> std::vector<T> getAll() const;
    template<class T> void getAll(std::vector<T>&) const;

    /*------------------------------------------------------------------------*/
    /** \brief  Accessor to the ids of incident nodes including null cells. In this case
     * 			cells are provided in the storage order
     */
    template<class T> std::vector<TCellID> getAllIDs() const;
    template<class T> void getAllIDs(std::vector<TCellID>&) const;


    /*------------------------------------------------------------------------*/
    /** \brief  update adjacency of type T to be equal to ACells.
     *  \param ACells cells to be added
     */
    template<class T> bool has(TCellID AId);
    template<class T> bool has(T& AElt);

  /*------------------------------------------------------------------------*/
    /** \brief  update adjacency of type T to be equal to ACells.
     *  \param ACells cells to be added
     */
    template<class T> void set(const std::vector<T>& ACells);


    template<class T> void set(const std::vector<TCellID>& ACells);
    /*------------------------------------------------------------------------*/
    /** \brief  update adjacency of type T to insert the new element AElt.
     *  \param AElt the element to be added
     */
    template<class T> void add(T& AElt);

    /*------------------------------------------------------------------------*/
    /** \brief  update adjacency of type T to insert the new element AElt.
     *  \param AElt the element to be added
     */
    template<class T> void add(TCellID AElt);

    /*------------------------------------------------------------------------*/
    /** \brief  update adjacency of type T to remove the element AElt.
     *          Does nothing if the element is not here
     *  \param AElt the element to be removed
     */
    template<class T>  void remove(T& AElt);

    /*------------------------------------------------------------------------*/
    /** \brief  update adjacency of type T to remove the element AElt.
     *  \param AElt the element to be removed
     */
    template<class T> void remove(TCellID AElt);

    /*------------------------------------------------------------------------*/
    /** \brief  update adjacency of type T to remove all the elements of type T
     */
    template<class T>  void removeAll();

    /*------------------------------------------------------------------------*/
    /** \brief  replace an incident T-typed cell AC1 by cell AC2 in the
     * 			incident elements
     *
     *  \param AC1 the cell to be replaced
     *  \param AC2 the new cell
     */
    template<class T>  void replace(T& AC1, T& AC2);

    /*------------------------------------------------------------------------*/
    /** \brief  replace the incident T-typed cell of ID AID11 by AID2 in the
     * 			incident elements
     *
     *  \param AID1 the id of the cell to be replaced
     *  \param AID2 the id of the new cell
     */
    template<class T>  void replace(TCellID AC1, TCellID AC2);

#ifdef GMDS_PARALLEL
    /*------------------------------------------------------------------------*/
    /** \brief  Gives the partition index of the mesh owner
     *
     *  \return Partition index
     */
    TInt getPartID() const;
#endif //GMDS_PARALLEL
  protected:

    /*------------------------------------------------------------------------*/
    /** \brief  Constructor.
     */
    Cell(Mesh* AMesh, const ECellType& AType, const TCellID& AID);

    /*------------------------------------------------------------------------*/
    /** \brief  Destructor.
     */
    virtual ~Cell(){ ; }
    /*------------------------------------------------------------------------*/
    /** \brief  Serializes the cell data into stream AStr
     *
     *  \param AStr an output stream where the cell data is written
     */
    void serializeCellData(std::ostream& AStr) const;
    /*------------------------------------------------------------------------*/
    /** \brief  Unserializes the cell data from stream AStr
     *
     *  \param AStr an input stream where the cell data is read from
     */
    void unserializeCellData(std::istream& AStr);

    /*------------------------------------------------------------------------*/
    /** \brief  Accessor to the incident cells. Only the non-null cells are
     * 			provided. T can be Node, Edge, Face or Region.
     */
    virtual void delegateGet(std::vector<Node>&   ACells) const = 0;
    virtual void delegateGet(std::vector<Edge>&   ACells) const = 0;
    virtual void delegateGet(std::vector<Face>&   ACells) const = 0;
    virtual void delegateGet(std::vector<Region>& ACells) const = 0;

    virtual void delegateGetNodeIDs(std::vector<TCellID>& ACells) const = 0;
    virtual void delegateGetEdgeIDs(std::vector<TCellID>& ACells) const = 0;
    virtual void delegateGetFaceIDs(std::vector<TCellID>& ACells) const = 0;
    virtual void delegateGetRegionIDs(std::vector<TCellID>& ACells) const = 0;

    virtual void delegateGetAll(std::vector<Node>&   ACells) const = 0;
    virtual void delegateGetAll(std::vector<Edge>&   ACells) const = 0;
    virtual void delegateGetAll(std::vector<Face>&   ACells) const = 0;
    virtual void delegateGetAll(std::vector<Region>& ACells) const = 0;

    virtual void delegateGetAllNodeIDs(std::vector<TCellID>& ACells) const = 0;
    virtual void delegateGetAllEdgeIDs(std::vector<TCellID>& ACells) const = 0;
    virtual void delegateGetAllFaceIDs(std::vector<TCellID>& ACells) const = 0;
    virtual void delegateGetAllRegionIDs(std::vector<TCellID>& ACells) const = 0;

    virtual void delegateSetNodeIDs(const std::vector<TCellID>& ACells) = 0;
    virtual void delegateSetEdgeIDs(const std::vector<TCellID>& ACells) = 0;
    virtual void delegateSetFaceIDs(const std::vector<TCellID>& ACells) = 0;
    virtual void delegateSetRegionIDs(const std::vector<TCellID>& ACells) = 0;


    virtual void delegateNodeAdd(TCellID AElt) = 0;
    virtual void delegateEdgeAdd(TCellID AElt) = 0;
    virtual void delegateFaceAdd(TCellID AElt) = 0;
    virtual void delegateRegionAdd(TCellID AElt) = 0;

    virtual void delegateNodeRemove(TCellID AElt) = 0;
    virtual void delegateEdgeRemove(TCellID AElt) = 0;
    virtual void delegateFaceRemove(TCellID AElt) = 0;
    virtual void delegateRegionRemove(TCellID AElt) = 0;

    virtual void delegateNodeReplace(TCellID AID1, TCellID AID2) = 0;
    virtual void delegateEdgeReplace(TCellID AID1, TCellID AID2) = 0;
    virtual void delegateFaceReplace(TCellID AID1, TCellID AID2) = 0;
    virtual void delegateRegionReplace(TCellID AID1, TCellID AID2) = 0;

    template<class T> std::vector<TCellID> convertCellToID(const std::vector<T>& ACells);
  protected:

    /** mesh containing *this*/
    Mesh* m_owner;

    /** cell type (node, edge, triangle, quad, polygon, tetrahedron, etc.)*/
    ECellType m_type;

    /** cell id locally to a part*/
    TCellID m_id;

    //
  };
  /*----------------------------------------------------------------------------*/
  // IMPLEMENTATION
  /*----------------------------------------------------------------------------*/
  template<class T> struct Type2Dim{};
  template<> struct Type2Dim<Node>  { static const int val = 0; };
  template<> struct Type2Dim<Edge>  { static const int val = 1; };
  template<> struct Type2Dim<Face>  { static const int val = 2; };
  template<> struct Type2Dim<Region>{ static const int val = 3; };
  /*----------------------------------------------------------------------------*/
  template<class T> std::vector<T> Cell::get() const{
    std::vector<T> vec;
    get(vec);
    return vec;
  }
  /*----------------------------------------------------------------------------*/
  template<class T> void Cell::get(std::vector<T>& ACells) const {
    delegateGet(ACells);
  }
  /*----------------------------------------------------------------------------*/
  template<class T> std::vector<TCellID> Cell::getIDs() const{
    std::vector<TCellID> vec;
    getIDs<T>(vec);
    return vec;
  }
  /*----------------------------------------------------------------------------*/
  template<class T> std::vector<T> Cell::getAll() const{
    std::vector<T> vec;
    getAll(vec);
    return vec;
  }
  /*----------------------------------------------------------------------------*/
  template<class T> void Cell::getAll(std::vector<T>& ACells) const {
    delegateGetAll(ACells);
  }
  /*----------------------------------------------------------------------------*/
  template<class T> std::vector<TCellID> Cell::getAllIDs() const{
    std::vector<TCellID> vec;
    getAllIDs<T>(vec);
    return vec;
  }
  /*----------------------------------------------------------------------------*/
  template<class T> void Cell::add(T& AElt) {
    add<T>(AElt.id());
  }
  /*----------------------------------------------------------------------------*/
  template<class T> bool Cell::has(T& AElt) {
    return has<T>(AElt.id());
  }
  /*----------------------------------------------------------------------------*/
  template<class T> bool Cell::has(TCellID AElt) {
    std::vector<TCellID> cellsIDs;
    getIDs<T>(cellsIDs);
    for(unsigned int iCell=0; iCell<cellsIDs.size(); iCell++) {
      if(cellsIDs[iCell] == AElt) {
        return true;
      }
    }
    return false;
  }
  /*----------------------------------------------------------------------------*/
  template<class T> void Cell::remove(T& AElt) {
    remove<T>(AElt.id());
  }
  /*----------------------------------------------------------------------------*/
  template<class T> void Cell::removeAll() {
    std::vector<TCellID> cellsIDs;
    getIDs<T>(cellsIDs);
    for(unsigned int iCell=0; iCell<cellsIDs.size(); iCell++) {
      remove<T>(cellsIDs[iCell]);
    }
  }
  /*----------------------------------------------------------------------------*/
  template<class T> void Cell::replace(T& AElt1, T& AElt2) {
    replace<T>(AElt1.id(), AElt2.id());
  }
  /*----------------------------------------------------------------------------*/
  template<class T> std::vector<TCellID>
    Cell::convertCellToID(const std::vector<T>& ACells)
    {
      std::vector<TCellID> cellIDs;
      cellIDs.resize(ACells.size());
      for (unsigned int i = 0; i < ACells.size(); i++){
	cellIDs[i] = ACells[i].id();
      }
      return cellIDs;
    }

  template<> GMDSIg_API void Cell::set<Node>(const std::vector<Node>& ACells);
  template<> GMDSIg_API void Cell::set<Edge>(const std::vector<Edge>& ACells);
  template<> GMDSIg_API void Cell::set<Face>(const std::vector<Face>& ACells);
  template<> GMDSIg_API void Cell::set<Region>(const std::vector<Region>& ACells);

  template<> GMDSIg_API void Cell::getIDs<Node>(std::vector<TCellID>& ACells) const;
  template<> GMDSIg_API void Cell::getIDs<Edge>(std::vector<TCellID>& ACells) const;
  template<> GMDSIg_API void Cell::getIDs<Face>(std::vector<TCellID>& ACells) const;
  template<> GMDSIg_API void Cell::getIDs<Region>(std::vector<TCellID>& ACells) const;

  template<> GMDSIg_API void Cell::getAllIDs<Node>(std::vector<TCellID>& ACells) const;
  template<> GMDSIg_API void Cell::getAllIDs<Edge>(std::vector<TCellID>& ACells) const;
  template<> GMDSIg_API void Cell::getAllIDs<Face>(std::vector<TCellID>& ACells) const;
  template<> GMDSIg_API void Cell::getAllIDs<Region>(std::vector<TCellID>& ACells) const;

  template<> GMDSIg_API void Cell::add<Node>(TCellID AElt);
  template<> GMDSIg_API void Cell::add<Edge>(TCellID AElt);
  template<> GMDSIg_API void Cell::add<Face>(TCellID AElt);
  template<> GMDSIg_API void Cell::add<Region>(TCellID AElt);

  template<> GMDSIg_API void Cell::remove<Node>(TCellID AElt);
  template<> GMDSIg_API void Cell::remove<Edge>(TCellID AElt);
  template<> GMDSIg_API void Cell::remove<Face>(TCellID AElt);
  template<> GMDSIg_API void Cell::remove<Region>(TCellID AElt);

  template<> GMDSIg_API void Cell::replace<Node>(TCellID AID1, TCellID AID2);
  template<> GMDSIg_API void Cell::replace<Edge>(TCellID AID1, TCellID AID2);
  template<> GMDSIg_API void Cell::replace<Face>(TCellID AID1, TCellID AID2);
  template<> GMDSIg_API void Cell::replace<Region>(TCellID AID1, TCellID AID2);

  /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /*GMDS_CELL_H_*/
/*----------------------------------------------------------------------------*/

