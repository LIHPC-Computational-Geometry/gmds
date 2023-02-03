/*----------------------------------------------------------------------------*/
/* Mesh_def.h
 *
 *  Created on: 20 mai 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MESH_H_
#	define GMDS_MESH_H_
/*----------------------------------------------------------------------------*/
#	include <list>
/*----------------------------------------------------------------------------*/
#include <gmds/utils/CommonTypes.h>
#include <gmds/utils/BitVector.h>
#include <gmds/utils/IndexedVector.h>
#include <gmds/utils/Variable.h>
#include <gmds/utils/VariableManager.h>
#include <gmds/utils/Marks32.h>
#include <gmds/ig/NodeContainer.h>
#include <gmds/ig/EdgeContainer.h>
#include <gmds/ig/FaceContainer.h>
#include <gmds/ig/RegionContainer.h>
#include <gmds/ig/CellGroup.h>
#include "GMDSIg_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*------------------------------------------------------------------------*/

namespace cad {
class GeomEntity;
}
/*------------------------------------------------------------------------*/
/** \class Mesh
    *
    *   \brief  this class represents meshes as general incidence graphs. It is
    *           possible to select the cells and connectivities.
 */
class GMDSIg_API Mesh
{
	friend class IGMeshIOService;

 public:
	/*------------------------------------------------------------------------*/
	/** \brief friend class to access to protected/private data
	 */
	friend class Region;
	friend class Face;
	friend class Edge;
	friend class Node;

	/*------------------------------------------------------------------------*/
	/** \brief Definition of local iterators
	 */
	using nodes_iterator = NodeContainer::iterator;
	using edges_iterator = EdgeContainer::iterator;
	using faces_iterator = FaceContainer::iterator;
	using regions_iterator = RegionContainer::iterator;

	using node = Node;
	using edge = Edge;
	using face = Face;
	using region = Region;

	template<typename T> using group_iterator = typename std::list<CellGroup<T> *>::iterator;
	/*------------------------------------------------------------------------*/
	/** \brief Constructor
         * \param AModel a mesh model defining available cells and  connectivities
	 */
	explicit Mesh(MeshModel model);

 public:
	NodeContainer &nodes() const
	{
		return *m_nodes_container;
	}
	EdgeContainer &edges() const
	{
		return *m_edges_container;
	}
	FaceContainer &faces() const
	{
		return *m_faces_container;
	}
	RegionContainer &regions() const
	{
		return *m_regions_container;
	}
	/*------------------------------------------------------------------------*/
	/** \brief Destructor
	 */
	virtual ~Mesh();

	/*------------------------------------------------------------------------*/
	/** \brief Accessor to the mesh model
         *
         * \return the mesh model
	 */
	MeshModel getModel() const;

	/*------------------------------------------------------------------------*/
	/** \brief Remove all the cells, groups and variables stored in this mesh
	 */
	void clear();

	/*------------------------------------------------------------------------*/
	/** \brief Change the mesh model
         *
         * \param the mesh model
         * \param a boolean that chooses whether to create new entites/adjacencies
	 */
	void changeModel(const MeshModel &AModel, const bool &ACallDoctor = false);

	/*------------------------------------------------------------------------*/
	/** \brief Accessor to the mesh dimension
         *
         * \return the mesh dimension
	 */
	inline TInt getDim() const
	{
		return m_model.getDim();
	}

	/*------------------------------------------------------------------------*/
	/** \brief return the number of cells in the mesh (per dimension)
	 */
	inline TInt getNbNodes() const
	{
		return m_nodes_container->getNbElements();
	}
	inline TInt getNbEdges() const
	{
		return m_edges_container->getNbElements();
	}
	inline TInt getNbFaces() const
	{
		return m_faces_container->getNbElements();
	}
	inline TInt getNbRegions() const
	{
		return m_regions_container->getNbElements();
	}

	inline TInt getNbTriangles() const
	{
		return m_faces_container->m_T2N->size();
	}
	inline TInt getNbQuadrilaterals() const
	{
		return m_faces_container->m_Q2N->size();
	}
	inline TInt getNbTetrahedra() const
	{
		return m_regions_container->m_T2N->size();
	}
	inline TInt getNbHexahedra() const
	{
		return m_regions_container->m_H2N->size();
	}

	/*------------------------------------------------------------------------*/
	/** \brief return the max id for dimension ADim
         *
         * \param ADim is the dimension of cells we want to retrieve the max id
	 */
	TCellID getMaxLocalID(const TInt &ADim) const;

	/*------------------------------------------------------------------------*/
	/** \brief retrieve a cell from its id
         *
         * \param AID the id of the cell we look for
	 */
	template<typename T> GMDSIg_API T get(const TCellID &AID) const;

	/*------------------------------------------------------------------------*/
	/** \brief provdes a STL vector containing all the T-type cells of the mesh
         *
         * \param AVec the vector of cells
	 */
	template<typename T> GMDSIg_API void getAll(std::vector<T> &AVec) const;

	/*------------------------------------------------------------------------*/
	/** \brief Tell if there is a T-type cell having id AID
         *
         * \param AID the id of the cell we look for
	 */
	template<typename T> GMDSIg_API bool has(const TCellID &AID) const;

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to create a node
         *
         * \param AX X coordinate
         * \param AY Y coordinate
         * \param AZ Z coordinate
         *
         * \return a node object that encapsulates access to the mesh node
	 */
	Node newNode(const TCoord &AX = 0.0, const TCoord &AY = 0.0, TCoord AZ = 0.0);
	Node newNode(const math::Point &APnt);

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to create an edge
         *
         * \param AN1 Node 1
         * \param AN2 Node 2
         * \return an edgece object that encapsulates access to the mesh edge
	 */
	Edge newEdge(const Node &AN1, const Node &AN2);
	Edge newEdge(const TCellID &AN1, const TCellID &AN2);

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to create a triangle
         *
         * \param AN1 Node 1
         * \param AN2 Node 2
         * \param AN3 Node 3
         * \return a face object that encapsulates access to the mesh face
	 */
	Face newTriangle(const Node &AN1, const Node &AN2, const Node &AN3);
	Face newTriangle(const TCellID &AN1 = NullID, const TCellID &AN2 = NullID, const TCellID &AN3 = NullID);

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to create a quad
         *
         * \param AN1 Node 1
         * \param AN2 Node 2
         * \param AN3 Node 3
         * \param AN4 Node 4
         * \return a face object that encapsulates access to the mesh face
	 */
	Face newQuad(const Node &AN1, const Node &AN2, const Node &AN3, const Node &AN4);
	Face newQuad(const TCellID &AN1 = NullID, const TCellID &AN2 = NullID, const TCellID &AN3 = NullID, const TCellID &AN4 = NullID);
	/*------------------------------------------------------------------------*/
	/** \brief Factory method to create a polygon
         *
         * \param ANodes ordered vector of nodes
         * \return a face object that encapsulates access to the mesh face
	 */
	Face newPolygon(const std::vector<Node> &ANodes);
	Face newPolygon(const std::vector<TCellID> &ANodes);
	/*------------------------------------------------------------------------*/
	/** \brief Factory method to create a face
         *
         * \param ANodes ordered vector of nodes
         * \return a face object that encapsulates access to the mesh face
	 */
	Face newFace(const std::vector<Node> &ANodes);
	Face newFace(const std::vector<TCellID> &ANodes);
	/*------------------------------------------------------------------------*/
	/** \brief Factory method to create a tetrahedral element
         *
         * \param AN1 Node 1
         * \param AN2 Node 2
         * \param AN3 Node 3
         * \param AN4 Node 4
         * \return a region object that encapsulates access to the mesh region
	 */
	Region newTet(const Node &AN1, const Node &AN2, const Node &AN3, const Node &AN4);
	Region newTet(const TCellID &AN1 = NullID, const TCellID &AN2 = NullID, const TCellID &AN3 = NullID, const TCellID &AN4 = NullID);

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to create a hexahedral element
         *
         * \param AN1 Node 1
         * \param AN2 Node 2
         * \param AN3 Node 3
         * \param AN4 Node 4
         * \param AN5 Node 5
         * \param AN6 Node 6
         * \param AN7 Node 7
         * \param AN8 Node 8
         * \return a region object that encapsulates access to the mesh region
	 */
	Region newHex(const Node &AN1, const Node &AN2, const Node &AN3, const Node &AN4, const Node &AN5, const Node &AN6, const Node &AN7, const Node &AN8);

	Region newHex(const TCellID &AN1 = NullID,
	              const TCellID &AN2 = NullID,
	              const TCellID &AN3 = NullID,
	              const TCellID &AN4 = NullID,
	              const TCellID &AN5 = NullID,
	              const TCellID &AN6 = NullID,
	              const TCellID &AN7 = NullID,
	              const TCellID &AN8 = NullID);
	/*------------------------------------------------------------------------*/
	/** \brief Factory method to create a pyramid element. The last node is
         * 		the top of the pyramid
         *
         * \param AN1 Node 1
         * \param AN2 Node 2
         * \param AN3 Node 3
         * \param AN4 Node 4
         * \param AN5 Node 5
         * \return a region object that encapsulates access to the mesh region
	 */
	Region newPyramid(const Node &AN1, const Node &AN2, const Node &AN3, const Node &AN4, const Node &AN5);
	Region newPyramid(
	   const TCellID &AN1 = NullID, const TCellID &AN2 = NullID, const TCellID &AN3 = NullID, const TCellID &AN4 = NullID, const TCellID &AN5 = NullID);

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to create a prism3 element.
         * \param AN1 Node 1
         * \param AN2 Node 2
         * \param AN3 Node 3
         * \param AN4 Node 4
         * \param AN5 Node 5
         * \param AN6 Node 6
         * \return a region object that encapsulates access to the mesh region
	 */
	Region newPrism3(const Node &AN1, const Node &AN2, const Node &AN3, const Node &AN4, const Node &AN5, const Node &AN6);
	Region newPrism3(const TCellID &AN1 = NullID,
	                 const TCellID &AN2 = NullID,
	                 const TCellID &AN3 = NullID,
	                 const TCellID &AN4 = NullID,
	                 const TCellID &AN5 = NullID,
	                 const TCellID &AN6 = NullID);

	/*------------------------------------------------------------------------*/
	/** \brief Deletion of the node AN from the mesh
         *  \param AN ANode
	 */
	void deleteNode(const Node &AN);
	/*------------------------------------------------------------------------*/
	/** \brief Deletion of a node from the mesh
         *  \param AID A node id
	 */
	void deleteNode(TCellID n);
	/*------------------------------------------------------------------------*/
	/** \brief Deletion of the edge AE from the mesh
         * \param AE An edge
	 */
	void deleteEdge(const Edge &e);

	/*------------------------------------------------------------------------*/
	/** \brief Deletion of an edge
         *  \param AID An edge id
	 */
	void deleteEdge(TCellID e);
	/*------------------------------------------------------------------------*/
	/** \brief Deletion of a face
         * \param AF AFace
	 */
	void deleteFace(const Face &f);
	/*------------------------------------------------------------------------*/
	/** \brief Deletion of the face AF from the mesh
         * \param AF A face id
	 */
	void deleteFace(TCellID f);
	/*------------------------------------------------------------------------*/
	/** \brief Deletion of a region
         * \param AR A region
	 */
	void deleteRegion(const Region &AR);

	/*------------------------------------------------------------------------*/
	/** \brief Deletion of a region
         *  \param AR A region id
	 */
	void deleteRegion(TCellID AR);

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to get an iterator on the first node
         *
         * \return a node iterator
	 */
	nodes_iterator nodes_begin()
	{
		return m_nodes_container->begin();
	}

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to get an iterator on the last node
         *
         * \return a node iterator
	 */
	nodes_iterator nodes_end()
	{
		return m_nodes_container->end();
	}

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to get an iterator on the first edge
         *
         * \return an edge  iterator
	 */
	edges_iterator edges_begin()
	{
		return m_edges_container->begin();
	}

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to get an iterator on the last edge
         *
         * \return an edge iterator
	 */
	edges_iterator edges_end()
	{
		return m_edges_container->end();
	}
	/*------------------------------------------------------------------------*/
	/** \brief Factory method to get an iterator on the first face
         *
         * \return a face iterator
	 */
	faces_iterator faces_begin()
	{
		return m_faces_container->begin();
	}

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to get an iterator on the last face
         *
         * \return a face iterator
	 */
	faces_iterator faces_end()
	{
		return m_faces_container->end();
	}
	/*------------------------------------------------------------------------*/
	/** \brief Factory method to get an iterator on the first region
         *
         * \return a region iterator
	 */
	regions_iterator regions_begin()
	{
		return m_regions_container->begin();
	}

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to get an iterator on the last region
         *
         * \return a region iterator
	 */
	regions_iterator regions_end()
	{
		return m_regions_container->end();
	}

	/*------------------------------------------------------------------------*/
	/** \brief  Reserve a mark of the mesh for the cell type T
         *
         *  \return A mark number
	 */
	template<typename T> GMDSIg_API TInt newMark();

	/*------------------------------------------------------------------------*/
	/** \brief Get the marks of a cell
         *
         *  \return A mark
	 */
	template<typename T> GMDSIg_API Marks32 getMarks(const T &ACell) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Free mark AMarkNumber which was previously reserved with
         * 			getNewMark()
         *
         *  \param AMarkNumber  A mark number
	 */
	template<typename T> GMDSIg_API void freeMark(TInt AMarkNumber);

	/*------------------------------------------------------------------------*/
	/** \brief  Invert the mark value for all the cell in the mesh in O(1). This
         * 			is useful for unmark all the cells of a mesh at the end of an
         * 			algorithm
         *
         *  \param AMarkNumber A mark number
	 */
	template<typename T> GMDSIg_API void negateMaskMark(TInt AMarkNumber);

	/*------------------------------------------------------------------------*/
	/** \brief  Invert the mark value for all the cell in the mesh in O(n) where
         * 			n is the number of cells in the mesh. It is the only way to
         * 			ensure to have a uniform state for this mark.
         *
         *  \param AMarkNumber A mark number
	 */
	template<typename T> GMDSIg_API void unmarkAll(TInt AMarkNumber);

	/*------------------------------------------------------------------------*/
	/** \brief  Test if ACell is marked with mark AMarkNumber.
         *
         *  \param ACell		a cell
         *  \param AMarkNumber 	A mark number
         *
         *  \return A boolean providing the value of mark AMarkNumber for cell ACell
	 */
	bool isMarked(const Node &ACell, TInt AMarkNumber) const;
	bool isMarked(const Edge &ACell, TInt AMarkNumber) const;
	bool isMarked(const Face &ACell, TInt AMarkNumber) const;
	bool isMarked(const Region &ACell, TInt AMarkNumber) const;
	template<typename T> GMDSIg_API bool isMarked(const TCellID &ACellID, TInt AMarkNumber) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Update value of mark AMarkNumber for cell ACell.
         *
         *  \param ACell		A cell
         *  \param AMarkNumber 	A mark number
         *  \param AState		The new value of mark AMarkNumber for ACell
	 */
	void markTo(const Node &ACell, TInt AMarkNumber, bool AState);
	void markTo(const Edge &ACell, TInt AMarkNumber, bool AState);
	void markTo(const Face &ACell, TInt AMarkNumber, bool AState);
	void markTo(const Region &ACell, TInt AMarkNumber, bool AState);
	template<typename T> GMDSIg_API void markTo(const TCellID &ACellID, TInt AMarkNumber, bool AState);

	/*------------------------------------------------------------------------*/
	/** \brief  Mark cell ACell with mark AMarkNumber.
         * 			Equivalent to markTo(ADart, AMarkNumber, true).
         *
         *  \param ACell		A cell
         *  \param AMarkNumber 	A mark number
	 */
	void mark(const Node &ACell, TInt AMarkNumber);
	void mark(const Edge &ACell, TInt AMarkNumber);
	void mark(const Face &ACell, TInt AMarkNumber);
	void mark(const Region &ACell, TInt AMarkNumber);
	template<typename T> GMDSIg_API void mark(const TCellID &ACellID, TInt AMarkNumber);
	/*------------------------------------------------------------------------*/
	/** \brief  Unmark cell ACell with mark AMarkNumber.
         * 			Equivalent to markTo(ADart, AMarkNumber, false).
         *
         *  \param ACell 		A cell
         *  \param AMarkNumber 	A mark number
	 */
	void unmark(const Node &ACell, TInt AMarkNumber);
	void unmark(const Edge &ACell, TInt AMarkNumber);
	void unmark(const Face &ACell, TInt AMarkNumber);
	void unmark(const Region &ACell, TInt AMarkNumber);
	template<typename T> GMDSIg_API void unmark(const TCellID &ACellID, TInt AMarkNumber);

	/*------------------------------------------------------------------------*/
	/** \brief  Create a new variable attached to a generic cell type
         * 			(GMDS_NODE, GMDS_EDGE, GMDS_FACE or GNDS_REGION). Raises an
         * 			exception if a variable on AType has already name AName
         *
         *  \param AType the cell type
         *  \param AName the name of the variable. If this name is already used, the
         *  	   variable is not created and an exception is thrown
         *
         *  \return A pointer on the variable
	 */
	template<typename T, ECellType E> Variable<T> *newVariable(const std::string &AName);

	/*------------------------------------------------------------------------*/
	/** \brief  Returns whether the variable attached to a generic cell type
         * 			(GMDS_NODE, GMDS_EDGE, GMDS_FACE or GNDS_REGION) exists.
         *
         *  \param AType the cell type
         *  \param AName the name of the queried variable.
         *
         *  \return A boolean
	 */
	bool hasVariable(ECellType AType, const std::string &AName) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to variable attached to a generic cell type
	 * 			(GMDS_NODE, GMDS_EDGE, GMDS_FACE or GNDS_REGION).Raises an
	 * 			exception if no variable on AType has  name AName

	 *
	 *  \param AType the cell type
	 *  \param AName the name of the variable. If this name does not exist, the
	 *  	   an exception is thrown
	 *
	 *  \return A pointer on the variable. This pointer can be null if the
	 *  		specified type is wrong
	 */
	template<typename T, ECellType E> Variable<T> *getVariable(const std::string &AName) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to variable attached to a generic cell type
         * 			(GMDS_NODE, GMDS_EDGE, GMDS_FACE or GNDS_REGION) If it does not
         * 			exist it is created.
         *
         *  \param AType the cell type
         *  \param AName the name of the variable. If this name does not exist, the
         *  	   an exception is thrown
         *
         *  \return A pointer on the variable. This pointer can be null if the
         *  		specified type is wrong
	 */
	template<typename T, ECellType E> Variable<T> *getOrCreateVariable(const std::string &AName);

	/*------------------------------------------------------------------------*/
	/** \brief  Access to variable attached to a generic cell type
         * 			(GMDS_NODE, GMDS_EDGE, GMDS_FACE or GNDS_REGION)
         *
         *  \param AType the cell type
         *  \param AName the name of the variable. If this name does not exist, the
         *  	   an exception is thrown
         *
         *  \return A pointer on the variable. This pointer can be null if the
         *  		specified type is wrong
	 */
	std::vector<VariableItf *> getAllVariables(ECellType AType) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Delete a variable attached to a cell type (GMDS_NODE, GMDS_EDGE,
         * 			GMDS_FACE or GNDS_REGION)
         *
         *  \param AType the cell type
         *  \param AName the name of the variable to be deleted
	 */
	void deleteVariable(ECellType AType, const std::string &AName);

	/** \brief  Delete a variable attached to a cell type (GMDS_NODE, GMDS_EDGE,
         * 			GMDS_FACE or GNDS_REGION)
         *
         *  \param AType the cell type
         *  \param AVar a pointer on the variable to be deleted
	 */
	void deleteVariable(ECellType AType, VariableItf *AVar);

	/*------------------------------------------------------------------------*/
	/** \brief  Create a new cloud inside the mesh. This cloud allows
         * 			users to gather nodes having common properties. Be careful when
         * 			you remove a node from the mesh. If this cell belongs to a
         * 			surface, this cloud will keep an invalid pointer onto this
         * 			node.
         *
         * 	\param AName the name of the new cloud
         *
         *  \return a cloud
	 */
	template<typename T> CellGroup<T> *newGroup(const std::string &AName);

	/*------------------------------------------------------------------------*/
	/** \brief  Delete a cloud already available in the mesh.
         *
         * 	\param ACloud the cloud to delete
	 */
	template<typename T> void deleteGroup(CellGroup<T> *ACloud);

	/*------------------------------------------------------------------------*/
	/** \brief  return the group of T cells named AName if it exists,
         *  \throw  an exception otherwise.
         *
         * 	\param AName the name of the group of T cells
         *
         *  \return a cloud
	 */
	template<typename T> CellGroup<T> *getGroup(const std::string &AName);

	/*------------------------------------------------------------------------*/
	/** \brief  return the AIndex-indexed group of T cells if it exists,
         *  \throw  an exception otherwise.
         *
         * 	\param AIndex the index of the group of T cells
         *
         *  \return a cloud
	 */
	template<typename T> CellGroup<T> *getGroup(unsigned int AIndex);

	/*------------------------------------------------------------------------*/
	/** \brief  return the number of groups of T cells stored in the mesh.
         *
         *  \return the number of groups of T cells
	 */
	template<typename T> unsigned int getNbGroups() const;

	/*------------------------------------------------------------------------*/
	/** \brief  return an iterator on the first group of T cells of the mesh
         *
         *  \return an iterator located on the first group of T cells of the mesh
	 */
	template<typename T> group_iterator<T> groups_begin();

	/*------------------------------------------------------------------------*/
	/** \brief  return an iterator on the last group of T cells of the mesh
         *
         *  \return an iterator located on the last group of T cells of the mesh
	 */
	template<typename T> group_iterator<T> groups_end();

	/*------------------------------------------------------------------------*/
	/** \brief  Returns the ids of nodes common to two faces AF1 and AF2
         *
         * 	\param AF1  a first face
         * 	\param AF2  a second face
         *  \param AVec a vector of node ids incident to both AF1 and AF2
	 */
	std::vector<TCellID> getCommonNodes(const Face &AF1, const Face &AF2);
	static void getCommonNodes(const Face &AF1, const Face &AF2, std::vector<TCellID> &Avec);
	/*------------------------------------------------------------------------*/
	/** \brief  Returns the ids of faces shared by two nodes AN1 and AN2.
         *          Be cautious about the fact that this method relies on the
         *          assumption that the N2F connection is available and filled
         *          correctly.
         *
         * 	\param AN1 a first node
         * 	\param AN2  a second node
         *  \param AVec a vector of the face ids incident to both AN1 and AN2
	 */
	std::vector<TCellID> getCommonFaces(const Node &AN1, const Node &AN2);

	static void getCommonFaces(const Node &AN1, const Node &AN2, std::vector<TCellID> &Avec);

	/*------------------------------------------------------------------------*/
	/** \brief  Initialize the geometry classification. In other words, this
         * 			method attaches a geometry variable for each cell dimension.
         * 			After that, each cell is responsible of the value of the
         * 			geometric classification.
	 */
	void initializeGeometryClassification();

	/*------------------------------------------------------------------------*/
	/** \brief  Returns whether the classification variable for dimension ADim
         * 			exists.
         *
         *  \param ADim the dimension of cells we want to know whether the classification
         *  			exists.
	 */
	bool doesGeometricClassificationExist(int ADim);

	/*------------------------------------------------------------------------*/
	/** \brief  Provide the classification variable for dimension ADim
         *
         *  \param ADim the dimension of cells we want the classification
	 */
	Variable<cad::GeomEntity *> *getGeometricClassification(int ADim);

	/*------------------------------------------------------------------------*/
	/** \brief  Provide the geometric entity associated to ACell
         *
         *  \param ACell the mesh cell we want to get the classification
	 */
	cad::GeomEntity *getGeometricClassification(const Node &ACell);
	cad::GeomEntity *getGeometricClassification(const Edge &ACell);
	cad::GeomEntity *getGeometricClassification(const Face &ACell);
	cad::GeomEntity *getGeometricClassification(const Region &ACell);
	/*------------------------------------------------------------------------*/
	/** \brief  Update the geometric entity associated to ACell
         *
         *  \param ACell the mesh cell we want to modify the classification
	 */
	void setGeometricClassification(const Node &ACell, cad::GeomEntity *e);
	void setGeometricClassification(const Edge &ACell, cad::GeomEntity *e);
	void setGeometricClassification(const Face &ACell, cad::GeomEntity *e);
	void setGeometricClassification(const Region &ACell, cad::GeomEntity *e);

	/*------------------------------------------------------------------------*/
	/** \brief  Serializes the mesh into stream AStr
         *
         *  \param AStr an output stream where the mesh data is written
	 */
	void serialize(std::ostream &AStr);
	/*------------------------------------------------------------------------*/
	/** \brief  Unserializes a mesh from stream AStr
         *
         *  \param AStr an input stream where the mesh data is read
	 */
	void unserialize(std::istream &AStr);

#	ifdef GMDS_PARALLEL

	TInt getPartID() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Set a cell as being (or not) a master cell. If this cell is not
     * 			shared, this operation makes it shared too.
     *
     * \param ACell 	the cell to transform as being a master cell
     * \para  AMaster	transforms the cell as a master (true) or not (false)
	 */
	template<typename T> void setMaster(const TCellID &AID, const bool &AMaster = true);
	void setMaster(const Node &ACell, const bool &AMaster = true);
	void setMaster(const Edge &ACell, const bool &AMaster = true);
	void setMaster(const Face &ACell, const bool &AMaster = true);
	void setMaster(const Region &ACell, const bool &AMaster = true);

	/*------------------------------------------------------------------------*/
	/** \brief  Check if a cell is shared (master or slave)
     *
     *  \param ACell a mesh cell
     *
     *  \return true if the cell is shared
	 */
	template<typename T> bool isShared(const TCellID &AID) const;
	bool isShared(const Node &ACell) const;
	bool isShared(const Edge &ACell) const;
	bool isShared(const Face &ACell) const;
	bool isShared(const Region &ACell) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Check if a cell is a master cell
     *
     *  \param ACell a mesh cell
     *
     *  \return true if the cell is a master cell
	 */
	template<typename T> bool isMaster(const TCellID &AID) const;
	bool isMaster(const Node &ACell) const;
	bool isMaster(const Edge &ACell) const;
	bool isMaster(const Face &ACell) const;
	bool isMaster(const Region &ACell) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Check if a cell is a slave value
     *
     *  \param ACell a mesh cell
     *
     *  \return true if the cell is a slave cell
	 */
	template<typename T> bool isSlave(const TCellID &AID) const;
	bool isSlave(const Node &ACell) const;
	bool isSlave(const Edge &ACell) const;
	bool isSlave(const face &ACell) const;
	bool isSlave(const Region &ACell) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Check if a cell is a slave cell of a master on partition APID.
     *
     *  \param ACell a local mesh cell
     *  \param APID  id of a partition
     *
     *  \return true if the cell master of ACell is on partition APID
	 */
	template<typename T> bool isSlaveOf(const TCellID &AID, const TInt &APID) const;
	bool isSlaveOf(Node &ACell, const TInt &APID) const;
	bool isSlaveOf(Edge &ACell, const TInt &APID) const;
	bool isSlaveOf(Face &ACell, const TInt &APID) const;
	bool isSlaveOf(Region &ACell, const TInt &APID) const;

	/*------------------------------------------------------------------------*/
	/** \brief make a cell a slaver of the master cell defined by (AMasterID,
     * 		   AMasterPartition)
     *
     *  \param ACell 		the cell that must be shared
     *  \param AMasterID 	the local of the master in AMasterPart
     *  \param AMasterPart 	the part id where the master is located
	 */
	template<typename T> void setSlave(const TCellID &ACell, const TCellID &AMasterID, const TInt &AMasterPart);
	void setSlave(Node &ACell, const TCellID &AMasterID, const TInt &AMasterPart);
	void setSlave(Edge &ACell, const TCellID &AMasterID, const TInt &AMasterPart);
	void setSlave(Face &ACell, const TCellID &AMasterID, const TInt &AMasterPart);
	void setSlave(Region &ACell, const TCellID &AMasterID, const TInt &AMasterPart);

	/*------------------------------------------------------------------------*/
	/** \brief make a cell not a slaver knowing the part of its master cell
     *
     *  \param ACell 		the cell that must be shared
     *  \param AMasterPart 	the part id where the master is located
	 */
	template<typename T> void unsetSlave(const TCellID &ACell, const TInt &AMasterPart);
	void unsetSlave(Node &ACell, const TInt &AMasterPart);
	void unsetSlave(Edge &ACell, const TInt &AMasterPart);
	void unsetSlave(Face &ACell, const TInt &AMasterPart);
	void unsetSlave(Region &ACell, const TInt &AMasterPart);

	/*------------------------------------------------------------------------*/
	/** \brief make a cell not a slaver without knowing the part of its
     * 		   master cell
     *
     *  \param ACell 		the cell that must be shared
	 */
	template<typename T> void unsetSlave(const TCellID &ACell);
	void unsetSlave(Node &ACell);
	void unsetSlave(Edge &ACell);
	void unsetSlave(Face &ACell);
	void unsetSlave(Region &ACell);

	/*------------------------------------------------------------------------*/
	/** \brief Get master data. If ACell is not a slave cell, this operation
     * 		   returns false
     *
     *  \param ACell  the cell that we want to collect the data
     *  \param AMPart the part the  master belongs to
     * 	\param AMID   the id of the master in AMPart
     *
     *  \return true if ACell is a slave cell, false otherwise
	 */
	template<typename T> bool getMasterData(const TCellID &AID, TInt &AMPart, TCellID &AMID);
	bool getMasterData(const Node &ACell, TInt &AMPart, TCellID &AMID);
	bool getMasterData(const Edge &ACell, TInt &AMPart, TCellID &AMID);
	bool getMasterData(const Face &ACell, TInt &AMPart, TCellID &AMID);
	bool getMasterData(const Region &ACell, TInt &AMPart, TCellID &AMID);

	/*------------------------------------------------------------------------*/
	/** \brief Get interface with partition APartID.
     *
     *  \param  ACells  the pair of ids we want to cell with
     *  							 (local slave id, distant master id)
     *  \param  APartID the partition we want the shared cells with
     *  \param  ADim	the dimension of the cells we work on
     *
     *  \return the collection of shared cells with this partition
	 */
	void getInterface(SmartVector<std::pair<TCellID, TCellID>> &ACells, const TInt &APartID, const int &ADim) const;

	/*------------------------------------------------------------------------*/
	/** \brief Clear slave interfaces.
     * */
	void clearInterface(const TInt &APartitionID, const TInt &ADimCell);
	void clearInterface(const TInt &ADimCell);

	void printShareInfo(std::ostream &str);

	/*------------------------------------------------------------------------*/
	/** \brief Add the part APart as being a part containing slave ADim-cell of
     * 		   whose master is on this mesh.
     * \param ADim  the dimension of slave cell
     * \param APart the id of the part containing the slave cell
     * */
	void addSlavePart(const int ADim, const int APart);

	/*------------------------------------------------------------------------*/
	/** \brief Get the parts containing slaves
     * \param ADim  the dimension of slave cell
	 */
	std::set<TInt> getSlaveParts(const int ADim) const;
	std::set<TInt> getSlaveParts() const;

#	endif     // GMDS_PARALLEL

 protected:
	/*------------------------------------------------------------------------*/
	/** \brief  Resize node container when the user exactly know the ids of
         * 		    cells he wants to create. Warning, in this case, the user must
         * 			take care of creating a coherent data structure. Otherwise,
         * 			segmentation faults might occur (not meet yet but it could
         * 			happen).
         *
         * 			Just before resizing all the ids are removed and all the cells
         * 			are erased in their respective memory allocators.
         *
         *  \param  AMaxID the max id stored for ADim-cells
	 */
	void clearAndResizeNodeIDContainer(TInt AMaxID);

	/*------------------------------------------------------------------------*/
	/** \brief  Resize edge container when the user exactly know the ids of
         * 		    cells he wants to create. Warning, in this case, the user must
         * 			take care of creating a coherent data structure. Otherwise,
         * 			segmentation faults might occur (not meet yet but it could
         * 			happen).
         *
         * 			Just before resizing all the ids are removed and all the cells
         * 			are erased in their respective memory allocators.
         *
         *  \param  AMaxID the max id stored for ADim-cells
	 */
	void clearAndResizeEdgeIDContainer(TInt AMaxID);

	/*------------------------------------------------------------------------*/
	/** \brief  Resize face container when the user exactly know the ids of
         * 		    cells he wants to create. Warning, in this case, the user must
         * 			take care of creating a coherent data structure. Otherwise,
         * 			segmentation faults might occur (not meet yet but it could
         * 			happen).
         *
         * 			Just before resizing all the ids are removed and all the cells
         * 			are erased in their respective memory allocators.
         *
         *  \param  AMaxID the max id stored for ADim-cells
	 */
	void clearAndResizeFaceIDContainer(TInt AMaxID);

	/*------------------------------------------------------------------------*/
	/** \brief  Resize region container when the user exactly know the ids of
         * 		    cells he wants to create. Warning, in this case, the user must
         * 			take care of creating a coherent data structure. Otherwise,
         * 			segmentation faults might occur (not meet yet but it could
         * 			happen).
         *
         * 			Just before resizing all the ids are removed and all the cells
         * 			are erased in their respective memory allocators.
         *
         *  \param  AMaxID the max id stored for ADim-cells
	 */
	void clearAndResizeRegionIDContainer(TInt AMaxID);

	/*------------------------------------------------------------------------*/
	/** \brief  Update the mesh id containers. This operations is NECESSARY to
         * 			keep valid meshes when we let the user to insert cells with
         * 			specified ids. Each time a new cell with its specified id is
         * 			added, the corresponding id container is not updated to keep
         * 			a valid index of its free spaces. As a consequence, this
         * 			operation is necessary as being a final step to get a valid
         * 			container.
	 */
	void updateIDContainers();

	/*------------------------------------------------------------------------*/
	/** \brief  Add a Node into the mesh. In 2D, only AX and AY are used.
         *
         * 			This operation must be called by reader class only. It is
         * 			particularly true in a distributed memory context. Otherwise,
         * 		    you encounter some troubles with the global ids.
         *
         * 			The size of the id container must be specified before adding
         * 		    nodes with ids. At the end, an update of the id container must
         * 			be done to have a coherent container.
         *
         *  \param AX X coordinate
         *  \param AY Y coordinate
         *  \param AZ Z coordinate
         *  \param AGID global id we want to assign to give to the new cell
         *
         *  \return a pointer on the new Node
	 */
	Node newNodeWithID(const TCoord &AX, const TCoord &AY, const TCoord &AZ, const TCellID &AGID);

	/*------------------------------------------------------------------------*/
	/** \brief  Add an edge defined by two vertice ids
         * 			This operation must be called by reader class only. It is
         * 			particularly true in a distributed memory context. Otherwise,
         * 		    you encounter some troubles with the global ids.
         *
         * 			The size of the id container must be specified before adding
         * 		    nodes with ids. At the end, an update of the id container must
         * 			be done to have a coherent container.
         *
         *  \param AV1 first  Node id
         *  \param AV2 second Node id
         *  \param AGID global id we want to assign to give to the new cell
         *
         *  \return a pointer on the new edge
	 */
	Edge newEdgeWithID(const TCellID &AV1, const TCellID &AV2, const TCellID &AGID);

	/*------------------------------------------------------------------------*/
	/** \brief  Add a triangle defined by three vertex ids.
         * 			This operation must be called by reader class only. It is
         * 			particularly true in a distributed memory context. Otherwise,
         * 		    you encounter some troubles with the global ids.
         *
         * 			The size of the id container must be specified before adding
         * 		    nodes with ids. At the end, an update of the id container must
         * 			be done to have a coherent container.
         *
         *  \param AN1 Node 1 id
         *  \param AN2 Node 2 id
         *  \param AN3 Node 3 id
         *  \param AGID global id we want to assign to give to the new cell
         *
         *  \return a pointer on the new face
	 */
	Face newTriangleWithID(const TCellID &AN1, const TCellID &AN2, const TCellID &AN3, const TCellID &AGID);

	/*------------------------------------------------------------------------*/
	/** \brief  Add a quad defined by four ordered vertices ids.
         * 			This operation must be called by reader class only. It is
         * 			particularly true in a distributed memory context. Otherwise,
         * 		    you encounter some troubles with the global ids.
         *
         * 			The size of the id container must be specified before adding
         * 		    nodes with ids. At the end, an update of the id container must
         * 			be done to have a coherent container.
         *
         *  \param AN1 Node 1 id
         *  \param AN2 Node 2 id
         *  \param AN3 Node 3 id
         *  \param AN4 Node 4 id
         *  \param AGID global id we want to assign to give to the new cell
         *
         *  \return a pointer on the new face
	 */
	Face newQuadWithID(const TCellID &AN1, const TCellID &AN2, const TCellID &AN3, const TCellID &AN4, const TCellID &AGID);

	/*------------------------------------------------------------------------*/
	/** \brief  Add a polygon defined by an ordered collection of vertex ids.
         * 			This operation must be called by reader class only. It is
         * 			particularly true in a distributed memory context. Otherwise,
         * 		    you encounter some troubles with the global ids.
         *
         * 			The size of the id container must be specified before adding
         * 		    nodes with ids. At the end, an update of the id container must
         * 			be done to have a coherent container.
         *
         *  \param ANodes a collection of vertex ids
         *  \param AGID global id we want to assign to give to the new cell
         *
         *  \return a pointer on the new face
	 */
	Face newPolygonWithID(const std::vector<TCellID> &ANodes, const TCellID &AGID);

	/*------------------------------------------------------------------------*/
	/** \brief  Add a face defined by an ordered collection of vertex ids.
         * 			Contrary to the polygon builder, the most appropriate face is
         * 			built in this case.
         *
         * 			This operation must be called by reader class only. It is
         * 			particularly true in a distributed memory context. Otherwise,
         * 		    you encounter some troubles with the global ids.
         *
         * 			The size of the id container must be specified before adding
         * 		    nodes with ids. At the end, an update of the id container must
         * 			be done to have a coherent container.
         *
         *  \param ANodes a collection of vertex ids
         *  \param AGID global id we want to assign to give to the new cell
         *
         *  \return a pointer on the new face
	 */
	Face newFaceWithID(const std::vector<TCellID> &AIDs, const TCellID &AGID);

	/*------------------------------------------------------------------------*/
	/** \brief  Add a tetrahedron defined by four vertices ids
         * 			This operation must be called by reader class only. It is
         * 			particularly true in a distributed memory context. Otherwise,
         * 		    you encounter some troubles with the global ids.
         *
         * 			The size of the id container must be specified before adding
         * 		    nodes with ids. At the end, an update of the id container must
         * 			be done to have a coherent container.
         *
         *  \param AN1 Node 1 id
         *  \param AN2 Node 2 id
         *  \param AN3 Node 3 id
         *  \param AN4 Node 4 id
         *  \param AGID global id we want to assign to give to the new cell
         *
         *  \return a pointer on the new region
	 */
	Region newTetWithID(const TCellID &AN1, const TCellID &AN2, const TCellID &AN3, const TCellID &AN4, const TCellID &AGID);

	/*------------------------------------------------------------------------*/
	/** \brief  Add a pyramid defined by 5 vertices id whose the fourth first
         * 			vertices define the square face.
         *
         *			              5
         *
         * 				  	2 ----------- 3
         * 				   /		     /
         *			      /             /
         *			     1 ----------- 4
         *
         * 			This operation must be called by reader class only. It is
         * 			particularly true in a distributed memory context. Otherwise,
         * 		    you encounter some troubles with the global ids.
         *
         * 			The size of the id container must be specified before adding
         * 		    nodes with ids. At the end, an update of the id container must
         * 			be done to have a coherent container.
         *
         *  \param AN1 Node 1 id
         *  \param AN2 Node 2 id
         *  \param AN3 Node 3 id
         *  \param AN4 Node 4 id
         *  \param AN5 Node 5 id
         *  \param AGID global id we want to assign to give to the new cell
         *
         *  \return a pointer on the new region
	 */
	Region newPyramidWithID(const TCellID &AN1, const TCellID &AN2, const TCellID &AN3, const TCellID &AN4, const TCellID &AN5, const TCellID &AGID);

	/*------------------------------------------------------------------------*/
	/** \brief  Add a prism defined by 6 vertices.
         *
         * 			This operation must be called by reader class only. It is
         * 			particularly true in a distributed memory context. Otherwise,
         * 		    you encounter some troubles with the global ids.
         *
         * 			The size of the id container must be specified before adding
         * 		    nodes with ids. At the end, an update of the id container must
         * 			be done to have a coherent container.
         *
         *  \param AN1 Node 1 id
         *  \param AN2 Node 2 id
         *  \param AN3 Node 3 id
         *  \param AN4 Node 4 id
         *  \param AN5 Node 5 id
         *  \param AN6 Node 6 id
         *  \param AGID global id we want to assign to give to the new cell
         *
         *  \return a pointer on the new region
	 */
	Region newPrism3WithID(
	   const TCellID &AN1, const TCellID &AN2, const TCellID &AN3, const TCellID &AN4, const TCellID &AN5, const TCellID &AN6, const TCellID &AGID);

	/*------------------------------------------------------------------------*/
	/** \brief  Add a hexahedron defined by eight vertices ids whose order is
         *
         * 					2 ----------- 3
         * 				   /|            /|
         *			      / |           / |
         *			     1 ----------- 4  |
         *			     |  |          |  |
         * 				 |	6 ---------|- 7
         * 				 | /		   | /
         *			     |/            |/
         *			     5 ----------- 8
         *
         * 			This operation must be called by reader class only. It is
         * 			particularly true in a distributed memory context. Otherwise,
         * 		    you encounter some troubles with the global ids.
         *
         * 			The size of the id container must be specified before adding
         * 		    nodes with ids. At the end, an update of the id container must
         * 			be done to have a coherent container.
         *
         *  \param AN1 Node 1 id
         *  \param AN2 Node 2 id
         *  \param AN3 Node 3 id
         *  \param AN4 Node 4 id
         *  \param AN5 Node 5 id
         *  \param AN6 Node 6 id
         *  \param AN7 Node 7 id
         *  \param AN8 Node 8 id
         *  \param AGID global id we want to assign to give to the new cell
         *
         *  \return a pointer on the new region
	 */
	Region newHexWithID(const TCellID &AN1,
	                    const TCellID &AN2,
	                    const TCellID &AN3,
	                    const TCellID &AN4,
	                    const TCellID &AN5,
	                    const TCellID &AN6,
	                    const TCellID &AN7,
	                    const TCellID &AN8,
	                    const TCellID &AGID);

	/*------------------------------------------------------------------------*/
	/** \brief Change the mesh model and build the new entities/adjacnecies
         *
         * \param the mesh model
	 */
	void changeModelWithDoctor(const MeshModel &AModel);

	/*------------------------------------------------------------------------*/
	/** \brief Change the mesh model but DO NOT build the new entities/adjacnecies
         *
         * \param the mesh model
	 */
	void changeModelWithoutDoctor(const MeshModel &AModel);

 protected:
	/** implemented mesh model */
	MeshModel m_model;

	/** Cells container */
	NodeContainer *m_nodes_container;
	EdgeContainer *m_edges_container;
	FaceContainer *m_faces_container;
	RegionContainer *m_regions_container;

	VariableManager *m_node_variable_manager;
	VariableManager *m_edge_variable_manager;
	VariableManager *m_face_variable_manager;
	VariableManager *m_region_variable_manager;

	std::list<CellGroup<Node> *> m_clouds;
	std::list<CellGroup<Edge> *> m_lines;
	std::list<CellGroup<Face> *> m_surfaces;
	std::list<CellGroup<Region> *> m_volumes;

	/** geometrical classification. Useful in some cases,
         *  WARNING it must be activated in the mesh interface.
	 */
	Variable<cad::GeomEntity *> *classification[4];

	/** Boolean marks management
	 */
	Variable<Marks32> *m_marks[4];

	/* marks currently used*/
	Marks32 m_usedMarks_nodes;
	Marks32 m_usedMarks_edges;
	Marks32 m_usedMarks_faces;
	Marks32 m_usedMarks_regions;
	/* mask of the marks (altered by negateMaskMark()) */
	Marks32 m_maskMarks_nodes;
	Marks32 m_maskMarks_edges;
	Marks32 m_maskMarks_faces;
	Marks32 m_maskMarks_regions;

	/* indicates free marks*/
	TInt m_marks_nodes[32] {};
	TInt m_marks_edges[32] {};
	TInt m_marks_faces[32] {};
	TInt m_marks_regions[32] {};
	/* number of used marks*/
	TInt m_nbUsedMarks_nodes;
	TInt m_nbUsedMarks_edges;
	TInt m_nbUsedMarks_faces;
	TInt m_nbUsedMarks_regions;

#	ifdef _DEBUG
	TInt m_maxNbUsedMarks_nodes;
	TInt m_maxNbUsedMarks_edges;
	TInt m_maxNbUsedMarks_faces;
	TInt m_maxNbUsedMarks_regions;
#	endif     // _DEBUG
};

template<typename T, ECellType E> Variable<T> *
Mesh::newVariable(const std::string &AName)
{
	Variable<T> *v;
	std::vector<TInt> ids;
	switch (E) {
	case GMDS_NODE: {
		for (auto i : nodes()) {
			ids.push_back(i);
		}
		v = m_node_variable_manager->newVariable<T>(AName, m_nodes_container->capacity(), &ids);
	} break;
	case GMDS_EDGE: {
		for (auto i : edges()) {
			ids.push_back(i);
		}
		v = m_edge_variable_manager->newVariable<T>(AName, m_edges_container->capacity(), &ids);
	} break;
	case GMDS_FACE: {
		for (auto i : faces()) {
			ids.push_back(i);
		}
		v = m_face_variable_manager->newVariable<T>(AName, m_faces_container->capacity(), &ids);
	} break;
	case GMDS_REGION: {
		for (auto i : regions()) {
			ids.push_back(i);
		}
		v = m_region_variable_manager->newVariable<T>(AName, m_regions_container->capacity(), &ids);
	} break;
	default: throw GMDSException("Unmanaged type of value -> impossible to create a variable");
	}

	return v;
}
/*----------------------------------------------------------------------------*/
template<typename T, ECellType E> Variable<T> *
Mesh::getVariable(const std::string &AName) const
{
	Variable<T> *v;
	std::vector<int> ids;
	switch (E) {
	case GMDS_NODE: {
		v = m_node_variable_manager->getVariable<T>(AName);
	} break;
	case GMDS_EDGE: {
		v = m_edge_variable_manager->getVariable<T>(AName);
	} break;
	case GMDS_FACE: {
		v = m_face_variable_manager->getVariable<T>(AName);
	} break;
	case GMDS_REGION: {
		v = m_region_variable_manager->getVariable<T>(AName);
	} break;
	default: throw GMDSException("Unmanaged type of value -> impossible to access to a variable");
	}
	return v;
}
/*----------------------------------------------------------------------------*/
template<typename T, ECellType E> Variable<T> *
Mesh::getOrCreateVariable(const std::string &AName)
{
	Variable<T> *v = NULL;
	try {
		v = getVariable<T, E>(AName);
	}
	catch (GMDSException &e) {
		v = newVariable<T, E>(AName);
	}
	return v;
}
/*----------------------------------------------------------------------------*/
template<> GMDSIg_API CellGroup<Node> *Mesh::newGroup<Node>(const std::string &AName);
template<> GMDSIg_API CellGroup<Edge> *Mesh::newGroup<Edge>(const std::string &AName);
template<> GMDSIg_API CellGroup<Face> *Mesh::newGroup<Face>(const std::string &AName);
template<> GMDSIg_API CellGroup<Region> *Mesh::newGroup<Region>(const std::string &AName);
/*----------------------------------------------------------------------------*/
template<> GMDSIg_API void Mesh::deleteGroup<Node>(CellGroup<Node> *ACloud);
template<> GMDSIg_API void Mesh::deleteGroup<Edge>(CellGroup<Edge> *ACloud);
template<> GMDSIg_API void Mesh::deleteGroup<Face>(CellGroup<Face> *ACloud);
template<> GMDSIg_API void Mesh::deleteGroup<Region>(CellGroup<Region> *ACloud);
/*----------------------------------------------------------------------------*/
template<> GMDSIg_API CellGroup<Node> *Mesh::getGroup<Node>(const std::string &AName);
template<> GMDSIg_API CellGroup<Edge> *Mesh::getGroup<Edge>(const std::string &AName);
template<> GMDSIg_API CellGroup<Face> *Mesh::getGroup<Face>(const std::string &AName);
template<> GMDSIg_API CellGroup<Region> *Mesh::getGroup<Region>(const std::string &AName);
/*----------------------------------------------------------------------------*/
template<> GMDSIg_API CellGroup<Node> *Mesh::getGroup<Node>(unsigned int AIndex);
template<> GMDSIg_API CellGroup<Edge> *Mesh::getGroup<Edge>(unsigned int AIndex);
template<> GMDSIg_API CellGroup<Face> *Mesh::getGroup<Face>(unsigned int AIndex);
template<> GMDSIg_API CellGroup<Region> *Mesh::getGroup<Region>(unsigned int AIndex);
/*----------------------------------------------------------------------------*/
template<> GMDSIg_API unsigned int Mesh::getNbGroups<Node>() const;
template<> GMDSIg_API unsigned int Mesh::getNbGroups<Edge>() const;
template<> GMDSIg_API unsigned int Mesh::getNbGroups<Face>() const;
template<> GMDSIg_API unsigned int Mesh::getNbGroups<Region>() const;
/*----------------------------------------------------------------------------*/
template<> GMDSIg_API Mesh::group_iterator<Node> Mesh::groups_begin<Node>();
template<> GMDSIg_API Mesh::group_iterator<Edge> Mesh::groups_begin<Edge>();
template<> GMDSIg_API Mesh::group_iterator<Face> Mesh::groups_begin<Face>();
template<> GMDSIg_API Mesh::group_iterator<Region> Mesh::groups_begin<Region>();
/*----------------------------------------------------------------------------*/
template<> GMDSIg_API Mesh::group_iterator<Node> Mesh::groups_end<Node>();
template<> GMDSIg_API Mesh::group_iterator<Edge> Mesh::groups_end<Edge>();
template<> GMDSIg_API Mesh::group_iterator<Face> Mesh::groups_end<Face>();
template<> GMDSIg_API Mesh::group_iterator<Region> Mesh::groups_end<Region>();

/*----------------------------------------------------------------------------*/
template<> GMDSIg_API void Mesh::getAll<Node>(std::vector<Node> &AVec) const;

template<> GMDSIg_API void Mesh::getAll<Edge>(std::vector<Edge> &AVec) const;

template<> GMDSIg_API void Mesh::getAll<Face>(std::vector<Face> &AVec) const;

template<> GMDSIg_API void Mesh::getAll<Region>(std::vector<Region> &AVec) const;

/*----------------------------------------------------------------------------*/
template<> GMDSIg_API bool Mesh::has<Node>(const TCellID &AID) const;
template<> GMDSIg_API bool Mesh::has<Edge>(const TCellID &AID) const;
template<> GMDSIg_API bool Mesh::has<Face>(const TCellID &AID) const;
template<> GMDSIg_API bool Mesh::has<Region>(const TCellID &AID) const;
/*----------------------------------------------------------------------------*/
template<> GMDSIg_API Node Mesh::get<Node>(const TCellID &AID) const;
template<> GMDSIg_API Edge Mesh::get<Edge>(const TCellID &AID) const;
template<> GMDSIg_API Face Mesh::get<Face>(const TCellID &AID) const;
template<> GMDSIg_API Region Mesh::get<Region>(const TCellID &AID) const;

/*----------------------------------------------------------------------------*/
template<> GMDSIg_API void Mesh::freeMark<Node>(TInt AMarkNumber);
template<> GMDSIg_API void Mesh::freeMark<Edge>(TInt AMarkNumber);
template<> GMDSIg_API void Mesh::freeMark<Face>(TInt AMarkNumber);
template<> GMDSIg_API void Mesh::freeMark<Region>(TInt AMarkNumber);

/*----------------------------------------------------------------------------*/
template<> GMDSIg_API void Mesh::unmarkAll<Node>(TInt AMarkNumber);
template<> GMDSIg_API void Mesh::unmarkAll<Edge>(TInt AMarkNumber);
template<> GMDSIg_API void Mesh::unmarkAll<Face>(TInt AMarkNumber);
template<> GMDSIg_API void Mesh::unmarkAll<Region>(TInt AMarkNumber);

/*----------------------------------------------------------------------------*/
template<> GMDSIg_API void Mesh::unmark<Node>(const TCellID &ACellID, TInt AMarkNumber);
template<> GMDSIg_API void Mesh::unmark<Edge>(const TCellID &ACellID, TInt AMarkNumber);
template<> GMDSIg_API void Mesh::unmark<Face>(const TCellID &ACellID, TInt AMarkNumber);
template<> GMDSIg_API void Mesh::unmark<Region>(const TCellID &ACellID, TInt AMarkNumber);

/*----------------------------------------------------------------------------*/
template<> GMDSIg_API void Mesh::mark<Node>(const TCellID &ACellID, TInt AMarkNumber);
template<> GMDSIg_API void Mesh::mark<Edge>(const TCellID &ACellID, TInt AMarkNumber);
template<> GMDSIg_API void Mesh::mark<Face>(const TCellID &ACellID, TInt AMarkNumber);
template<> GMDSIg_API void Mesh::mark<Region>(const TCellID &ACellID, TInt AMarkNumber);
/*----------------------------------------------------------------------------*/
template<> GMDSIg_API void Mesh::markTo<Node>(const TCellID &ACellID, TInt AMarkNumber, bool AState);
template<> GMDSIg_API void Mesh::markTo<Edge>(const TCellID &ACellID, TInt AMarkNumber, bool AState);
template<> GMDSIg_API void Mesh::markTo<Face>(const TCellID &ACellID, TInt AMarkNumber, bool AState);
template<> GMDSIg_API void Mesh::markTo<Region>(const TCellID &ACellID, TInt AMarkNumber, bool AState);

/*----------------------------------------------------------------------------*/
template<> GMDSIg_API bool Mesh::isMarked<Node>(const TCellID &ACellID, TInt AMarkNumber) const;
template<> GMDSIg_API bool Mesh::isMarked<Edge>(const TCellID &ACellID, TInt AMarkNumber) const;
template<> GMDSIg_API bool Mesh::isMarked<Face>(const TCellID &ACellID, TInt AMarkNumber) const;
template<> GMDSIg_API bool Mesh::isMarked<Region>(const TCellID &ACellID, TInt AMarkNumber) const;
/*----------------------------------------------------------------------------*/
template<> GMDSIg_API void Mesh::negateMaskMark<Node>(TInt AMarkNumber);
template<> GMDSIg_API void Mesh::negateMaskMark<Edge>(TInt AMarkNumber);
template<> GMDSIg_API void Mesh::negateMaskMark<Face>(TInt AMarkNumber);
template<> GMDSIg_API void Mesh::negateMaskMark<Region>(TInt AMarkNumber);
/*----------------------------------------------------------------------------*/
template<> GMDSIg_API TInt Mesh::newMark<Node>();
template<> GMDSIg_API TInt Mesh::newMark<Edge>();
template<> GMDSIg_API TInt Mesh::newMark<Face>();
template<> GMDSIg_API TInt Mesh::newMark<Region>();
/*----------------------------------------------------------------------------*/
template<> GMDSIg_API Marks32 Mesh::getMarks<Node>(const Node &ACell) const;
template<> GMDSIg_API Marks32 Mesh::getMarks<Edge>(const Edge &ACell) const;
template<> GMDSIg_API Marks32 Mesh::getMarks<Face>(const Face &ACell) const;
template<> GMDSIg_API Marks32 Mesh::getMarks<Region>(const Region &ACell) const;
/*----------------------------------------------------------------------------*/
}     // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_IGMESH_H_ */
/*----------------------------------------------------------------------------*/
