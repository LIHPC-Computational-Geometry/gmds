/*----------------------------------------------------------------------------*/
/*
 * Node.h
 *
 *  Created on: 5 févr. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_NODE_H_
#define GMDS_NODE_H_
/*----------------------------------------------------------------------------*/
#include <gmds/math/Point.h>
#include "Cell.h"
#include "GMDSIg_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
class Mesh;
class NodeContainer;


/*----------------------------------------------------------------------------*/
/** \class Node
 *
 *  \brief A node instance is an object that provided an object-type access to
 *  	   the data representing a mesh node.
 *
 */
class GMDSIg_API Node : public Cell{
public:
	/*------------------------------------------------------------------------*/
	/** \brief Default constructor. Used for stl container initialization
	 */
	Node();

	/*------------------------------------------------------------------------*/
	/** \brief Constructor
	 *
	 * \param AMesh the mesh containing this cell
	 * \param AID the cell id
	 * \param AX X coordinate
	 * \param AY Y coordinate
	 * \param AZ Z coordinate
	 */
	Node(Mesh* AMesh, const TCellID& AID,
			const TCoord& AX, const TCoord& AY, const TCoord& AZ);

	/*------------------------------------------------------------------------*/
	/** \brief Constructor
	 *
	 * \param AMesh the mesh containing this cell
	 * \param AID the cell id
	 * \param APt a point
	 */
	Node(Mesh* AMesh,const TCellID& AID, const math::Point& APt);

	/*------------------------------------------------------------------------*/
	/** \brief Copy constructor
	 */
	Node(const Node&);

	/*------------------------------------------------------------------------*/
	/** \brief Destructor
	 */
	virtual ~Node();

	void operator=(const Node& ANode);

        /*------------------------------------------------------------------------*/
        /** \brief Overide operator==. It is id-based.
         *
         * \param ANode a node
         */
	bool operator==(const Node& ANode) const;

        /*------------------------------------------------------------------------*/
        /** \brief Overide operator!=. It is id-based.
         *
         * \param ANode a node
         */
        bool operator!=(const Node& ANode) const;

        /*------------------------------------------------------------------------*/
        /** \brief Overide operator<. It is id-based.
 	 *
 	 * \param ANode a node	 
         */
	bool operator<(const Node& ANode) const {return this->m_id < ANode.m_id;}

	/*------------------------------------------------------------------------*/
    /** \brief  Accesor to the cell dim.
     */
	virtual int dim() const {return 0;}

	/*------------------------------------------------------------------------*/
    /** \brief Accessor th the number of incident nodes, edges, faces and
     * 		   adjacent regions
     */
	virtual TInt nbNodes()   const;
	virtual TInt nbEdges()   const;
	virtual TInt nbFaces()   const;
	virtual TInt nbRegions() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Compute the center of the region
 	 *
 	 * \return the center of the region
         */
         math::Point center() const;
	/*------------------------------------------------------------------------*/
    /** \brief  Accessor to the incident cells. Only the non-null cells are
     * 			provided. T can be Node, Edge, Face or Region.
     */
	virtual void delegateGet(std::vector<Node>&   ACells) const;
	virtual void delegateGet(std::vector<Edge>&   ACells) const;
	virtual void delegateGet(std::vector<Face>&   ACells) const;
	virtual void delegateGet(std::vector<Region>& ACells) const;

	virtual void delegateGetNodeIDs  (std::vector<TCellID>& ACells) const;
	virtual void delegateGetEdgeIDs  (std::vector<TCellID>& ACells) const;
	virtual void delegateGetFaceIDs  (std::vector<TCellID>& ACells) const;
	virtual void delegateGetRegionIDs(std::vector<TCellID>& ACells) const;


	virtual void delegateGetAll(std::vector<Node>&   ACells) const;
	virtual void delegateGetAll(std::vector<Edge>&   ACells) const;
	virtual void delegateGetAll(std::vector<Face>&   ACells) const;
	virtual void delegateGetAll(std::vector<Region>& ACells) const;

	virtual void delegateGetAllNodeIDs  (std::vector<TCellID>& ACells) const;
        virtual void delegateGetAllEdgeIDs  (std::vector<TCellID>& ACells) const;
        virtual void delegateGetAllFaceIDs  (std::vector<TCellID>& ACells) const;
        virtual void delegateGetAllRegionIDs(std::vector<TCellID>& ACells) const;

	virtual void delegateSetNodeIDs(const std::vector<TCellID>& ACells) ;
	virtual void delegateSetEdgeIDs(const std::vector<TCellID>& ACells) ;
	virtual void delegateSetFaceIDs(const std::vector<TCellID>& ACells) ;
	virtual void delegateSetRegionIDs(const std::vector<TCellID>& ACells) ;

	virtual void delegateNodeAdd(TCellID AElt);
	virtual void delegateEdgeAdd(TCellID AElt);
	virtual void delegateFaceAdd(TCellID AElt);
	virtual void delegateRegionAdd(TCellID AElt);

	virtual void delegateNodeRemove(TCellID AElt);
	virtual void delegateEdgeRemove(TCellID AElt);
	virtual void delegateFaceRemove(TCellID AElt);
	virtual void delegateRegionRemove(TCellID AElt);

	virtual void delegateNodeReplace  (TCellID AID1, TCellID AID2);
	virtual void delegateEdgeReplace  (TCellID AID1, TCellID AID2);
	virtual void delegateFaceReplace  (TCellID AID1, TCellID AID2);
	virtual void delegateRegionReplace(TCellID AID1, TCellID AID2);


	/*------------------------------------------------------------------------*/
    /** \brief  Accesor to the cell point location.
     */
	math::Point point() const;

	/*------------------------------------------------------------------------*/
    /** \brief  Accesor to the cell point location.
     */
	void setPoint(const math::Point& APnt);
	/*------------------------------------------------------------------------*/
    /** \brief  Accesor to the node coordinates
     */
	TCoord X() const;
	TCoord Y() const;
	TCoord Z() const;
	TCoord& X();
	TCoord& Y();
	TCoord& Z();
	void setX(const TCoord AVal);
	void setY(const TCoord AVal);
	void setZ(const TCoord AVal);
	void setXYZ(const TCoord AX, const TCoord AY, const TCoord AZ);

	GMDSIg_API friend std::ostream& operator<<(std::ostream& AStream, const Node& AN);

protected:

	/** A link to the generic node container of the owner*/
	NodeContainer* m_nodes_container;
};
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_NODE_H_ */
/*----------------------------------------------------------------------------*/
