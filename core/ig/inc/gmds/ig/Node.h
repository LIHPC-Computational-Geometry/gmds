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
	~Node() override;

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
	int dim() const override {return 0;}

	/*------------------------------------------------------------------------*/
    /** \brief Accessor th the number of incident nodes, edges, faces and
     * 		   adjacent regions
     */
	TInt nbNodes()   const override;
	TInt nbEdges()   const override;
	TInt nbFaces()   const override;
	TInt nbRegions() const override;

	/*------------------------------------------------------------------------*/
	/** \brief  Compute the center of the region
 	 *
 	 * \return the center of the region
         */
         math::Point center() const override;
	/*------------------------------------------------------------------------*/
    /** \brief  Accessor to the incident cells. Only the non-null cells are
     * 			provided. T can be Node, Edge, Face or Region.
     */
	void delegateGet(std::vector<Node>&   ACells) const override;
	void delegateGet(std::vector<Edge>&   ACells) const override;
	void delegateGet(std::vector<Face>&   ACells) const override;
	void delegateGet(std::vector<Region>& ACells) const override;

	void delegateGetNodeIDs  (std::vector<TCellID>& ACells) const override;
	void delegateGetEdgeIDs  (std::vector<TCellID>& ACells) const override;
	void delegateGetFaceIDs  (std::vector<TCellID>& ACells) const override;
	void delegateGetRegionIDs(std::vector<TCellID>& ACells) const override;


	void delegateGetAll(std::vector<Node>&   ACells) const override;
	void delegateGetAll(std::vector<Edge>&   ACells) const override;
	void delegateGetAll(std::vector<Face>&   ACells) const override;
	void delegateGetAll(std::vector<Region>& ACells) const override;

	void delegateGetAllNodeIDs  (std::vector<TCellID>& ACells) const override;
        void delegateGetAllEdgeIDs  (std::vector<TCellID>& ACells) const override;
        void delegateGetAllFaceIDs  (std::vector<TCellID>& ACells) const override;
        void delegateGetAllRegionIDs(std::vector<TCellID>& ACells) const override;

	void delegateSetNodeIDs(const std::vector<TCellID>& ACells) override ;
	void delegateSetEdgeIDs(const std::vector<TCellID>& ACells) override ;
	void delegateSetFaceIDs(const std::vector<TCellID>& ACells) override ;
	void delegateSetRegionIDs(const std::vector<TCellID>& ACells) override ;

	void delegateNodeAdd(TCellID AElt) override;
	void delegateEdgeAdd(TCellID AElt) override;
	void delegateFaceAdd(TCellID AElt) override;
	void delegateRegionAdd(TCellID AElt) override;

	void delegateNodeRemove(TCellID AElt) override;
	void delegateEdgeRemove(TCellID AElt) override;
	void delegateFaceRemove(TCellID AElt) override;
	void delegateRegionRemove(TCellID AElt) override;

	void delegateNodeReplace  (TCellID AID1, TCellID AID2) override;
	void delegateEdgeReplace  (TCellID AID1, TCellID AID2) override;
	void delegateFaceReplace  (TCellID AID1, TCellID AID2) override;
	void delegateRegionReplace(TCellID AID1, TCellID AID2) override;


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
	void setX(TCoord AVal);
	void setY(TCoord AVal);
	void setZ(TCoord AVal);
	void setXYZ(TCoord AX, TCoord AY, TCoord AZ);

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
