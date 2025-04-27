/*----------------------------------------------------------------------------*/
/*
 * Edge.h
 *
 *  Created on: 19 May 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_EDGE_H_
#define GMDS_EDGE_H_
/*----------------------------------------------------------------------------*/
#include "Cell.h"
#include "GMDSIg_export.h"
#include <gmds/math/Point.h>
#include <gmds/math/Segment.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
class Mesh;
class EdgeContainer;
class FaceContainer;
/*----------------------------------------------------------------------------*/
/** \class Edge
 *
 *  \brief An edge instance is an object that provided an object-type access to
 *  	   the data representing a mesh edge.
 *
 */
class GMDSIg_API Edge : public Cell{
public:

	/*------------------------------------------------------------------------*/
	/** \brief Default constructor. Used for stl container initialization
	 */
	Edge();

	/*------------------------------------------------------------------------*/
	/** \brief Constructor
	 *
	 * \param AMesh the mesh containing this cell
	 * \param AType the type of cell
	 * \param AID the cell id
	 */
	Edge(Mesh* AMesh, const TCellID& AID);

	/*------------------------------------------------------------------------*/
	/** \brief Copy Constructor
	 *
	 * \param AE the edge to be copied
	 */
	Edge(const Edge& AE);

	/*------------------------------------------------------------------------*/
	/** \brief Destructor
	 */
	~Edge() override;

        /*------------------------------------------------------------------------*/
        /** \brief Overide operator==. It is id-based.
         *
         * \param AEdge an edge
         */
        bool operator==(const Edge& AEdge) const;

        /*------------------------------------------------------------------------*/
        /** \brief Overide operator==. It is id-based.
         *
         * \param AEdge an edge
         */
        bool operator!=(const Edge& AEdge) const;


        /*------------------------------------------------------------------------*/
        /** \brief Overide operator<. It is id-based.
         *
         * \param AEdge an edge
         */
        bool operator<(const Edge& AEdge) const {return this->m_id < AEdge.m_id;}

	/*------------------------------------------------------------------------*/
    /** \brief  Accesor to the cell dim.
     */
	int dim() const override {return 1;}

	/*------------------------------------------------------------------------*/
    /** \brief Accessor th the number of incident nodes, edges, faces and
     * 		   adjacent regions
     */
	TInt nbNodes()   const override;
	TInt nbEdges()   const override;
	TInt nbFaces()   const override;
	TInt nbRegions() const override;

	/*------------------------------------------------------------------------*/
    /** \brief Provides the edge length
     */
	TCoord length() const;

		/*------------------------------------------------------------------------*/
		/** \brief Get the opposite node
		 */
		Node getOppositeNode(const Node& ANode) const;
		Node getOppositeNode(const TCellID& ANodeID) const;
	   TCellID getOppositeNodeId(const Node& ANode)const;
	   TCellID getOppositeNodeId(const TCellID& ANodeID)const;
	/*------------------------------------------------------------------------*/
    /** \brief Provides the middle point
     */
	math::Point center() const override;

	/*------------------------------------------------------------------------*/
    /** \brief Provides the edge segment
     */
    math::Segment segment() const;

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

	friend std::ostream & operator << (std::ostream & AStream, const Edge& AN);

private:

	/** A link to the generic edge container of the owner*/
	EdgeContainer* m_edges_container;

};
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_EDGE_H_ */
/*----------------------------------------------------------------------------*/



