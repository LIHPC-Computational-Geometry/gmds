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
	virtual ~Edge();

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
	virtual int dim() const {return 1;}

	/*------------------------------------------------------------------------*/
    /** \brief Accessor th the number of incident nodes, edges, faces and
     * 		   adjacent regions
     */
	virtual TInt nbNodes()   const;
	virtual TInt nbEdges()   const;
	virtual TInt nbFaces()   const;
	virtual TInt nbRegions() const;

	/*------------------------------------------------------------------------*/
    /** \brief Provides the edge length
     */
	TCoord length() const;

	/*------------------------------------------------------------------------*/
    /** \brief Provides the middle point
     */
	math::Point center() const;

	/*------------------------------------------------------------------------*/
    /** \brief Provides the edge segment
     */
    math::Segment segment() const;

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



