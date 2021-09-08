/*----------------------------------------------------------------------------*/
/*
 * Region.h
 *
 *  Created on: 19 May 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_REGION_H_
#define GMDS_REGION_H_
/*----------------------------------------------------------------------------*/
#include <set>
/*----------------------------------------------------------------------------*/
#include "Cell.h"
#include "GMDSIg_export.h"
#include <gmds/math/Point.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
class IGMesh;
class EdgeContainer;
class RegionContainer;
class FaceContainer;
/*----------------------------------------------------------------------------*/
/** \class Region
 *
 *  \brief A region instance is an object that provided an object-type access to
 *  	   the data representing a mesh region (3-cell).
 *
 */
class GMDSIg_API Region : public Cell{
public:

	/*------------------------------------------------------------------------*/
	/** \brief Default constructor. Used for stl container initialization
	 */
	Region();

	/*------------------------------------------------------------------------*/
	/** \brief Constructor
	 *
	 * \param AMesh the mesh containing this cell
	 * \param AType the type of cell
	 * \param AID the cell id
	 */
	Region(Mesh* AMesh, const ECellType AType, const TCellID& AID);

	/*------------------------------------------------------------------------*/
	/** \brief Copy Constructor
	 *
	 * \param AReg the region to be copied
	 */
	Region(const Region& AReg);

	/*------------------------------------------------------------------------*/
	/** \brief Destructor
	 */
	virtual ~Region();

        /*------------------------------------------------------------------------*/
        /** \brief Overide operator==. It is id-based.
         *
         * \param ARegion a region
         */
        bool operator==(const Region& ARegion) const;

        /*------------------------------------------------------------------------*/
        /** \brief Overide operator==. It is id-based.
         *
         * \param ARegion a region
         */
        bool operator!=(const Region& ARegion) const;

        /*------------------------------------------------------------------------*/
        /** \brief Overide operator<. It is id-based.
         *
         * \param ARegion a region
         */
        bool operator<(const Region& ARegion) const {return this->m_id < ARegion.m_id;}


	/*------------------------------------------------------------------------*/
    /** \brief Accessor th the number of incident nodes, edges, faces and
     * 		   adjacent regions
     */
	virtual TInt nbNodes()   const;
	virtual TInt nbEdges()   const;
	virtual TInt nbFaces()   const;
	virtual TInt nbRegions() const;

	/*------------------------------------------------------------------------*/
    /** \brief  Accesor to the cell dim.
     */
	virtual int dim() const {return 3;}
    
    /*------------------------------------------------------------------------*/
    /** \brief  Compute the center of the region
     *
     * \return the center of the region
     */
    math::Point center() const;
    
    
    /*------------------------------------------------------------------------*/
    /** \brief  Compute the volume of the region. Only works for tet, hex, 
     *         prism3 and pyramid. The computation is simply based on the
     *         decomposition of the region into tetrahedral elements.
     *
     * \return the center of the region
     */
    TCoord volume() const;
    
        /*------------------------------------------------------------------------*/
        /** \brief  Compute interpolation points of the region. 
 	 *	RIGHT NOW ONLY RETURNS THE NODES OF THE REGION
         *
         * \return the interpolation points
         */
	std::vector<math::Point> computeNGLLPoints(int ADegree) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Compute the scaled jacobian of the region.
         *
         * \return the scaled jacobian
         */
        double computeScaledJacobian() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Compute the scaled jacobian of the region
	 *          normalized between [-1., 1.]
         *
         * \return the scaled jacobian
         */
        double computeNormalizedScaledJacobian() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the mean ratio of the region.
         *
         * \return the mean ratio
         */
        double computeMeanRatio() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Return the faces composing the region in the form of 
 	 * 	vectors of ordered nodes oriented outwards
         *
         * \return the ordered nodes of the faces
         */
        std::vector<std::vector<Node> > getOrderedNodesFaces() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Return the faces composing the region in the form of 
         *      vectors of ordered nodes oriented outwards
         *
         * \return the ordered nodes ids of the faces
         */
        std::vector<std::vector<TCellID> > getOrderedNodesFacesIDs() const;

	/*------------------------------------------------------------------------*/
		/** \brief  Return the edges composing the region in the form of
   	 	 *      vectors of nodes
   	 	 *
     	 * \return the ordered nodes ids of the faces
     	 */
		std::vector<std::vector<TCellID> > getNodesEdgesIDs() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Return whether the face is outward oriented.
         *
         * \param the nodes of the face
         *
         * \return whether the face is outward oriented
         */
	bool isFaceOrientedOutward(std::vector<Node> ANodes) const;

	/*------------------------------------------------------------------------*/
        /** \brief  Return whether the face is outward oriented.
         *
         * \param the ids of the nodes of the face
         *
         * \return whether the face is outward oriented
         */
        bool isFaceOrientedOutward(std::vector<TCellID> AIDs) const;

	/*------------------------------------------------------------------------*/
        /** \brief  Return the fake faces of the region.
         *
         * \return the fake faces
         */
        std::vector<VirtualFace> getFakeFaces() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Return the fake edges of the region
         *
         * \return the fake edges
         */
        std::vector<VirtualEdge> getFakeEdges() const;

    /*------------------------------------------------------------------------*/
    /** \brief
     */
    double findTheNorm(double ATargetVolume, std::set<VirtualFace> freeFaces) const;

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


	friend std::ostream & operator << (std::ostream & AStream, const Region& AN);

private:

	/** A link to the generic region container of the owner*/
	RegionContainer* m_regions_container;

	/** type-id */
	TCellID m_type_id;
};
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_REGION_H_ */
/*----------------------------------------------------------------------------*/





