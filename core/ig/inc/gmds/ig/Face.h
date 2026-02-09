/*----------------------------------------------------------------------------*/
/*
 * Face.h
 *
 *  Created on: 20 févr. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_FACE_H_
#define GMDS_FACE_H_
/*----------------------------------------------------------------------------*/
#include "Cell.h"
#include "GMDSIg_export.h"
#include <gmds/math/Vector.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
	/*----------------------------------------------------------------------------*/
	class Mesh;
	class FaceContainer;
	/*----------------------------------------------------------------------------*/
	/** \class Face
	 *
	 *  \brief A face instance is an object that provided an object-type access to
	 *  	   the data representing a mesh face.
	 *
	 */
	class GMDSIg_API Face : public Cell{
	public:

		/*------------------------------------------------------------------------*/
		/** \brief Default constructor. Used for stl container initialization
		 */
		Face();

		/*------------------------------------------------------------------------*/
		/** \brief Constructor
		 *
		 * \param AMesh the mesh containing this cell
		 * \param AType the type of cell
		 * \param AID the cell id
		 */
		Face(Mesh* AMesh, ECellType AType, const TCellID& AID);

		/*------------------------------------------------------------------------*/
		/** \brief Copy Constructor
		 *
		 * \param AF the face to be copied
		 */
		Face(const Face& AF);


		void operator=(const Face& AFace);

		/*------------------------------------------------------------------------*/
		/** \brief Overide operator==. It is id-based.
		 *
		 * \param AFace a face
		 */
		bool operator==(const Face& AFace) const;

		/*------------------------------------------------------------------------*/
		/** \brief Overide operator==. It is id-based.
		 *
		 * \param AFace a face
		 */
		bool operator!=(const Face& AFace) const;

		/*------------------------------------------------------------------------*/
		/** \brief Overide operator<. It is id-based.
		 *
		 * \param AFace a face
		 */
		bool operator<(const Face& AFace) const { return this->m_id < AFace.m_id; }

		/*------------------------------------------------------------------------*/
		/** \brief Destructor
		 */
		~Face() override;

		inline TInt getTypeID() const { return m_type_id; }
		/*------------------------------------------------------------------------*/
		/** \brief  Accesor to the cell dim.
		 */
		int dim() const override { return 2; }

		/*------------------------------------------------------------------------*/
		/** \brief Accessor th the number of incident nodes, edges, faces and
		 * 		   adjacent regions
		 */
		TInt nbNodes()   const override;
		TInt nbEdges()   const override;
		TInt nbFaces()   const override;
		TInt nbRegions() const override;

		/*------------------------------------------------------------------------*/
                /** \brief Compute the normal of the face.
                 */
		math::Vector3d normal() const;

		/*------------------------------------------------------------------------*/
                /** \brief Compute the normal of the face (triangle that contains ANode).
		 *
		 * \param ANode the node from which the normal will be computed
                 */
                math::Vector3d normal(const Node& ANode) const;

		/*------------------------------------------------------------------------*/
                /** \brief Compute the center of the face
                 */
		math::Point center() const override;

		/*------------------------------------------------------------------------*/
                /** \brief Compute the area of the face
		 * 
		 * \return the area of the face
                 */
		TCoord area() const;

		/*------------------------------------------------------------------------*/
                /** \brief Compute the signed area of the face
		 * 
		 * \return the signed area of the face
                 */
		TCoord signedArea() const;

		/*------------------------------------------------------------------------*/
	        /** \brief  Compute the scaled jacobian of the face
        	 *
	         * \return the scaled jacobian
        	 */
	        double computeScaledJacobian2D() const;

		/*------------------------------------------------------------------------*/
		/** \brief  compute the shortest distance from a point to this.
		 *
		 * 			WARNING: only implemented for triangular faces.
		 *
		 * 			We use the orthogonal projection to the plane containing
		 * 			(*this). If the projected point is in the triangle, we keep
		 * 			it to compute the distance to AP. Otherwise, we project AP
		 * 			onto the closest edge of (*this).
		 *
		 * \param AP a point
		 *
		 * \return the distance
		 */
		TCoord distance(const math::Point& AP) const;
		/*------------------------------------------------------------------------*/
		/** \brief  Orthogonal projection of a point AP onto this.
		 *
		 * 			WARNING: only implemented for triangular faces.
		 *
		 * \param AP a point
		 *
		 * \return the orthogonal projection of AP onto T
		 */
		math::Point project(const math::Point& AP) const;

		/*------------------------------------------------------------------------*/
		/** \brief  Provide the nodes adjacent to ANode1
		 *
		 *  \param ANode1 incoming node
		 *  \param ANode2 first adjacent node in this
		 *  \param ANode3 second adjacent node in this
		 */
		void getAdjacentNodes(const Node& ANode1, Node& ANode2, Node& ANode3);
        void getAdjacentNodes(const TCellID & ANode1, TCellID& ANode2, TCellID& ANode3);

		/*------------------------------------------------------------------------*/
                /** \brief  Provide the edges in an ordered fashion. Same order as the nodes
                 *
                 *  \param AEdges the container of the ordered edges
                 */
                void getOrderedEdges(std::vector<Edge>& AEdges) const;

		/*------------------------------------------------------------------------*/
		/** \brief  Accessor to the incident cells. Only the non-null cells are
		 * 			provided. T can be Node, Edge, Face or Region.
		 */
		void delegateGet(std::vector<Node>&   ACells) const override;
		void delegateGet(std::vector<Edge>&   ACells) const override;
		void delegateGet(std::vector<Face>&   ACells) const override;
		void delegateGet(std::vector<Region>& ACells) const override;

		void delegateGetNodeIDs(std::vector<TCellID>& ACells) const override;
		void delegateGetEdgeIDs(std::vector<TCellID>& ACells) const override;
		void delegateGetFaceIDs(std::vector<TCellID>& ACells) const override;
		void delegateGetRegionIDs(std::vector<TCellID>& ACells) const override;

		void delegateGetAll(std::vector<Node>&   ACells) const override;
		void delegateGetAll(std::vector<Edge>&   ACells) const override;
		void delegateGetAll(std::vector<Face>&   ACells) const override;
		void delegateGetAll(std::vector<Region>& ACells) const override;

		void delegateGetAllNodeIDs  (std::vector<TCellID>& ACells) const override;
        void delegateGetAllEdgeIDs  (std::vector<TCellID>& ACells) const override;
        void delegateGetAllFaceIDs  (std::vector<TCellID>& ACells) const override;
        void delegateGetAllRegionIDs(std::vector<TCellID>& ACells) const override;

		void delegateSetNodeIDs(const std::vector<TCellID>& ACells) override;
		void delegateSetEdgeIDs(const std::vector<TCellID>& ACells) override;
		void delegateSetFaceIDs(const std::vector<TCellID>& ACells) override;
		void delegateSetRegionIDs(const std::vector<TCellID>& ACells) override;

		void delegateNodeAdd(TCellID AElt) override;
		void delegateEdgeAdd(TCellID AElt) override;
		void delegateFaceAdd(TCellID AElt) override;
		void delegateRegionAdd(TCellID AElt) override;

		void delegateNodeRemove(TCellID AElt) override;
		void delegateEdgeRemove(TCellID AElt) override;
		void delegateFaceRemove(TCellID AElt) override;
		void delegateRegionRemove(TCellID AElt) override;

		void delegateNodeReplace(TCellID AID1, TCellID AID2) override;
		void delegateEdgeReplace(TCellID AID1, TCellID AID2) override;
		void delegateFaceReplace(TCellID AID1, TCellID AID2) override;
		void delegateRegionReplace(TCellID AID1, TCellID AID2) override;

		GMDSIg_API friend std::ostream & operator << (std::ostream & AStream, const Face& AN);

	protected:

		/** A link to the generic face container of the owner*/
		FaceContainer* m_faces_container;

		/** type-id */
		TCellID m_type_id;
	};
	/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_FACE_H_ */
/*----------------------------------------------------------------------------*/

