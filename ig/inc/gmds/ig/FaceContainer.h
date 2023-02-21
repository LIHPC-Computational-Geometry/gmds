/*----------------------------------------------------------------------------*/
/*
 * FaceContainer.h
 *
 *  Created on: 26 mars 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_FACECONTAINER_H_
#define GMDS_FACECONTAINER_H_
/*----------------------------------------------------------------------------*/
#include "Node.h"
#include "Face.h"
#include "GMDSIg_export.h"
#include <gmds/utils/BitVector.h>
#include <gmds/utils/IndexedVector.h>
#include <gmds/utils/SmartVector.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
	class Mesh;
/*----------------------------------------------------------------------------*/
/** \class FaceContainer
 *
 *  \brief A face container manages the storage of faces for a mesh object
 */
	class GMDSIg_API FaceContainer {

	public:

		friend class Node;
		friend class Edge;
		friend class Face;
		friend class Region;
		friend class Mesh;

		/*------------------------------------------------------------------------*/
		/** \brief Constructor.
         *
         * \param AMesh the mesh this face container builds faces for
         */
		explicit FaceContainer(Mesh* AMesh);

		/*------------------------------------------------------------------------*/
		/** \brief Destructor
         */
		virtual ~FaceContainer();

		/*------------------------------------------------------------------------*/
		/** \brief Update of the mesh model
         *
         *  \param AModel the new model
         */
		inline void setModel(const MeshModel& AModel) {
			m_model = AModel;
		}
		/*------------------------------------------------------------------------*/
		/** \brief add the containers for accessing to ADim-dimensional cells
         *
         *  \param ADim the dimension of the cells we want to store the adjacency
         *  			to
         */
		void addConnectivityContainers(TInt ADim);

		/*------------------------------------------------------------------------*/
		/** \brief remove the containers for accessing to ADim-dimensional cells
         *
         *  \param ADim the dimension of the cells we want to suppress the
         *  			adjacency to
         */
		void removeConnectivityContainers(TInt ADim);
		/*------------------------------------------------------------------------*/
		/** \brief Creation of a triangle
         *
         * \param AN1 Node 1
         * \param AN2 Node 2
         * \param AN3 Node 3
         * \return a face object that encapsulates access to the mesh face
         */
		Face addTriangle(const TCellID& AN1,const TCellID& AN2,const TCellID& AN3);
		Face addTriangle();

		/*------------------------------------------------------------------------*/
		/** \brief Creation of a quad
         *
         * \param AN1 Node 1
         * \param AN2 Node 2
         * \param AN3 Node 3
         * \param AN4 Node 4
         * \return a face object that encapsulates access to the mesh face
         */
		Face addQuad(const TCellID& AN1,const TCellID& AN2,const TCellID& AN3, const TCellID& AN4);

		/*------------------------------------------------------------------------*/
		/** \brief Creation of a polygon
         *
         * \param ANodes ordered vector of nodes
         * \return a face object that encapsulates access to the mesh face
         */
		Face addPolygon(const std::vector<TCellID>& ANodes);

		/*------------------------------------------------------------------------*/
		/** \brief Returns the number of elements stored in the container
         */
		TInt getNbElements()const {return m_face_ids.size();}

		/*------------------------------------------------------------------------*/
		/** \brief Returns the max id of a stored element
         */
		TCellID getMaxID()const {return m_face_ids.top()-1;}


		class GMDSIg_API iterator {
		public:


			using self_type = iterator;
			using value_type = TCellID;
			using pointer = TCellID*;
			using reference = TCellID&;

			iterator(FaceContainer* AContainer, bool ABegin):m_it(AContainer->m_face_ids.begin()) {
				if(!ABegin)
					m_it=AContainer->m_face_ids.end();
			}
			iterator(const iterator& AIt)= default;


			self_type operator++() {
				++m_it;
                return *this;
			}

			reference operator*()  { return *m_it; }
			pointer operator->() { return &(*m_it); }
			bool operator==(const self_type& rhs) { return (m_it==rhs.m_it);}
			bool operator!=(const self_type& rhs) { return (m_it!=rhs.m_it);}
		private:
			BitVector::iterator m_it;
		};


		/*------------------------------------------------------------------------*/
		/** \brief Provide an iterator onto the first element of this container
         */
		iterator begin() {return {this,true};};

		/*------------------------------------------------------------------------*/
		/** \brief Provide an iterator onto the last element of this container
         */
		iterator end() {return {this,false};};

		/*------------------------------------------------------------------------*/
		/** \brief Get the node infos for the face of id AID
         *
         * \param AID the id of the face we want to get infos
         * \param ANbNodes the number of nodes of the face
         */
		void getNodesData(const TCellID& AID, TInt& ANbNodes) const;

		/*------------------------------------------------------------------------*/
		/** \brief Get the edge infos for the face of id AID
         * \param AID the id of the face we want to get infos
         * \param ANbEdges the number of adjacent edges
         */
		void getEdgesData(const TCellID& AID, TInt& ANbEdges) const;

		/*------------------------------------------------------------------------*/
		/** \brief Get the face infos for the face of id AID
         * \param AID the id of the face we want to get infos
         * \param ANbFaces the number of adjacent faces
         */
		void getFacesData(const TCellID& AID, TInt& ANbFaces) const;

		/*------------------------------------------------------------------------*/
		/** \brief Get the region infos for the face of id AID
         *
         * \param AID the id of the face we want to get infos
         * \param ANbRegions the number of regions adjacent to the face
         */
		void getRegionsData(const TCellID& AID, TInt& ANbRegions) const;

		/*------------------------------------------------------------------------*/
		/** \brief Removes the face AIndex from the container
         * \param AIndex index of the face to be removed
         */
		void remove(TInt AIndex);

		void clear();
		void resize(TInt);

		/*------------------------------------------------------------------------*/
		/** \brief  This method is necessary when you want to regularize the
         * 			container after several assignements. Indeed assignement method
         * 			has the advantage to check nothing before insertion. The problem
         * 			is then that the container is no more coherent. This method fix
         * 			the container.
         */
		void update();

		/*------------------------------------------------------------------------*/
		/** \brief Indicate if this container contains a cell of id AID
         *
         *  \param AID a cell id
         */
		bool has(const TCellID& AID) const;



		/*------------------------------------------------------------------------*/
		/** \brief serialize (*this) in AStr
         *
         * \param AStr an output streammap
         */
		void serialize(std::ostream& AStr);

		/*------------------------------------------------------------------------*/
		/** \brief unserialize (*this) from AStr
         *
         * \param AStr an input stream
         */
		void unserialize(std::istream& AStr);

	protected:


		/*------------------------------------------------------------------------*/
		/** \brief Creation of a triangle
         *
         * \param AN1 Node 1
         * \param AN2 Node 2
         * \param AN3 Node 3
         * \param AID the id of the face
         * \return a face object that encapsulates access to the mesh face
         */
		Face addTriangle(const TCellID& AN1,const TCellID& AN2,
						 const TCellID& AN3, const TCellID& AID);

		/*------------------------------------------------------------------------*/
		/** \brief Creation of a quad
         *
         * \param AN1 Node 1
         * \param AN2 Node 2
         * \param AN3 Node 3
         * \param AN4 Node 4
         * \param AID the id of the face
         * \return a face object that encapsulates access to the mesh face
         */
		Face addQuad(const TCellID& AN1,const TCellID& AN2,const TCellID& AN3,
					 const TCellID& AN4, const TCellID& AID);

		/*------------------------------------------------------------------------*/
		/** \brief Creation of a polygon
         *
         * \param ANodes ordered vector of nodes
         * \param AID the id of the face
         * \return a face object that encapsulates access to the mesh face
         */
		Face addPolygon(const std::vector<TCellID>& ANodes, const TCellID& AID);

		/*------------------------------------------------------------------------*/
		/** \brief Provide the type if of a face from its global face id
         */
		inline TCellID  getTypeID(const TCellID& AID) const	{
			return m_face_types[AID].type_id;
		}


		/*------------------------------------------------------------------------*/
		/** \brief Build a face object
         * \param AIndex index of the face to be built
         */
		Face buildFace(TInt AIndex) const;

		/*------------------------------------------------------------------------*/
		/** \brief Initialization of the containers depending of the mesh model
         */
		void setConnectivityContainers();

		/*------------------------------------------------------------------------*/
		/** \brief return the storage capacity
         *
         * \return the storage capacity
         */
		TInt capacity() const {return m_face_ids.capacity();}

	protected:

		Mesh* m_mesh;
		/** supported mesh model */
		MeshModel m_model;

		/** bit set container to access to faces */
		BitVector m_face_ids;

		/** \struct FaceInfo
         * \brief Nested structure to handle some face informations like cell
         * 		  type, and the face id for this type
         */
		struct FaceInfo{
			ECellType type;    ///face type : triangle, quad, polygon
			TInt	  type_id; ///id of the typed face
			explicit FaceInfo(ECellType t=GMDS_TRIANGLE, TInt i=1):type(t), type_id(i){}
		};

		/** Indexed collection of face types*/
		IndexedVector<FaceInfo> m_face_types;
#ifdef GMDS_PARALLEL
		IndexedVector<DistributedCellData> m_distributed_data;
#endif //GMDS_PARALLEL

		/** containers of connectivity depending of the cell type */
		SmartVector<TabCellID<3> >* m_T2N;
		SmartVector<TabCellID<3> >* m_T2E;
		SmartVector<TabCellID<3> >* m_T2F;
		SmartVector<TabCellID<2> >* m_T2R;

		SmartVector<TabCellID<4> >* m_Q2N;
		SmartVector<TabCellID<4> >* m_Q2E;
		SmartVector<TabCellID<4> >* m_Q2F;
		SmartVector<TabCellID<2> >* m_Q2R;

		SmartVector<TabCellID<size_undef> >* m_P2N;
		SmartVector<TabCellID<size_undef> >* m_P2E;
		SmartVector<TabCellID<size_undef> >* m_P2F;
		SmartVector<TabCellID<2> >* 		 m_P2R;


		/** \struct Nested structure to handle some face informations
         */

		template<int N> struct AdjUpdate
		{
			SmartVector<TabCellID<N> >* m_adj;

			explicit AdjUpdate(SmartVector<TabCellID<N> >* adj):m_adj(adj){}

			size_t select(){
				return m_adj->selectNewIndex();
			}
			size_t getSizeOfAnElement(){
				return N;
			}
			// TCellID* getElements(TCellID type_id){return &((*m_adj)[type_id]).val[0];}
		};


		/** \struct TAccessor
         * \param instanciate generic pointers to specialize the access to cells and
         * 		  adjacency relations for triangles.
         */
		struct TAccessor{
		public:
			FaceContainer* m_owner;
			SmartVector<TabCellID<3> >* m_N;
			SmartVector<TabCellID<3> >* m_E;
			SmartVector<TabCellID<3> >* m_F;
			SmartVector<TabCellID<2> >* m_R;
			AdjUpdate<3>* adj_N;
			AdjUpdate<3>* adj_E;
			AdjUpdate<3>* adj_F;
			AdjUpdate<2>* adj_R;

			TAccessor(FaceContainer* AOwner, const MeshModel& AModel);
			virtual ~TAccessor();
			TInt getID() const;
		};
		/** \struct QAccessor
         * \param instanciate generic pointers to specialize the access to cells and
         * 		  adjacency relations for quads.
         */
		struct QAccessor{
		public:
			FaceContainer* m_owner;
			SmartVector<TabCellID<4> >* m_N;
			SmartVector<TabCellID<4> >* m_E;
			SmartVector<TabCellID<4> >* m_F;
			SmartVector<TabCellID<2> >* m_R;
			AdjUpdate<4>* adj_N;
			AdjUpdate<4>* adj_E;
			AdjUpdate<4>* adj_F;
			AdjUpdate<2>* adj_R;

			QAccessor(FaceContainer* AOwner, const MeshModel& AModel);
			virtual ~QAccessor();
			TInt getID() const;
		};
		/** \struct PAccessor
         * \param instanciate generic pointers to specialize the access to cells and
         * 		  adjacency relations for polygons.
         */
		struct PAccessor{
		public:
			FaceContainer* m_owner;
			SmartVector<TabCellID<size_undef> >* m_N;
			SmartVector<TabCellID<size_undef> >* m_E;
			SmartVector<TabCellID<size_undef> >* m_F;
			SmartVector<TabCellID<2> >* m_R;
			AdjUpdate<size_undef>* adj_N;
			AdjUpdate<size_undef>* adj_E;
			AdjUpdate<size_undef>* adj_F;
			AdjUpdate<2>* adj_R;

			PAccessor(FaceContainer* AOwner, const MeshModel& AModel);
			virtual ~PAccessor();
			TInt getID() const;
		};

		/** accessor to triangles */
		TAccessor* m_triangles;
		/** accessor to quads */
		QAccessor* m_quads;
		/** accessor to polygons */
		PAccessor* m_polygons;
	};
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_FACECONTAINER_H_ */
/*----------------------------------------------------------------------------*/