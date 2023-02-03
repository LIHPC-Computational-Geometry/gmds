/*----------------------------------------------------------------------------*/
/*
 * RegionContainer.h
 *
 *  Created on: 20 mai 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_REGIONCONTAINER_H_
#define GMDS_REGIONCONTAINER_H_
/*----------------------------------------------------------------------------*/
#include "Node.h"
#include "Edge.h"
#include "Face.h"
#include "Region.h"
#include "GMDSIg_export.h"
/*----------------------------------------------------------------------------*/
#include <gmds/utils/BitVector.h>
#include <gmds/utils/IndexedVector.h>
#include <gmds/utils/SmartVector.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
	class Mesh;
/*----------------------------------------------------------------------------*/
/** \class RegionContainer
 *
 *  \brief A region container manages the storage of regions for a mesh object
 */
	class GMDSIg_API RegionContainer {

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
		RegionContainer(Mesh* AMesh);

		/*------------------------------------------------------------------------*/
		/** \brief Destructor
         */
		virtual ~RegionContainer();
		/*------------------------------------------------------------------------*/
		/** \brief Update of the mesh model
         *
         *  \param AModel the new model
         */
		inline void setModel(const MeshModel& AModel) {
			m_model = AModel;
		}

		/*------------------------------------------------------------------------*/
		/** \brief Indicate if this container contains a cell of id AID
         *
         *  \param AID a cell id
         */
		bool has(const TCellID& AID) const;
		/*------------------------------------------------------------------------*/
		/** \brief Creation of a tetrahedron
         *
         * \param AN1 Node 1
         * \param AN2 Node 2
         * \param AN3 Node 3
         * \param AN4 Node 4
         * \return a region object that encapsulates access to the mesh region
         */
		Region addTet(const TCellID& AN1,const TCellID& AN2,const TCellID& AN3,
					  const TCellID& AN4);


		/*------------------------------------------------------------------------*/
		/** \brief Creation of a pyramid
         *
         * \param AN1 Node 1
         * \param AN2 Node 2
         * \param AN3 Node 3
         * \param AN4 Node 4
         * \param AN5 Node 5
         * \return a region object that encapsulates access to the mesh region
         */
		Region addPyramid(const TCellID& AN1,const TCellID& AN2,
						  const TCellID& AN3, const TCellID& AN4, const TCellID& AN5);
		/*------------------------------------------------------------------------*/
		/** \brief Creation of a prism3
		 *
		 * \param AN1 Node 1
		 * \param AN2 Node 2
		 * \param AN3 Node 3
		 * \param AN4 Node 4
		 * \param AN5 Node 5
		 * \param AN6 Node 6
		 * \return a region object that encapsulates access to the mesh region
		 */
		Region addPrism3(const TCellID& AN1,const TCellID& AN2,
						 const TCellID& AN3, const TCellID& AN4, const TCellID& AN5,
						 const TCellID& AN6);

		/*------------------------------------------------------------------------*/
		/** \brief Creation of a hexahedral element
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
		Region addHex(const TCellID& AN1, const TCellID& AN2, const TCellID& AN3,
					  const TCellID& AN4, const TCellID& AN5, const TCellID& AN6,
					  const TCellID& AN7, const TCellID& AN8);


		/*------------------------------------------------------------------------*/
		/** \brief Returns the number of elements stored in the container
         */
		TInt getNbElements()const {return m_region_ids.size();}

		/*------------------------------------------------------------------------*/
		/** \brief Returns the max id of a stored element
         */
		TCellID getMaxID()const {return m_region_ids.top()-1;}


		class GMDSIg_API iterator {
		public:


			using self_type = iterator;
			using value_type = TCellID;
			using pointer = TCellID*;
			using reference = TCellID&;

			iterator(RegionContainer* AContainer, bool ABegin):m_it(AContainer->m_region_ids.begin()) {
				if(!ABegin)
					m_it=AContainer->m_region_ids.end();
			}
			iterator(const iterator& AIt):m_it(AIt.m_it){}


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
		iterator begin() {return iterator(this,true);};

		/*------------------------------------------------------------------------*/
		/** \brief Provide an iterator onto the last element of this container
         */
		iterator end() {return iterator(this,false);};

		/*------------------------------------------------------------------------*/
		/** \brief Get the node infos for the region of id AID
         *
         * \param AID the id of the region we want to get infos
         * \param ANbNodes the number of nodes of the face
         */
		void getNodesData(const TCellID& AID, int& ANbNodes) const;

		/*------------------------------------------------------------------------*/
		/** \brief Get the edge infos for the region of id AID
         * \param AID the id of the region we want to get infos
         * \param ANbEdges the number of adjacent edges
         */
		void getEdgesData(const TCellID& AID, int& ANbEdges) const;

		/*------------------------------------------------------------------------*/
		/** \brief Get the face infos for the region of id AID
         * \param AID the id of the region we want to get infos
         * \param ANbFaces the number of adjacent faces
         */
		void getFacesData(const TCellID& AID, int& ANbFaces) const;

		/*------------------------------------------------------------------------*/
		/** \brief Get the region infos for the region of id AID
         * \param AID the id of the region we want to get infos
         * \param ANbRegions the number of regions adjacent to the edge
         */
		void getRegionsData(const TCellID& AID, int& ANbRegions) const;

		/*------------------------------------------------------------------------*/
		/** \brief Removes the face AIndex from the container
         * \param AIndex index of the face to be removed
         */
		void remove(TInt AIndex);
		void clear();
		void resize(const TInt);

		/*------------------------------------------------------------------------*/
		/** \brief  This method is necessary when you want to regularize the
         * 			container after several assignements. Indeed assignement method
         * 			has the advantage to check nothing before insertion. The problem
         * 			is then that the container is no more coherent. This method fix
         * 			the container.
         */
		void update();


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

		/*------------------------------------------------------------------------*/
		/** \brief return the storage capacity
         *
         * \return the storage capacity
         */
		TInt capacity() const {return m_region_ids.capacity();}

	protected:

		/*------------------------------------------------------------------------*/
		/** \brief Provide the type if of a face from its global face id
         */
		inline TCellID  getTypeID(const TCellID& AID) const	{
			return m_region_types[AID].type_id;
		}


		/*------------------------------------------------------------------------*/
		/** \brief Build a region object
         * \param AIndex index of the region to be built
         */
		Region buildRegion(const TInt AIndex) const;

		/*------------------------------------------------------------------------*/
		/** \brief Initialization of the containers depending of the mesh model
         */
		void setConnectivityContainers();

		/*------------------------------------------------------------------------*/
		/** \brief add the containers for accessing to ADim-dimensional cells
         *
         *  \param ADim the dimension of the cells we want to store the adjacency
         *  			to
         */
		void addConnectivityContainers(const TInt ADim);

		/*------------------------------------------------------------------------*/
		/** \brief remove the containers for accessing to ADim-dimensional cells
         *
         *  \param ADim the dimension of the cells we want to suppress the
         *  			adjacency to
         */
		void removeConnectivityContainers(const TInt ADim);

		/*------------------------------------------------------------------------*/
		/** \brief Creation of a tetrahedron
         *
         * \param AN1 Node 1
         * \param AN2 Node 2
         * \param AN3 Node 3
         * \param AN4 Node 4
         * \param AID id of the region
         * \return a region object that encapsulates access to the mesh region
         */
		Region addTet(const TCellID& AN1,const TCellID& AN2,const TCellID& AN3,
					  const TCellID& AN4, const TCellID& AID);


		/*------------------------------------------------------------------------*/
		/** \brief Creation of a pyramid
         *
         * \param AN1 Node 1
         * \param AN2 Node 2
         * \param AN3 Node 3
         * \param AN4 Node 4
         * \param AN5 Node 5
         * \param AID id of the region
         * \return a region object that encapsulates access to the mesh region
         */
		Region addPyramid(const TCellID& AN1,const TCellID& AN2,
						  const TCellID& AN3, const TCellID& AN4, const TCellID& AN5,
						  const TCellID& AID);
		/*------------------------------------------------------------------------*/
		/** \brief Creation of a prism3
         *
         * \param AN1 Node 1
         * \param AN2 Node 2
         * \param AN3 Node 3
         * \param AN4 Node 4
         * \param AN5 Node 5
         * \param AN6 Node 6
         * \param AID id of the region
         * \return a region object that encapsulates access to the mesh region
         */
		Region addPrism3(const TCellID& AN1,const TCellID& AN2,
						 const TCellID& AN3, const TCellID& AN4, const TCellID& AN5,
						 const TCellID& AN6, const TCellID& AID);

		/*------------------------------------------------------------------------*/
		/** \brief Creation of a hexahedral element
         *
         * \param AN1 Node 1
         * \param AN2 Node 2
         * \param AN3 Node 3
         * \param AN4 Node 4
         * \param AN5 Node 5
         * \param AN6 Node 6
         * \param AN7 Node 7
         * \param AN8 Node 8
         * \param AID id of the region
         * \return a region object that encapsulates access to the mesh region
         */
		Region addHex(const TCellID& AN1, const TCellID& AN2, const TCellID& AN3,
					  const TCellID& AN4, const TCellID& AN5, const TCellID& AN6,
					  const TCellID& AN7, const TCellID& AN8, const TCellID& AID);




	protected:

		Mesh* m_mesh;
		/** supported mesh model */
		MeshModel m_model;

		/** bit set container to access to regions */
		BitVector m_region_ids;

		/** \struct RegionInfo
         * \brief Nested structure to handle some region informations like cell
         * 		  type, and the face id for this type
         */
		struct RegionInfo{
			ECellType type;    ///region type : tet, hex, prism, pyramid,...
			TInt	  type_id; ///id of the typed face
			RegionInfo(ECellType t=GMDS_TETRA, TInt i=1):type(t), type_id(i){}
		};

		/** Indexed collection of region types*/
		IndexedVector<RegionInfo> m_region_types;
#ifdef GMDS_PARALLEL
		IndexedVector<DistributedCellData> m_distributed_data;
#endif //GMDS_PARALLEL

		/** containers of connectivity depending of the cell type */
		SmartVector<TabCellID<4> >* m_T2N;
		SmartVector<TabCellID<6> >* m_T2E;
		SmartVector<TabCellID<4> >* m_T2F;
		SmartVector<TabCellID<4> >* m_T2R;

		SmartVector<TabCellID<5> >* m_PY2N;
		SmartVector<TabCellID<8> >* m_PY2E;
		SmartVector<TabCellID<5> >* m_PY2F;
		SmartVector<TabCellID<5> >* m_PY2R;

		SmartVector<TabCellID<6> >* m_PR2N;
		SmartVector<TabCellID<9> >* m_PR2E;
		SmartVector<TabCellID<5> >* m_PR2F;
		SmartVector<TabCellID<5> >* m_PR2R;

		SmartVector<TabCellID<8> >*  m_H2N;
		SmartVector<TabCellID<12> >* m_H2E;
		SmartVector<TabCellID<6> >*  m_H2F;
		SmartVector<TabCellID<6> >*  m_H2R;

		SmartVector<TabCellID<size_undef> >* m_P2N;
		SmartVector<TabCellID<size_undef> >* m_P2E;
		SmartVector<TabCellID<size_undef> >* m_P2F;
		SmartVector<TabCellID<size_undef> >* m_P2R;


		/** \struct Nested structure to handle some region informations
         */

		template<int N> struct AdjUpdate
		{
			SmartVector<TabCellID<N> >* m_adj;

			AdjUpdate(SmartVector<TabCellID<N> >* adj):m_adj(adj){}

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
         * 		  adjacency relations for tetrahedra.
         */
		struct TAccessor{
			RegionContainer* m_owner;
			SmartVector<TabCellID<4> >* m_N;
			SmartVector<TabCellID<6> >* m_E;
			SmartVector<TabCellID<4> >* m_F;
			SmartVector<TabCellID<4> >* m_R;
			AdjUpdate<4>* adj_N;
			AdjUpdate<6>* adj_E;
			AdjUpdate<4>* adj_F;
			AdjUpdate<4>* adj_R;

			TAccessor(RegionContainer* AOwner, const MeshModel& AModel);
			~TAccessor();
			TInt getID();
		};

		/** \struct PyAccessor
         * \param instanciate generic pointers to specialize the access to cells and
         * 		  adjacency relations for pyramids.
         */
		struct PyAccessor{
			RegionContainer* m_owner;
			SmartVector<TabCellID<5> >* m_N;
			SmartVector<TabCellID<8> >* m_E;
			SmartVector<TabCellID<5> >* m_F;
			SmartVector<TabCellID<5> >* m_R;
			AdjUpdate<5>* adj_N;
			AdjUpdate<8>* adj_E;
			AdjUpdate<5>* adj_F;
			AdjUpdate<5>* adj_R;

			PyAccessor(RegionContainer* AOwner, const MeshModel& AModel);
			~PyAccessor();
			TInt getID();
		};
		/** \struct PrAccessor
         * \param instanciate generic pointers to specialize the access to cells and
         * 		  adjacency relations for prisms 3.
         */
		struct PrAccessor{
			RegionContainer* m_owner;
			SmartVector<TabCellID<6> >* m_N;
			SmartVector<TabCellID<9> >* m_E;
			SmartVector<TabCellID<5> >* m_F;
			SmartVector<TabCellID<5> >* m_R;
			AdjUpdate<6>* adj_N;
			AdjUpdate<9>* adj_E;
			AdjUpdate<5>* adj_F;
			AdjUpdate<5>* adj_R;

			PrAccessor(RegionContainer* AOwner, const MeshModel& AModel);
			~PrAccessor();
			TInt getID();
		};
		/** \struct HAccessor
         * \param instanciate generic pointers to specialize the access to cells and
         * 		  adjacency relations for hexahedra.
         */
		struct HAccessor{
			RegionContainer* m_owner;
			SmartVector<TabCellID<8 > >* m_N;
			SmartVector<TabCellID<12> >* m_E;
			SmartVector<TabCellID<6 > >* m_F;
			SmartVector<TabCellID<6 > >* m_R;
			AdjUpdate<8>* adj_N;
			AdjUpdate<12>* adj_E;
			AdjUpdate<6>* adj_F;
			AdjUpdate<6>* adj_R;

			HAccessor(RegionContainer* AOwner, const MeshModel& AModel);
			~HAccessor();
			TInt getID();
		};
		/** \struct PAccessor
         * \param instanciate generic pointers to specialize the access to cells and
         * 		  adjacency relations for polyhedra.
         */
		struct PAccessor{
			RegionContainer* m_owner;
			SmartVector<TabCellID<size_undef> >* m_N;
			SmartVector<TabCellID<size_undef> >* m_E;
			SmartVector<TabCellID<size_undef> >* m_F;
			SmartVector<TabCellID<size_undef> >* m_R;
			AdjUpdate<size_undef>* adj_N;
			AdjUpdate<size_undef>* adj_E;
			AdjUpdate<size_undef>* adj_F;
			AdjUpdate<size_undef>* adj_R;

			PAccessor(RegionContainer* AOwner, const MeshModel& AModel);
			~PAccessor();
			TInt getID();
		};

		/** accessor to tetrahedral elements */
		TAccessor* m_tet;
		/** accessor to pyramid elements */
		PyAccessor* m_pyra;
		/** accessor to prism3 elements */
		PrAccessor* m_prism3;
		/** accessor to hexahedral elements */
		HAccessor* m_hex;
		/** accessor to polyhedral elements */
		PAccessor* m_poly;
	};
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_REGIONCONTAINER_H_ */
/*----------------------------------------------------------------------------*/
