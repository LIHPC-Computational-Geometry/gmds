/*----------------------------------------------------------------------------*/
/*
 * EdgeContainer.h
 *
 *  Created on: 19 mai 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_EDGECONTAINER_H_
#define GMDS_EDGECONTAINER_H_
/*----------------------------------------------------------------------------*/
#include <gmds/utils/BitVector.h>
#include <gmds/utils/IndexedVector.h>
/*----------------------------------------------------------------------------*/
#include "Node.h"
#include "Edge.h"
#include "Face.h"
#include "GMDSIg_export.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
    class Mesh;
/*----------------------------------------------------------------------------*/
    class GMDSIg_API EdgeContainer {

        friend class Node;
        friend class Edge;
        friend class Face;
        friend class Region;
        friend class Mesh;
    public:

        explicit EdgeContainer(Mesh* AMesh);
        virtual ~EdgeContainer();

        /*------------------------------------------------------------------------*/
        /** \brief Add an edge in the container
         *  \param AN1 first node id
         *  \param AN2 second node id
         */
        Edge add(const TCellID& AN1,const TCellID& AN2);

        /*------------------------------------------------------------------------*/
        /** \brief Indicate if this container contains a cell of id AID
         *
         *  \param AID a cell id
         */
        bool has(const TCellID& AID) const;
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

        TInt getNbElements() const {return m_edge_ids.size();}
        TCellID getMaxID()   const {return m_edge_ids.top()-1;}

        class iterator {
        public:


            using self_type = iterator;
            using value_type = TCellID;
            using pointer = TCellID*;
            using reference = TCellID&;

            iterator(EdgeContainer* AContainer, bool ABegin):m_it(AContainer->m_edge_ids.begin()) {
                if(!ABegin)
                    m_it=AContainer->m_edge_ids.end();
            }
            iterator(const iterator& AIt) = default;


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
        /** \brief Get the node infos for the edge of id AID
         *
         * \param AID the id of the edge we want to get infos
         * \param ANbNodes the number of nodes of the edge
         */
        void getNodesData(const TCellID& AID, TInt& ANbNodes) const;

        /*------------------------------------------------------------------------*/
        /** \brief Get the edge infos for the edge of id AID
         * \param AID the id of the edge we want to get infos
         * \param ANbEdges the number of adjacent edges
         */
        void getEdgesData(const TCellID& AID, TInt& ANbEdges) const;

        /*------------------------------------------------------------------------*/
        /** \brief Get the face infos for the edge of id AID
         * \param AID the id of the edge we want to get infos
         * \param ANbFaces the number of adjacent faces
         */
        void getFacesData(const TCellID& AID, TInt& ANbFaces) const;


        /*------------------------------------------------------------------------*/
        /** \brief Get the region infos for the edge of id AID
         *
         * \param AID the id of the edge we want to get infos
         * \param ANbRegions the number of regions adjacent to the edge
         */
        void getRegionsData(const TCellID& AID, TInt& ANbRegions) const;

        inline void remove(TInt index)
        {
            m_edge_ids.unselect(index);
            //connectivities container are indexed on m_edge_ids and so not updated
        }
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
        TInt capacity() const {return m_edge_ids.capacity();}

    protected:
        Edge buildEdge(TInt) const;
        /*------------------------------------------------------------------------*/
        /** \brief Add an edge in the container
         *  \param AN1 first node id
         *  \param AN2 second node id
         *  \param AID the id of the edge
         */
        Edge add(const TCellID& AN1,const TCellID& AN2, const TCellID& AID);
    protected:

        Mesh* m_mesh;
        /** supported mesh model */
        MeshModel m_model;

        BitVector m_edge_ids;

        IndexedVector<TabCellID<2> >* m_E2N;
        IndexedVector<TabCellID<size_undef> >* m_E2E;
        IndexedVector<TabCellID<size_undef> >* m_E2F;
        IndexedVector<TabCellID<size_undef> >* m_E2R;

    };
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_EDGECONTAINER_H_ */
/*----------------------------------------------------------------------------*/




