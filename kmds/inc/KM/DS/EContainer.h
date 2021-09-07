/*----------------------------------------------------------------------------*/
/*
 * NodeContainer.h
 *
 *  Created on: 14 feb 2018
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#ifndef KMDS_ECONTAINER_H_
#define KMDS_ECONTAINER_H_
/*----------------------------------------------------------------------------*/
// STL headers
#include <iostream>
/*----------------------------------------------------------------------------*/
// Kokkos headers
#include <Kokkos_Core.hpp>
/*----------------------------------------------------------------------------*/
// KMDS headers
#include <KM/Utils/KTypes.h>
/*----------------------------------------------------------------------------*/
namespace kmds {
/*----------------------------------------------------------------------------*/
    class Mesh;
/*----------------------------------------------------------------------------*/
    class EXPORT_KMDS EContainer
    {
    public:
        /*------------------------------------------------------------------------*/
        /** \brief Constructor
                *
                * \param[in] AOwner the mesh that owns this edge container
                * \param[in] ACapacity initial capacity of the container
                                */
        EContainer(Mesh* AOwner, const TInt32 ACapacity = 16);

        /*------------------------------------------------------------------------*/
        /** \brief  Destructor
         */
        virtual ~EContainer();

        /*------------------------------------------------------------------------*/
        /** \brief Indicate if this container contains a node of id AID
         *
         *  \param AID a cell id
         */
        bool has(const TCellID AID) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Add/allocate \p ANb new cells
         *
         * \param[in] ANb the number of cells to be created
         * \return    The index of the first created cell. Others are contiguous righ behind.
         */
        TCellID add(const TInt32 ANb);

        /*------------------------------------------------------------------------*/
        /** \brief  Add/allocate one
         *
         * \return    The index of the created cell
         */
        TCellID add();
        TCellID add_unsafe();

        TCellID newEdge(const TCellID AN1, const TCellID AN2);
        TCellID newEdge_unsafe(const TCellID AN1, const TCellID AN2);

        /*------------------------------------------------------------------------*/
        /** \brief  Remove/Deallocate cell of id \p AId.
                 *
         * \param[in] AId The index of the cell to be removed.
         */
        void remove(const TCellID AId);

        /*------------------------------------------------------------------------*/
        /** \brief  Set the nodes of cell \p AId
         *
         *  \param[in] AI cell id we want to get the nodes
         *  \param[in] AN new nodes of the cell
         */
        void set(const TCellID AId, const TCellID AN1, const TCellID AN2);

        /*------------------------------------------------------------------------*/
        /** \brief  Get the nodes of cell \p AI
         *
         *  \param[in] AI cell id we want to get the nodes of
         *  \param[in] AN nodes of the cell
         */
        void get(const TCellID AI, Kokkos::View<TCellID*>& AV) const;
        void get(const TCellID AId, TCellID& AN1, TCellID& AN2) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the capacity in terms of array size (including holes)
         */
        TSize capacity() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the top id (including holes)
         */
        TSize top() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Returns the number of cells really stored
         */
        TSize nbCells() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Resize the container
         */
        void resize(const TInt32 ASize);

        /*------------------------------------------------------------------------*/
        /** \brief  Empty the container
         */
        void removeAll();

    private:
        /** mesh owner of this cell container */
        Mesh* m_mesh;

        /** store cell info*/
        Kokkos::View<bool*> m_cells;
        Kokkos::View<TCellID * [2]> m_E2N;
        //, Kokkos::MemoryTraits<Kokkos::Atomic> I remove the memory traits since I cannot resize otherwise

        /** top of the heap including holes so */
        TSize m_top;
    };
/*----------------------------------------------------------------------------*/
}  // namespace KMDS
/*----------------------------------------------------------------------------*/
#endif /* KMDS_ECONTAINER_H_ */
/*----------------------------------------------------------------------------*/
