/*----------------------------------------------------------------------------*/
/*
 * NodeContainer.h
 *
 *  Created on: 26 mars 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef KMDS_NCONTAINER_H_
#define KMDS_NCONTAINER_H_
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
class EXPORT_KMDS NContainer
{
 public:
        /*------------------------------------------------------------------------*/
        /** \brief Constructor
                *
                * \param[in] AOwner the mesh that owns this node container
                * \param[in] ACapacity initial capacity of the container
                                */
        NContainer(Mesh* AOwner, const TInt32 ACapacity = 16);

        /*------------------------------------------------------------------------*/
        /** \brief  Destructor
         */
        virtual ~NContainer();

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

        TCellID newNode(const TCoord AX, const TCoord AY, const TCoord AZ);

        /*------------------------------------------------------------------------*/
        /** \brief  Remove/Desallocate cell of id \p AId.
                 *
         * \param[in] AId The index of the cell to be removed.
         */
        void remove(const TCellID AId);

        /*------------------------------------------------------------------------*/
        /** \brief  Set the coordinates of node \p AId.
         * \param[in] AId The index of the cell to be removed.
         * \param[in] AX  X coordinate
         * \param[in] AY  Y coordinate
         * \param[in] AZ  Z coordinate */
        void set(const TCellID AId, const TCoord AX, const TCoord AY, const TCoord AZ);
        /*------------------------------------------------------------------------*/
        /** \brief  Get the coordinate of node of id \p AId.
         * \param[in] AId The index of the cell to be removed.
         * \param[out] AX  X coordinate
         * \param[out] AY  Y coordinate
         * \param[out] AZ  Z coordinate */
        void get(const TCellID AId, TCoord& AX, TCoord& AY, TCoord& AZ) const;

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
        Kokkos::View<double * [3]> m_xyz;
        //, Kokkos::MemoryTraits<Kokkos::Atomic> I remove the memory traits since I cannot resize otherwise

        /** top of the heap including holes so */
        TSize m_top;
};
/*----------------------------------------------------------------------------*/
}  // namespace KMDS
/*----------------------------------------------------------------------------*/
#endif /* KMDS_NCONTAINER_H_ */
/*----------------------------------------------------------------------------*/
